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


#include "SimdElectronRepulsionVrrRecGK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_gk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

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
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
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
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
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
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
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
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_41 = buffer.data(gh0 + 41);
    const auto *gh0_42 = buffer.data(gh0 + 42);
    const auto *gh0_43 = buffer.data(gh0 + 43);
    const auto *gh0_44 = buffer.data(gh0 + 44);
    const auto *gh0_45 = buffer.data(gh0 + 45);
    const auto *gh0_46 = buffer.data(gh0 + 46);
    const auto *gh0_47 = buffer.data(gh0 + 47);
    const auto *gh0_48 = buffer.data(gh0 + 48);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_74 = buffer.data(gh0 + 74);
    const auto *gh0_75 = buffer.data(gh0 + 75);
    const auto *gh0_76 = buffer.data(gh0 + 76);
    const auto *gh0_77 = buffer.data(gh0 + 77);

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
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_42 = buffer.data(gh1 + 42);
    const auto *gh1_43 = buffer.data(gh1 + 43);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);
    const auto *gh1_48 = buffer.data(gh1 + 48);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_62 = buffer.data(gh1 + 62);
    const auto *gh1_63 = buffer.data(gh1 + 63);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_74 = buffer.data(gh1 + 74);
    const auto *gh1_75 = buffer.data(gh1 + 75);
    const auto *gh1_76 = buffer.data(gh1 + 76);
    const auto *gh1_77 = buffer.data(gh1 + 77);

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
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, \
                         gi_0, gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = pb_y[k] * gi_2[k];

        t_5[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, \
                         gh1_2, gh1_3, gi_3, gi_4, gi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_7[k] = pb_z[k] * gi_3[k];

        t_8[k] = pb_y[k] * gi_4[k];

        t_9[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_10[k] = f_7 * gh0_3[k]
                  - f_8 * gh1_3[k]
                  + pb_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, gh0_4, gh0_5, gh1_4, \
                         gh1_5, gi_5, gi_6, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gi_5[k];

        t_12[k] = f_3 * gh0_4[k]
                  - f_4 * gh1_4[k]
                  + pb_y[k] * gi_6[k];

        t_13[k] = pb_y[k] * gi_7[k];

        t_14[k] = f_7 * gh0_4[k]
                  - f_8 * gh1_4[k]
                  + pb_z[k] * gi_7[k];

        t_15[k] = f_9 * gh0_5[k]
                  - f_10 * gh1_5[k]
                  + pb_y[k] * gi_8[k];

        t_16[k] = pb_z[k] * gi_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, gh0_6, gh0_7, gh1_6, gh1_7, gi_9, \
                         gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * gh0_6[k]
                  - f_6 * gh1_6[k]
                  + pb_y[k] * gi_9[k];

        t_18[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_19[k] = pb_y[k] * gi_11[k];

        t_20[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, fi_14, fi_16, fi_17, fi_18, \
                         gi_12, gi_14, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * fi_14[k]
                  + pb_x[k] * gi_14[k];

        t_22[k] = pb_z[k] * gi_12[k];

        t_23[k] = f_0 * fi_16[k]
                  + pb_x[k] * gi_15[k];

        t_24[k] = f_0 * fi_17[k]
                  + pb_x[k] * gi_16[k];

        t_25[k] = f_0 * fi_18[k]
                  + pb_x[k] * gi_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, fi_20, gh0_8, gh1_8, gi_13, \
                         gi_14, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * gi_13[k];

        t_27[k] = f_0 * fi_20[k]
                  + pb_x[k] * gi_19[k];

        t_28[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_14[k];

        t_29[k] = pb_z[k] * gi_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, gh0_9, gh0_10, gh0_11, gh1_9, gh1_10, gh1_11, \
                         gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_9[k]
                  + pb_y[k] * gi_15[k];

        t_31[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_10[k]
                  + pb_y[k] * gi_16[k];

        t_32[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_11[k]
                  + pb_y[k] * gi_17[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, fi_0, fk_0, \
                         gh0_12, gh1_12, gi_18, gi_19, gi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_12[k]
                  + pb_y[k] * gi_18[k];

        t_34[k] = pb_y[k] * gi_19[k];

        t_35[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_12[k]
                  + pb_z[k] * gi_19[k];

        t_36[k] = pa_y[k] * fk_0[k];

        t_37[k] = f_11 * fi_0[k]
                  + pb_y[k] * gi_20[k];

        t_38[k] = pb_z[k] * gi_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, fi_1, fi_3, fk_1, fk_2, \
                         fk_3, gi_21, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * fi_1[k]
                  + pa_y[k] * fk_1[k];

        t_40[k] = pb_z[k] * gi_21[k];

        t_41[k] = pa_y[k] * fk_2[k];

        t_42[k] = f_13 * fi_3[k]
                  + pa_y[k] * fk_3[k];

        t_43[k] = pb_z[k] * gi_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, fi_4, fi_5, fi_7, \
                         fk_4, fk_5, fk_6, gi_23, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * fi_4[k]
                  + pb_y[k] * gi_23[k];

        t_45[k] = pa_y[k] * fk_4[k];

        t_46[k] = f_0 * fi_5[k]
                  + pa_y[k] * fk_5[k];

        t_47[k] = pb_z[k] * gi_24[k];

        t_48[k] = f_12 * fi_7[k]
                  + pa_y[k] * fk_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, fi_8, fi_9, fi_11, \
                         fk_7, fk_8, fk_9, gi_25, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * fi_8[k]
                  + pb_y[k] * gi_25[k];

        t_50[k] = pa_y[k] * fk_7[k];

        t_51[k] = f_14 * fi_9[k]
                  + pa_y[k] * fk_8[k];

        t_52[k] = pb_z[k] * gi_26[k];

        t_53[k] = f_13 * fi_11[k]
                  + pa_y[k] * fk_9[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, fi_12, fi_13, fi_28, fk_10, \
                         fk_11, gi_27, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * fi_12[k]
                  + pa_y[k] * fk_10[k];

        t_55[k] = f_11 * fi_13[k]
                  + pb_y[k] * gi_27[k];

        t_56[k] = pa_y[k] * fk_11[k];

        t_57[k] = f_13 * fi_28[k]
                  + pb_x[k] * gi_29[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, fi_29, fi_30, fi_31, fi_32, \
                         gi_28, gi_30, gi_31, gi_32, gi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * gi_28[k];

        t_59[k] = f_13 * fi_29[k]
                  + pb_x[k] * gi_30[k];

        t_60[k] = f_13 * fi_30[k]
                  + pb_x[k] * gi_31[k];

        t_61[k] = f_13 * fi_31[k]
                  + pb_x[k] * gi_32[k];

        t_62[k] = f_13 * fi_32[k]
                  + pb_x[k] * gi_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, fi_14, fi_16, fi_17, fk_13, \
                         fk_14, fk_15, fk_16, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * fk_13[k];

        t_64[k] = f_15 * fi_14[k]
                  + pa_y[k] * fk_14[k];

        t_65[k] = pb_z[k] * gi_29[k];

        t_66[k] = f_14 * fi_16[k]
                  + pa_y[k] * fk_15[k];

        t_67[k] = f_0 * fi_17[k]
                  + pa_y[k] * fk_16[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, fi_18, fi_19, fi_20, \
                         fk_0, fk_17, fk_18, fk_19, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_13 * fi_18[k]
                  + pa_y[k] * fk_17[k];

        t_69[k] = f_12 * fi_19[k]
                  + pa_y[k] * fk_18[k];

        t_70[k] = f_11 * fi_20[k]
                  + pb_y[k] * gi_34[k];

        t_71[k] = pa_y[k] * fk_19[k];

        t_72[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, fi_0, fi_2, \
                         fk_1, fk_2, fk_3, gi_35, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * gi_35[k];

        t_74[k] = f_11 * fi_0[k]
                  + pb_z[k] * gi_35[k];

        t_75[k] = pa_z[k] * fk_1[k];

        t_76[k] = pb_y[k] * gi_36[k];

        t_77[k] = f_12 * fi_2[k]
                  + pa_z[k] * fk_2[k];

        t_78[k] = pa_z[k] * fk_3[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, fi_3, fi_4, fi_5, \
                         fk_4, fk_5, gi_37, gi_38, gi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * fi_3[k]
                  + pb_z[k] * gi_37[k];

        t_80[k] = pb_y[k] * gi_38[k];

        t_81[k] = f_13 * fi_4[k]
                  + pa_z[k] * fk_4[k];

        t_82[k] = pa_z[k] * fk_5[k];

        t_83[k] = f_11 * fi_5[k]
                  + pb_z[k] * gi_39[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, fi_6, fi_8, fi_9, \
                         fk_6, fk_7, fk_8, gi_40, gi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * fi_6[k]
                  + pa_z[k] * fk_6[k];

        t_85[k] = pb_y[k] * gi_40[k];

        t_86[k] = f_0 * fi_8[k]
                  + pa_z[k] * fk_7[k];

        t_87[k] = pa_z[k] * fk_8[k];

        t_88[k] = f_11 * fi_9[k]
                  + pb_z[k] * gi_41[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, fi_10, fi_11, fi_13, fk_9, \
                         fk_10, fk_11, fk_12, gi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * fi_10[k]
                  + pa_z[k] * fk_9[k];

        t_90[k] = f_13 * fi_11[k]
                  + pa_z[k] * fk_10[k];

        t_91[k] = pb_y[k] * gi_42[k];

        t_92[k] = f_14 * fi_13[k]
                  + pa_z[k] * fk_11[k];

        t_93[k] = pa_z[k] * fk_12[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, fi_46, fi_47, fi_48, fi_49, \
                         gi_43, gi_45, gi_46, gi_47, gi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_13 * fi_46[k]
                  + pb_x[k] * gi_45[k];

        t_95[k] = f_13 * fi_47[k]
                  + pb_x[k] * gi_46[k];

        t_96[k] = f_13 * fi_48[k]
                  + pb_x[k] * gi_47[k];

        t_97[k] = f_13 * fi_49[k]
                  + pb_x[k] * gi_48[k];

        t_98[k] = pb_y[k] * gi_43[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, fi_14, fi_15, fi_51, \
                         fk_14, fk_15, gi_44, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * fi_51[k]
                  + pb_x[k] * gi_49[k];

        t_100[k] = pa_z[k] * fk_14[k];

        t_101[k] = f_11 * fi_14[k]
                   + pb_z[k] * gi_44[k];

        t_102[k] = f_12 * fi_15[k]
                   + pa_z[k] * fk_15[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, fi_16, fi_17, fi_18, \
                         fi_20, fk_16, fk_17, fk_18, fk_19, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_13 * fi_16[k]
                   + pa_z[k] * fk_16[k];

        t_104[k] = f_0 * fi_17[k]
                   + pa_z[k] * fk_17[k];

        t_105[k] = f_14 * fi_18[k]
                   + pa_z[k] * fk_18[k];

        t_106[k] = pb_y[k] * gi_49[k];

        t_107[k] = f_15 * fi_20[k]
                   + pa_z[k] * fk_19[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_y, pb_z, dk0_0, dk1_0, fi_21, fk_20, \
                         gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_16 * dk0_0[k]
                   - f_17 * dk1_0[k]
                   + pa_y[k] * fk_20[k];

        t_109[k] = f_12 * fi_21[k]
                   + pb_y[k] * gi_50[k];

        t_110[k] = pb_z[k] * gi_50[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pb_z, fi_53, gh0_13, gh0_15, gh1_13, \
                         gh1_15, gi_51, gi_52, gi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * fi_53[k]
                   + f_9 * gh0_15[k]
                   - f_10 * gh1_15[k]
                   + pb_x[k] * gi_53[k];

        t_112[k] = pb_z[k] * gi_51[k];

        t_113[k] = f_3 * gh0_13[k]
                   - f_4 * gh1_13[k]
                   + pb_z[k] * gi_52[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pb_y, pb_z, fi_23, fi_55, gh0_14, \
                         gh0_17, gh1_14, gh1_17, gi_53, gi_54, gi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_12 * fi_55[k]
                   + f_7 * gh0_17[k]
                   - f_8 * gh1_17[k]
                   + pb_x[k] * gi_55[k];

        t_115[k] = pb_z[k] * gi_53[k];

        t_116[k] = f_12 * fi_23[k]
                   + pb_y[k] * gi_54[k];

        t_117[k] = f_5 * gh0_14[k]
                   - f_6 * gh1_14[k]
                   + pb_z[k] * gi_54[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_z, fi_57, gh0_15, gh0_20, gh1_15, \
                         gh1_20, gi_55, gi_56, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * fi_57[k]
                   + f_5 * gh0_20[k]
                   - f_6 * gh1_20[k]
                   + pb_x[k] * gi_58[k];

        t_119[k] = pb_z[k] * gi_55[k];

        t_120[k] = f_3 * gh0_15[k]
                   - f_4 * gh1_15[k]
                   + pb_z[k] * gi_56[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, fi_25, fi_59, gh0_16, \
                         gh0_21, gh1_16, gh1_21, gi_57, gi_58, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * fi_25[k]
                   + pb_y[k] * gi_57[k];

        t_122[k] = f_7 * gh0_16[k]
                   - f_8 * gh1_16[k]
                   + pb_z[k] * gi_57[k];

        t_123[k] = f_12 * fi_59[k]
                   + f_3 * gh0_21[k]
                   - f_4 * gh1_21[k]
                   + pb_x[k] * gi_62[k];

        t_124[k] = pb_z[k] * gi_58[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, fi_27, gh0_17, gh0_18, \
                         gh0_19, gh1_17, gh1_18, gh1_19, gi_59, gi_60, \
                         gi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * gh0_17[k]
                   - f_4 * gh1_17[k]
                   + pb_z[k] * gi_59[k];

        t_126[k] = f_5 * gh0_18[k]
                   - f_6 * gh1_18[k]
                   + pb_z[k] * gi_60[k];

        t_127[k] = f_12 * fi_27[k]
                   + pb_y[k] * gi_61[k];

        t_128[k] = f_9 * gh0_19[k]
                   - f_10 * gh1_19[k]
                   + pb_z[k] * gi_61[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, fi_60, fi_61, fi_62, \
                         fi_63, gi_62, gi_63, gi_65, gi_66, gi_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * fi_60[k]
                   + pb_x[k] * gi_63[k];

        t_130[k] = pb_z[k] * gi_62[k];

        t_131[k] = f_12 * fi_61[k]
                   + pb_x[k] * gi_65[k];

        t_132[k] = f_12 * fi_62[k]
                   + pb_x[k] * gi_66[k];

        t_133[k] = f_12 * fi_63[k]
                   + pb_x[k] * gi_67[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_z, dk0_1, dk1_1, fi_64, \
                         fi_65, fk_50, gi_63, gi_68, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_12 * fi_64[k]
                   + pb_x[k] * gi_68[k];

        t_135[k] = f_12 * fi_65[k]
                   + pb_x[k] * gi_69[k];

        t_136[k] = f_16 * dk0_1[k]
                   - f_17 * dk1_1[k]
                   + pa_x[k] * fk_50[k];

        t_137[k] = pb_z[k] * gi_63[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_z, gh0_21, gh0_22, gh0_23, gh1_21, gh1_22, \
                         gh1_23, gi_64, gi_65, gi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * gh0_21[k]
                   - f_4 * gh1_21[k]
                   + pb_z[k] * gi_64[k];

        t_139[k] = f_5 * gh0_22[k]
                   - f_6 * gh1_22[k]
                   + pb_z[k] * gi_65[k];

        t_140[k] = f_7 * gh0_23[k]
                   - f_8 * gh1_23[k]
                   + pb_z[k] * gi_66[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, fi_33, fk_28, gh0_24, \
                         gh0_25, gh1_24, gh1_25, gi_67, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * gh0_24[k]
                   - f_10 * gh1_24[k]
                   + pb_z[k] * gi_67[k];

        t_142[k] = f_12 * fi_33[k]
                   + pb_y[k] * gi_69[k];

        t_143[k] = f_1 * gh0_25[k]
                   - f_2 * gh1_25[k]
                   + pb_z[k] * gi_69[k];

        t_144[k] = pa_y[k] * fk_28[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, fi_35, \
                         fk_21, fk_22, fk_23, fk_29, fk_30, gi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * fk_21[k];

        t_146[k] = pa_y[k] * fk_29[k];

        t_147[k] = pa_z[k] * fk_22[k];

        t_148[k] = f_11 * fi_35[k]
                   + pb_y[k] * gi_70[k];

        t_149[k] = pa_y[k] * fk_30[k];

        t_150[k] = pa_z[k] * fk_23[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, fi_22, fi_37, \
                         fk_24, fk_31, gi_71, gi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * fi_22[k]
                   + pb_z[k] * gi_71[k];

        t_152[k] = f_11 * fi_37[k]
                   + pb_y[k] * gi_72[k];

        t_153[k] = pa_y[k] * fk_31[k];

        t_154[k] = pa_z[k] * fk_24[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, fi_24, fi_39, fi_40, \
                         fk_32, fk_33, gi_73, gi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * fi_24[k]
                   + pb_z[k] * gi_73[k];

        t_156[k] = f_12 * fi_39[k]
                   + pa_y[k] * fk_32[k];

        t_157[k] = f_11 * fi_40[k]
                   + pb_y[k] * gi_74[k];

        t_158[k] = pa_y[k] * fk_33[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, fi_26, fi_42, fi_43, \
                         fk_25, fk_34, fk_35, gi_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * fk_25[k];

        t_160[k] = f_11 * fi_26[k]
                   + pb_z[k] * gi_75[k];

        t_161[k] = f_13 * fi_42[k]
                   + pa_y[k] * fk_34[k];

        t_162[k] = f_12 * fi_43[k]
                   + pa_y[k] * fk_35[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, fi_44, fi_73, \
                         fk_26, fk_36, gi_76, gi_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * fi_44[k]
                   + pb_y[k] * gi_76[k];

        t_164[k] = pa_y[k] * fk_36[k];

        t_165[k] = pa_z[k] * fk_26[k];

        t_166[k] = f_12 * fi_73[k]
                   + pb_x[k] * gi_78[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, fi_74, fi_75, fi_76, \
                         fi_77, fk_37, gi_79, gi_80, gi_81, gi_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_12 * fi_74[k]
                   + pb_x[k] * gi_79[k];

        t_168[k] = f_12 * fi_75[k]
                   + pb_x[k] * gi_80[k];

        t_169[k] = f_12 * fi_76[k]
                   + pb_x[k] * gi_81[k];

        t_170[k] = f_12 * fi_77[k]
                   + pb_x[k] * gi_82[k];

        t_171[k] = pa_y[k] * fk_37[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pa_z, pb_z, fi_28, fi_47, fi_48, \
                         fk_27, fk_38, fk_39, gi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_z[k] * fk_27[k];

        t_173[k] = f_11 * fi_28[k]
                   + pb_z[k] * gi_77[k];

        t_174[k] = f_14 * fi_47[k]
                   + pa_y[k] * fk_38[k];

        t_175[k] = f_0 * fi_48[k]
                   + pa_y[k] * fk_39[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, fi_49, fi_50, fi_51, fk_40, \
                         fk_41, fk_42, gi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * fi_49[k]
                   + pa_y[k] * fk_40[k];

        t_177[k] = f_12 * fi_50[k]
                   + pa_y[k] * fk_41[k];

        t_178[k] = f_11 * fi_51[k]
                   + pb_y[k] * gi_83[k];

        t_179[k] = pa_y[k] * fk_42[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_y, pb_z, dk0_0, dk1_0, fi_34, \
                         fk_28, gh0_26, gh1_26, gi_84, gi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_16 * dk0_0[k]
                   - f_17 * dk1_0[k]
                   + pa_z[k] * fk_28[k];

        t_181[k] = pb_y[k] * gi_84[k];

        t_182[k] = f_12 * fi_34[k]
                   + pb_z[k] * gi_84[k];

        t_183[k] = f_3 * gh0_26[k]
                   - f_4 * gh1_26[k]
                   + pb_y[k] * gi_85[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pb_y, pb_z, fi_36, fi_81, gh0_27, \
                         gh0_29, gh1_27, gh1_29, gi_86, gi_87, gi_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * gi_86[k];

        t_185[k] = f_12 * fi_81[k]
                   + f_9 * gh0_29[k]
                   - f_10 * gh1_29[k]
                   + pb_x[k] * gi_88[k];

        t_186[k] = f_5 * gh0_27[k]
                   - f_6 * gh1_27[k]
                   + pb_y[k] * gi_87[k];

        t_187[k] = f_12 * fi_36[k]
                   + pb_z[k] * gi_87[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pb_x, pb_y, pb_z, fi_38, fi_83, gh0_28, \
                         gh0_32, gh1_28, gh1_32, gi_88, gi_89, gi_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_y[k] * gi_88[k];

        t_189[k] = f_12 * fi_83[k]
                   + f_7 * gh0_32[k]
                   - f_8 * gh1_32[k]
                   + pb_x[k] * gi_91[k];

        t_190[k] = f_7 * gh0_28[k]
                   - f_8 * gh1_28[k]
                   + pb_y[k] * gi_89[k];

        t_191[k] = f_12 * fi_38[k]
                   + pb_z[k] * gi_89[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, fi_85, gh0_29, gh0_33, gh1_29, \
                         gh1_33, gi_90, gi_91, gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * gh0_29[k]
                   - f_4 * gh1_29[k]
                   + pb_y[k] * gi_90[k];

        t_193[k] = pb_y[k] * gi_91[k];

        t_194[k] = f_12 * fi_85[k]
                   + f_5 * gh0_33[k]
                   - f_6 * gh1_33[k]
                   + pb_x[k] * gi_95[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, fi_41, gh0_30, gh0_31, \
                         gh0_32, gh1_30, gh1_31, gh1_32, gi_92, gi_93, \
                         gi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_9 * gh0_30[k]
                   - f_10 * gh1_30[k]
                   + pb_y[k] * gi_92[k];

        t_196[k] = f_12 * fi_41[k]
                   + pb_z[k] * gi_92[k];

        t_197[k] = f_5 * gh0_31[k]
                   - f_6 * gh1_31[k]
                   + pb_y[k] * gi_93[k];

        t_198[k] = f_3 * gh0_32[k]
                   - f_4 * gh1_32[k]
                   + pb_y[k] * gi_94[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_x, pb_y, fi_86, fi_87, fi_88, gh0_38, \
                         gh1_38, gi_95, gi_96, gi_97, gi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * gi_95[k];

        t_200[k] = f_12 * fi_86[k]
                   + f_3 * gh0_38[k]
                   - f_4 * gh1_38[k]
                   + pb_x[k] * gi_96[k];

        t_201[k] = f_12 * fi_87[k]
                   + pb_x[k] * gi_97[k];

        t_202[k] = f_12 * fi_88[k]
                   + pb_x[k] * gi_98[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, pb_y, fi_89, fi_90, fi_91, \
                         fi_92, gi_96, gi_99, gi_100, gi_101, gi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_12 * fi_89[k]
                   + pb_x[k] * gi_99[k];

        t_204[k] = f_12 * fi_90[k]
                   + pb_x[k] * gi_100[k];

        t_205[k] = f_12 * fi_91[k]
                   + pb_x[k] * gi_101[k];

        t_206[k] = pb_y[k] * gi_96[k];

        t_207[k] = f_12 * fi_92[k]
                   + pb_x[k] * gi_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pb_z, fi_45, gh0_34, gh0_35, \
                         gh0_36, gh1_34, gh1_35, gh1_36, gi_97, gi_99, \
                         gi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * gh0_34[k]
                   - f_2 * gh1_34[k]
                   + pb_y[k] * gi_97[k];

        t_209[k] = f_12 * fi_45[k]
                   + pb_z[k] * gi_97[k];

        t_210[k] = f_9 * gh0_35[k]
                   - f_10 * gh1_35[k]
                   + pb_y[k] * gi_99[k];

        t_211[k] = f_7 * gh0_36[k]
                   - f_8 * gh1_36[k]
                   + pb_y[k] * gi_100[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pb_y, dk0_2, dk1_2, fk_58, gh0_37, \
                         gh0_38, gh1_37, gh1_38, gi_101, gi_102, \
                         gi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * gh0_37[k]
                   - f_6 * gh1_37[k]
                   + pb_y[k] * gi_101[k];

        t_213[k] = f_3 * gh0_38[k]
                   - f_4 * gh1_38[k]
                   + pb_y[k] * gi_102[k];

        t_214[k] = pb_y[k] * gi_103[k];

        t_215[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_2[k]
                   + pa_x[k] * fk_58[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pa_x, pb_y, pb_z, fi_52, fi_93, \
                         fi_95, fk_59, fk_61, gi_104, gi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_15 * fi_93[k]
                   + pa_x[k] * fk_59[k];

        t_217[k] = f_13 * fi_52[k]
                   + pb_y[k] * gi_104[k];

        t_218[k] = pb_z[k] * gi_104[k];

        t_219[k] = f_14 * fi_95[k]
                   + pa_x[k] * fk_61[k];

        t_220[k] = pb_z[k] * gi_105[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_x, pb_y, pb_z, fi_54, fi_96, fi_97, \
                         fk_62, fk_63, gi_106, gi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_14 * fi_96[k]
                   + pa_x[k] * fk_62[k];

        t_222[k] = f_0 * fi_97[k]
                   + pa_x[k] * fk_63[k];

        t_223[k] = pb_z[k] * gi_106[k];

        t_224[k] = f_13 * fi_54[k]
                   + pb_y[k] * gi_107[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pb_z, fi_99, fi_100, fi_102, fk_64, \
                         fk_65, fk_66, gi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * fi_99[k]
                   + pa_x[k] * fk_64[k];

        t_226[k] = f_13 * fi_100[k]
                   + pa_x[k] * fk_65[k];

        t_227[k] = pb_z[k] * gi_108[k];

        t_228[k] = f_13 * fi_102[k]
                   + pa_x[k] * fk_66[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_x, pb_y, pb_z, fi_56, fi_103, fi_104, \
                         fk_67, fk_68, gi_109, gi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * fi_56[k]
                   + pb_y[k] * gi_109[k];

        t_230[k] = f_13 * fi_103[k]
                   + pa_x[k] * fk_67[k];

        t_231[k] = f_12 * fi_104[k]
                   + pa_x[k] * fk_68[k];

        t_232[k] = pb_z[k] * gi_110[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pb_y, fi_58, fi_105, fi_106, \
                         fi_107, fk_69, fk_70, fk_71, gi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_12 * fi_105[k]
                   + pa_x[k] * fk_69[k];

        t_234[k] = f_12 * fi_106[k]
                   + pa_x[k] * fk_70[k];

        t_235[k] = f_13 * fi_58[k]
                   + pb_y[k] * gi_111[k];

        t_236[k] = f_12 * fi_107[k]
                   + pa_x[k] * fk_71[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_x, pb_z, fi_108, fi_110, \
                         fi_111, fi_112, gi_112, gi_113, gi_114, gi_115, \
                         gi_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_11 * fi_108[k]
                   + pb_x[k] * gi_113[k];

        t_238[k] = pb_z[k] * gi_112[k];

        t_239[k] = f_11 * fi_110[k]
                   + pb_x[k] * gi_114[k];

        t_240[k] = f_11 * fi_111[k]
                   + pb_x[k] * gi_115[k];

        t_241[k] = f_11 * fi_112[k]
                   + pb_x[k] * gi_116[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pa_x, pb_x, pb_z, fi_113, fi_114, \
                         fk_72, fk_73, gi_113, gi_117, gi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_11 * fi_113[k]
                   + pb_x[k] * gi_117[k];

        t_243[k] = f_11 * fi_114[k]
                   + pb_x[k] * gi_118[k];

        t_244[k] = pa_x[k] * fk_72[k];

        t_245[k] = pb_z[k] * gi_113[k];

        t_246[k] = pa_x[k] * fk_73[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, pa_x, pa_z, fk_43, \
                         fk_44, fk_74, fk_75, fk_76, fk_77, fk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = pa_x[k] * fk_74[k];

        t_248[k] = pa_x[k] * fk_75[k];

        t_249[k] = pa_x[k] * fk_76[k];

        t_250[k] = pa_x[k] * fk_77[k];

        t_251[k] = pa_x[k] * fk_78[k];

        t_252[k] = pa_z[k] * fk_43[k];

        t_253[k] = pa_z[k] * fk_44[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pa_x, pa_z, pb_y, pb_z, fi_52, fi_66, \
                         fi_118, fk_45, fk_79, gi_119, gi_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_11 * fi_52[k]
                   + pb_z[k] * gi_119[k];

        t_255[k] = pa_z[k] * fk_45[k];

        t_256[k] = f_12 * fi_66[k]
                   + pb_y[k] * gi_120[k];

        t_257[k] = f_14 * fi_118[k]
                   + pa_x[k] * fk_79[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pa_x, pa_z, pb_y, pb_z, fi_53, fi_68, \
                         fi_120, fk_46, fk_80, gi_121, gi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_z[k] * fk_46[k];

        t_259[k] = f_11 * fi_53[k]
                   + pb_z[k] * gi_121[k];

        t_260[k] = f_12 * fi_68[k]
                   + pb_y[k] * gi_122[k];

        t_261[k] = f_0 * fi_120[k]
                   + pa_x[k] * fk_80[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, pa_x, pa_z, pb_y, pb_z, fi_55, fi_70, \
                         fi_122, fk_47, fk_81, gi_123, gi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = pa_z[k] * fk_47[k];

        t_263[k] = f_11 * fi_55[k]
                   + pb_z[k] * gi_123[k];

        t_264[k] = f_13 * fi_122[k]
                   + pa_x[k] * fk_81[k];

        t_265[k] = f_12 * fi_70[k]
                   + pb_y[k] * gi_124[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_x, pa_z, pb_z, fi_57, fi_123, fi_124, \
                         fk_48, fk_82, fk_83, gi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_13 * fi_123[k]
                   + pa_x[k] * fk_82[k];

        t_267[k] = pa_z[k] * fk_48[k];

        t_268[k] = f_11 * fi_57[k]
                   + pb_z[k] * gi_125[k];

        t_269[k] = f_12 * fi_124[k]
                   + pa_x[k] * fk_83[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pa_x, pa_z, pb_y, fi_72, fi_125, fi_126, \
                         fk_49, fk_84, fk_85, gi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_12 * fi_125[k]
                   + pa_x[k] * fk_84[k];

        t_271[k] = f_12 * fi_72[k]
                   + pb_y[k] * gi_126[k];

        t_272[k] = f_12 * fi_126[k]
                   + pa_x[k] * fk_85[k];

        t_273[k] = pa_z[k] * fk_49[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_x, fi_128, fi_129, fi_130, \
                         fi_131, fi_132, gi_127, gi_128, gi_129, gi_130, \
                         gi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_11 * fi_128[k]
                   + pb_x[k] * gi_127[k];

        t_275[k] = f_11 * fi_129[k]
                   + pb_x[k] * gi_128[k];

        t_276[k] = f_11 * fi_130[k]
                   + pb_x[k] * gi_129[k];

        t_277[k] = f_11 * fi_131[k]
                   + pb_x[k] * gi_130[k];

        t_278[k] = f_11 * fi_132[k]
                   + pb_x[k] * gi_131[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, t_283, t_284, pa_x, pb_x, fi_133, fk_86, \
                         fk_87, fk_88, fk_89, fk_90, gi_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_11 * fi_133[k]
                   + pb_x[k] * gi_132[k];

        t_280[k] = pa_x[k] * fk_86[k];

        t_281[k] = pa_x[k] * fk_87[k];

        t_282[k] = pa_x[k] * fk_88[k];

        t_283[k] = pa_x[k] * fk_89[k];

        t_284[k] = pa_x[k] * fk_90[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, t_290, pa_x, pa_y, pb_y, fi_78, \
                         fk_51, fk_52, fk_91, fk_92, fk_93, gi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = pa_x[k] * fk_91[k];

        t_286[k] = pa_x[k] * fk_92[k];

        t_287[k] = pa_x[k] * fk_93[k];

        t_288[k] = pa_y[k] * fk_51[k];

        t_289[k] = f_11 * fi_78[k]
                   + pb_y[k] * gi_133[k];

        t_290[k] = pa_y[k] * fk_52[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pa_x, pa_y, pb_y, fi_79, fi_136, fi_138, \
                         fk_53, fk_94, fk_95, gi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_14 * fi_136[k]
                   + pa_x[k] * fk_94[k];

        t_292[k] = f_11 * fi_79[k]
                   + pb_y[k] * gi_134[k];

        t_293[k] = pa_y[k] * fk_53[k];

        t_294[k] = f_0 * fi_138[k]
                   + pa_x[k] * fk_95[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pa_x, pa_y, pb_y, pb_z, fi_67, fi_81, \
                         fi_140, fk_54, fk_96, gi_135, gi_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_12 * fi_67[k]
                   + pb_z[k] * gi_135[k];

        t_296[k] = f_11 * fi_81[k]
                   + pb_y[k] * gi_136[k];

        t_297[k] = pa_y[k] * fk_54[k];

        t_298[k] = f_13 * fi_140[k]
                   + pa_x[k] * fk_96[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pa_x, pa_y, pb_y, pb_z, fi_69, fi_83, \
                         fi_141, fk_55, fk_97, gi_137, gi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_12 * fi_69[k]
                   + pb_z[k] * gi_137[k];

        t_300[k] = f_13 * fi_141[k]
                   + pa_x[k] * fk_97[k];

        t_301[k] = f_11 * fi_83[k]
                   + pb_y[k] * gi_138[k];

        t_302[k] = pa_y[k] * fk_55[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pb_z, fi_71, fi_143, fi_144, \
                         fi_145, fk_98, fk_99, fk_100, gi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_12 * fi_143[k]
                   + pa_x[k] * fk_98[k];

        t_304[k] = f_12 * fi_71[k]
                   + pb_z[k] * gi_139[k];

        t_305[k] = f_12 * fi_144[k]
                   + pa_x[k] * fk_99[k];

        t_306[k] = f_12 * fi_145[k]
                   + pa_x[k] * fk_100[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_y, pb_x, pb_y, fi_85, fi_146, fi_147, \
                         fk_56, gi_140, gi_141, gi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_11 * fi_85[k]
                   + pb_y[k] * gi_140[k];

        t_308[k] = pa_y[k] * fk_56[k];

        t_309[k] = f_11 * fi_146[k]
                   + pb_x[k] * gi_141[k];

        t_310[k] = f_11 * fi_147[k]
                   + pb_x[k] * gi_142[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, t_315, pa_y, pb_x, fi_148, fi_149, \
                         fi_150, fi_151, fk_57, gi_143, gi_144, gi_145, \
                         gi_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_11 * fi_148[k]
                   + pb_x[k] * gi_143[k];

        t_312[k] = f_11 * fi_149[k]
                   + pb_x[k] * gi_144[k];

        t_313[k] = f_11 * fi_150[k]
                   + pb_x[k] * gi_145[k];

        t_314[k] = f_11 * fi_151[k]
                   + pb_x[k] * gi_146[k];

        t_315[k] = pa_y[k] * fk_57[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, t_320, t_321, t_322, pa_x, fk_101, \
                         fk_102, fk_103, fk_104, fk_105, fk_106, \
                         fk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_x[k] * fk_101[k];

        t_317[k] = pa_x[k] * fk_102[k];

        t_318[k] = pa_x[k] * fk_103[k];

        t_319[k] = pa_x[k] * fk_104[k];

        t_320[k] = pa_x[k] * fk_105[k];

        t_321[k] = pa_x[k] * fk_106[k];

        t_322[k] = pa_x[k] * fk_107[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, pa_x, pb_y, pb_z, fi_78, fi_153, \
                         fi_156, fk_108, fk_109, fk_111, gi_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = pa_x[k] * fk_108[k];

        t_324[k] = f_15 * fi_153[k]
                   + pa_x[k] * fk_109[k];

        t_325[k] = pb_y[k] * gi_147[k];

        t_326[k] = f_13 * fi_78[k]
                   + pb_z[k] * gi_147[k];

        t_327[k] = f_14 * fi_156[k]
                   + pa_x[k] * fk_111[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, pa_x, pb_y, pb_z, fi_80, fi_157, \
                         fi_158, fk_112, fk_113, gi_148, gi_149, \
                         gi_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pb_y[k] * gi_148[k];

        t_329[k] = f_14 * fi_157[k]
                   + pa_x[k] * fk_112[k];

        t_330[k] = f_0 * fi_158[k]
                   + pa_x[k] * fk_113[k];

        t_331[k] = f_13 * fi_80[k]
                   + pb_z[k] * gi_149[k];

        t_332[k] = pb_y[k] * gi_150[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pa_x, pb_z, fi_82, fi_160, fi_161, \
                         fi_162, fk_114, fk_115, fk_116, gi_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_0 * fi_160[k]
                   + pa_x[k] * fk_114[k];

        t_334[k] = f_13 * fi_161[k]
                   + pa_x[k] * fk_115[k];

        t_335[k] = f_13 * fi_82[k]
                   + pb_z[k] * gi_151[k];

        t_336[k] = f_13 * fi_162[k]
                   + pa_x[k] * fk_116[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_x, pb_y, pb_z, fi_84, fi_164, fi_165, \
                         fk_117, fk_118, gi_152, gi_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = pb_y[k] * gi_152[k];

        t_338[k] = f_13 * fi_164[k]
                   + pa_x[k] * fk_117[k];

        t_339[k] = f_12 * fi_165[k]
                   + pa_x[k] * fk_118[k];

        t_340[k] = f_13 * fi_84[k]
                   + pb_z[k] * gi_153[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_x, pb_y, fi_166, fi_167, fi_168, \
                         fk_119, fk_120, fk_121, gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_12 * fi_166[k]
                   + pa_x[k] * fk_119[k];

        t_342[k] = f_12 * fi_167[k]
                   + pa_x[k] * fk_120[k];

        t_343[k] = pb_y[k] * gi_154[k];

        t_344[k] = f_12 * fi_168[k]
                   + pa_x[k] * fk_121[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pb_x, fi_169, fi_170, fi_171, \
                         fi_172, fi_173, gi_156, gi_157, gi_158, gi_159, \
                         gi_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_11 * fi_169[k]
                   + pb_x[k] * gi_156[k];

        t_346[k] = f_11 * fi_170[k]
                   + pb_x[k] * gi_157[k];

        t_347[k] = f_11 * fi_171[k]
                   + pb_x[k] * gi_158[k];

        t_348[k] = f_11 * fi_172[k]
                   + pb_x[k] * gi_159[k];

        t_349[k] = f_11 * fi_173[k]
                   + pb_x[k] * gi_160[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, pa_x, pb_x, pb_y, fi_175, \
                         fk_122, fk_123, fk_124, fk_125, gi_155, \
                         gi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_y[k] * gi_155[k];

        t_351[k] = f_11 * fi_175[k]
                   + pb_x[k] * gi_161[k];

        t_352[k] = pa_x[k] * fk_122[k];

        t_353[k] = pa_x[k] * fk_123[k];

        t_354[k] = pa_x[k] * fk_124[k];

        t_355[k] = pa_x[k] * fk_125[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, t_360, pa_x, pb_x, pb_y, fk_126, fk_127, \
                         fk_128, gh0_39, gh1_39, gi_161, gi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * fk_126[k];

        t_357[k] = pa_x[k] * fk_127[k];

        t_358[k] = pb_y[k] * gi_161[k];

        t_359[k] = pa_x[k] * fk_128[k];

        t_360[k] = f_1 * gh0_39[k]
                   - f_2 * gh1_39[k]
                   + pb_x[k] * gi_162[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pb_x, pb_y, pb_z, fi_93, gh0_40, gh1_40, \
                         gi_162, gi_163, gi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_0 * fi_93[k]
                   + pb_y[k] * gi_162[k];

        t_362[k] = pb_z[k] * gi_162[k];

        t_363[k] = f_9 * gh0_40[k]
                   - f_10 * gh1_40[k]
                   + pb_x[k] * gi_164[k];

        t_364[k] = pb_z[k] * gi_163[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pb_x, pb_y, pb_z, fi_96, gh0_41, gh0_42, \
                         gh1_41, gh1_42, gi_164, gi_165, gi_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_9 * gh0_41[k]
                   - f_10 * gh1_41[k]
                   + pb_x[k] * gi_165[k];

        t_366[k] = f_7 * gh0_42[k]
                   - f_8 * gh1_42[k]
                   + pb_x[k] * gi_166[k];

        t_367[k] = pb_z[k] * gi_164[k];

        t_368[k] = f_0 * fi_96[k]
                   + pb_y[k] * gi_165[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pb_x, pb_z, gh0_43, gh0_44, gh0_45, \
                         gh1_43, gh1_44, gh1_45, gi_166, gi_167, gi_168, \
                         gi_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_7 * gh0_43[k]
                   - f_8 * gh1_43[k]
                   + pb_x[k] * gi_167[k];

        t_370[k] = f_5 * gh0_44[k]
                   - f_6 * gh1_44[k]
                   + pb_x[k] * gi_168[k];

        t_371[k] = pb_z[k] * gi_166[k];

        t_372[k] = f_5 * gh0_45[k]
                   - f_6 * gh1_45[k]
                   + pb_x[k] * gi_169[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_x, pb_y, pb_z, fi_99, gh0_46, gh0_47, \
                         gh1_46, gh1_47, gi_167, gi_168, gi_170, \
                         gi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * fi_99[k]
                   + pb_y[k] * gi_167[k];

        t_374[k] = f_5 * gh0_46[k]
                   - f_6 * gh1_46[k]
                   + pb_x[k] * gi_170[k];

        t_375[k] = f_3 * gh0_47[k]
                   - f_4 * gh1_47[k]
                   + pb_x[k] * gi_171[k];

        t_376[k] = pb_z[k] * gi_168[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pb_x, pb_y, fi_103, gh0_49, gh0_50, gh1_49, \
                         gh1_50, gi_170, gi_172, gi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_3 * gh0_49[k]
                   - f_4 * gh1_49[k]
                   + pb_x[k] * gi_172[k];

        t_378[k] = f_3 * gh0_50[k]
                   - f_4 * gh1_50[k]
                   + pb_x[k] * gi_173[k];

        t_379[k] = f_0 * fi_103[k]
                   + pb_y[k] * gi_170[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, t_385, pb_x, gh0_51, gh1_51, \
                         gi_174, gi_175, gi_176, gi_177, gi_178, \
                         gi_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_3 * gh0_51[k]
                   - f_4 * gh1_51[k]
                   + pb_x[k] * gi_174[k];

        t_381[k] = pb_x[k] * gi_175[k];

        t_382[k] = pb_x[k] * gi_176[k];

        t_383[k] = pb_x[k] * gi_177[k];

        t_384[k] = pb_x[k] * gi_178[k];

        t_385[k] = pb_x[k] * gi_179[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, pb_x, pb_y, pb_z, fi_108, gh0_47, \
                         gh1_47, gi_175, gi_176, gi_180, gi_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = pb_x[k] * gi_180[k];

        t_387[k] = pb_x[k] * gi_181[k];

        t_388[k] = f_0 * fi_108[k]
                   + f_1 * gh0_47[k]
                   - f_2 * gh1_47[k]
                   + pb_y[k] * gi_175[k];

        t_389[k] = pb_z[k] * gi_175[k];

        t_390[k] = f_3 * gh0_47[k]
                   - f_4 * gh1_47[k]
                   + pb_z[k] * gi_176[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pb_z, gh0_48, gh0_49, gh0_50, gh1_48, gh1_49, \
                         gh1_50, gi_177, gi_178, gi_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_5 * gh0_48[k]
                   - f_6 * gh1_48[k]
                   + pb_z[k] * gi_177[k];

        t_392[k] = f_7 * gh0_49[k]
                   - f_8 * gh1_49[k]
                   + pb_z[k] * gi_178[k];

        t_393[k] = f_9 * gh0_50[k]
                   - f_10 * gh1_50[k]
                   + pb_z[k] * gi_179[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pa_z, pb_y, pb_z, fi_93, fi_114, \
                         fk_59, fk_60, gh0_51, gh1_51, gi_181, gi_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_0 * fi_114[k]
                   + pb_y[k] * gi_181[k];

        t_395[k] = f_1 * gh0_51[k]
                   - f_2 * gh1_51[k]
                   + pb_z[k] * gi_181[k];

        t_396[k] = pa_z[k] * fk_59[k];

        t_397[k] = pa_z[k] * fk_60[k];

        t_398[k] = f_11 * fi_93[k]
                   + pb_z[k] * gi_182[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, pa_z, pb_y, pb_z, fi_94, fi_95, \
                         fi_116, fk_61, fk_62, fk_63, gi_183, gi_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pa_z[k] * fk_61[k];

        t_400[k] = f_13 * fi_116[k]
                   + pb_y[k] * gi_183[k];

        t_401[k] = f_12 * fi_94[k]
                   + pa_z[k] * fk_62[k];

        t_402[k] = pa_z[k] * fk_63[k];

        t_403[k] = f_11 * fi_95[k]
                   + pb_z[k] * gi_184[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_z, pb_y, pb_z, fi_96, fi_97, fi_118, \
                         fk_64, fk_65, gi_185, gi_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_13 * fi_118[k]
                   + pb_y[k] * gi_185[k];

        t_405[k] = f_13 * fi_96[k]
                   + pa_z[k] * fk_64[k];

        t_406[k] = pa_z[k] * fk_65[k];

        t_407[k] = f_11 * fi_97[k]
                   + pb_z[k] * gi_186[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pa_z, pb_y, fi_98, fi_99, fi_120, fk_66, \
                         fk_67, fk_68, gi_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_12 * fi_98[k]
                   + pa_z[k] * fk_66[k];

        t_409[k] = f_13 * fi_120[k]
                   + pb_y[k] * gi_187[k];

        t_410[k] = f_0 * fi_99[k]
                   + pa_z[k] * fk_67[k];

        t_411[k] = pa_z[k] * fk_68[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_z, pb_y, pb_z, fi_100, fi_101, fi_102, \
                         fi_123, fk_69, fk_70, gi_188, gi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_11 * fi_100[k]
                   + pb_z[k] * gi_188[k];

        t_413[k] = f_12 * fi_101[k]
                   + pa_z[k] * fk_69[k];

        t_414[k] = f_13 * fi_102[k]
                   + pa_z[k] * fk_70[k];

        t_415[k] = f_13 * fi_123[k]
                   + pb_y[k] * gi_189[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, t_421, pa_z, pb_x, fi_103, fk_71, \
                         gi_190, gi_191, gi_192, gi_193, gi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_14 * fi_103[k]
                   + pa_z[k] * fk_71[k];

        t_417[k] = pb_x[k] * gi_190[k];

        t_418[k] = pb_x[k] * gi_191[k];

        t_419[k] = pb_x[k] * gi_192[k];

        t_420[k] = pb_x[k] * gi_193[k];

        t_421[k] = pb_x[k] * gi_194[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, pa_z, pb_x, pb_z, fi_108, fi_109, \
                         fk_72, fk_73, gi_190, gi_195, gi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = pb_x[k] * gi_195[k];

        t_423[k] = pb_x[k] * gi_196[k];

        t_424[k] = pa_z[k] * fk_72[k];

        t_425[k] = f_11 * fi_108[k]
                   + pb_z[k] * gi_190[k];

        t_426[k] = f_12 * fi_109[k]
                   + pa_z[k] * fk_73[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pa_z, pb_y, fi_110, fi_111, fi_112, \
                         fi_133, fk_74, fk_75, fk_76, gi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_13 * fi_110[k]
                   + pa_z[k] * fk_74[k];

        t_428[k] = f_0 * fi_111[k]
                   + pa_z[k] * fk_75[k];

        t_429[k] = f_14 * fi_112[k]
                   + pa_z[k] * fk_76[k];

        t_430[k] = f_13 * fi_133[k]
                   + pb_y[k] * gi_196[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pa_z, pb_x, pb_y, pb_z, fi_114, fi_115, \
                         fi_134, fk_78, gh0_52, gh1_52, gi_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_15 * fi_114[k]
                   + pa_z[k] * fk_78[k];

        t_432[k] = f_1 * gh0_52[k]
                   - f_2 * gh1_52[k]
                   + pb_x[k] * gi_197[k];

        t_433[k] = f_12 * fi_134[k]
                   + pb_y[k] * gi_197[k];

        t_434[k] = f_12 * fi_115[k]
                   + pb_z[k] * gi_197[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pb_x, pb_y, fi_135, gh0_53, gh0_54, gh1_53, \
                         gh1_54, gi_198, gi_199, gi_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_9 * gh0_53[k]
                   - f_10 * gh1_53[k]
                   + pb_x[k] * gi_199[k];

        t_436[k] = f_12 * fi_135[k]
                   + pb_y[k] * gi_198[k];

        t_437[k] = f_9 * gh0_54[k]
                   - f_10 * gh1_54[k]
                   + pb_x[k] * gi_200[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pb_x, pb_y, pb_z, fi_117, fi_137, gh0_55, \
                         gh1_55, gi_199, gi_200, gi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_7 * gh0_55[k]
                   - f_8 * gh1_55[k]
                   + pb_x[k] * gi_201[k];

        t_439[k] = f_12 * fi_117[k]
                   + pb_z[k] * gi_199[k];

        t_440[k] = f_12 * fi_137[k]
                   + pb_y[k] * gi_200[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pb_x, pb_z, fi_119, gh0_56, gh0_57, gh1_56, \
                         gh1_57, gi_201, gi_202, gi_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_7 * gh0_56[k]
                   - f_8 * gh1_56[k]
                   + pb_x[k] * gi_202[k];

        t_442[k] = f_5 * gh0_57[k]
                   - f_6 * gh1_57[k]
                   + pb_x[k] * gi_203[k];

        t_443[k] = f_12 * fi_119[k]
                   + pb_z[k] * gi_201[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_x, pb_y, fi_139, gh0_58, gh0_59, gh1_58, \
                         gh1_59, gi_202, gi_204, gi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_5 * gh0_58[k]
                   - f_6 * gh1_58[k]
                   + pb_x[k] * gi_204[k];

        t_445[k] = f_12 * fi_139[k]
                   + pb_y[k] * gi_202[k];

        t_446[k] = f_5 * gh0_59[k]
                   - f_6 * gh1_59[k]
                   + pb_x[k] * gi_205[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pb_x, pb_z, fi_121, gh0_60, gh0_61, gh1_60, \
                         gh1_61, gi_203, gi_206, gi_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_3 * gh0_60[k]
                   - f_4 * gh1_60[k]
                   + pb_x[k] * gi_206[k];

        t_448[k] = f_12 * fi_121[k]
                   + pb_z[k] * gi_203[k];

        t_449[k] = f_3 * gh0_61[k]
                   - f_4 * gh1_61[k]
                   + pb_x[k] * gi_207[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pb_x, pb_y, fi_142, gh0_62, gh0_64, \
                         gh1_62, gh1_64, gi_205, gi_208, gi_209, \
                         gi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_3 * gh0_62[k]
                   - f_4 * gh1_62[k]
                   + pb_x[k] * gi_208[k];

        t_451[k] = f_12 * fi_142[k]
                   + pb_y[k] * gi_205[k];

        t_452[k] = f_3 * gh0_64[k]
                   - f_4 * gh1_64[k]
                   + pb_x[k] * gi_209[k];

        t_453[k] = pb_x[k] * gi_210[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, t_459, pb_x, gi_211, gi_212, \
                         gi_213, gi_214, gi_215, gi_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pb_x[k] * gi_211[k];

        t_455[k] = pb_x[k] * gi_212[k];

        t_456[k] = pb_x[k] * gi_213[k];

        t_457[k] = pb_x[k] * gi_214[k];

        t_458[k] = pb_x[k] * gi_215[k];

        t_459[k] = pb_x[k] * gi_216[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pa_z, pb_y, pb_z, dk0_1, dk1_1, fi_127, fi_148, \
                         fk_86, gh0_61, gh1_61, gi_210, gi_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_16 * dk0_1[k]
                   - f_17 * dk1_1[k]
                   + pa_z[k] * fk_86[k];

        t_461[k] = f_12 * fi_127[k]
                   + pb_z[k] * gi_210[k];

        t_462[k] = f_12 * fi_148[k]
                   + f_9 * gh0_61[k]
                   - f_10 * gh1_61[k]
                   + pb_y[k] * gi_212[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pb_y, fi_149, fi_150, fi_151, gh0_62, gh0_63, \
                         gh0_64, gh1_62, gh1_63, gh1_64, gi_213, gi_214, \
                         gi_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_12 * fi_149[k]
                   + f_7 * gh0_62[k]
                   - f_8 * gh1_62[k]
                   + pb_y[k] * gi_213[k];

        t_464[k] = f_12 * fi_150[k]
                   + f_5 * gh0_63[k]
                   - f_6 * gh1_63[k]
                   + pb_y[k] * gi_214[k];

        t_465[k] = f_12 * fi_151[k]
                   + f_3 * gh0_64[k]
                   - f_4 * gh1_64[k]
                   + pb_y[k] * gi_215[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, pa_y, pb_y, dk0_2, dk1_2, fi_152, \
                         fi_153, fk_108, fk_109, fk_110, gi_216, \
                         gi_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_12 * fi_152[k]
                   + pb_y[k] * gi_216[k];

        t_467[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_2[k]
                   + pa_y[k] * fk_108[k];

        t_468[k] = pa_y[k] * fk_109[k];

        t_469[k] = f_11 * fi_153[k]
                   + pb_y[k] * gi_217[k];

        t_470[k] = pa_y[k] * fk_110[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pa_y, pb_y, fi_154, fi_155, fi_156, \
                         fk_111, fk_112, fk_113, gi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_12 * fi_154[k]
                   + pa_y[k] * fk_111[k];

        t_472[k] = f_11 * fi_155[k]
                   + pb_y[k] * gi_218[k];

        t_473[k] = pa_y[k] * fk_112[k];

        t_474[k] = f_13 * fi_156[k]
                   + pa_y[k] * fk_113[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_y, pb_y, pb_z, fi_136, fi_157, fi_158, \
                         fk_114, fk_115, gi_219, gi_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_13 * fi_136[k]
                   + pb_z[k] * gi_219[k];

        t_476[k] = f_11 * fi_157[k]
                   + pb_y[k] * gi_220[k];

        t_477[k] = pa_y[k] * fk_114[k];

        t_478[k] = f_0 * fi_158[k]
                   + pa_y[k] * fk_115[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_y, pb_y, pb_z, fi_138, fi_159, fi_160, \
                         fk_116, fk_117, gi_221, gi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * fi_138[k]
                   + pb_z[k] * gi_221[k];

        t_480[k] = f_12 * fi_159[k]
                   + pa_y[k] * fk_116[k];

        t_481[k] = f_11 * fi_160[k]
                   + pb_y[k] * gi_222[k];

        t_482[k] = pa_y[k] * fk_117[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_y, pb_z, fi_140, fi_161, fi_162, \
                         fi_163, fk_118, fk_119, fk_120, gi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_14 * fi_161[k]
                   + pa_y[k] * fk_118[k];

        t_484[k] = f_13 * fi_140[k]
                   + pb_z[k] * gi_223[k];

        t_485[k] = f_13 * fi_162[k]
                   + pa_y[k] * fk_119[k];

        t_486[k] = f_12 * fi_163[k]
                   + pa_y[k] * fk_120[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, t_492, pa_y, pb_x, pb_y, fi_164, \
                         fk_121, gi_224, gi_225, gi_226, gi_227, \
                         gi_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_11 * fi_164[k]
                   + pb_y[k] * gi_224[k];

        t_488[k] = pa_y[k] * fk_121[k];

        t_489[k] = pb_x[k] * gi_225[k];

        t_490[k] = pb_x[k] * gi_226[k];

        t_491[k] = pb_x[k] * gi_227[k];

        t_492[k] = pb_x[k] * gi_228[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, pa_y, pb_x, pb_z, fi_146, fi_169, \
                         fk_122, gi_225, gi_229, gi_230, gi_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = pb_x[k] * gi_229[k];

        t_494[k] = pb_x[k] * gi_230[k];

        t_495[k] = pb_x[k] * gi_231[k];

        t_496[k] = f_15 * fi_169[k]
                   + pa_y[k] * fk_122[k];

        t_497[k] = f_13 * fi_146[k]
                   + pb_z[k] * gi_225[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pa_y, fi_171, fi_172, fi_173, fi_174, \
                         fk_124, fk_125, fk_126, fk_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_14 * fi_171[k]
                   + pa_y[k] * fk_124[k];

        t_499[k] = f_0 * fi_172[k]
                   + pa_y[k] * fk_125[k];

        t_500[k] = f_13 * fi_173[k]
                   + pa_y[k] * fk_126[k];

        t_501[k] = f_12 * fi_174[k]
                   + pa_y[k] * fk_127[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, pa_y, pb_x, pb_y, pb_z, fi_153, \
                         fi_175, fk_128, gh0_65, gh1_65, gi_231, \
                         gi_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_11 * fi_175[k]
                   + pb_y[k] * gi_231[k];

        t_503[k] = pa_y[k] * fk_128[k];

        t_504[k] = f_1 * gh0_65[k]
                   - f_2 * gh1_65[k]
                   + pb_x[k] * gi_232[k];

        t_505[k] = pb_y[k] * gi_232[k];

        t_506[k] = f_0 * fi_153[k]
                   + pb_z[k] * gi_232[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pb_x, pb_y, gh0_66, gh0_67, gh0_68, \
                         gh1_66, gh1_67, gh1_68, gi_233, gi_234, gi_235, \
                         gi_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_9 * gh0_66[k]
                   - f_10 * gh1_66[k]
                   + pb_x[k] * gi_234[k];

        t_508[k] = pb_y[k] * gi_233[k];

        t_509[k] = f_9 * gh0_67[k]
                   - f_10 * gh1_67[k]
                   + pb_x[k] * gi_235[k];

        t_510[k] = f_7 * gh0_68[k]
                   - f_8 * gh1_68[k]
                   + pb_x[k] * gi_236[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pb_x, pb_y, pb_z, fi_156, gh0_69, gh0_70, \
                         gh1_69, gh1_70, gi_234, gi_235, gi_237, \
                         gi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_0 * fi_156[k]
                   + pb_z[k] * gi_234[k];

        t_512[k] = pb_y[k] * gi_235[k];

        t_513[k] = f_7 * gh0_69[k]
                   - f_8 * gh1_69[k]
                   + pb_x[k] * gi_237[k];

        t_514[k] = f_5 * gh0_70[k]
                   - f_6 * gh1_70[k]
                   + pb_x[k] * gi_238[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_x, pb_y, pb_z, fi_158, gh0_71, gh0_72, \
                         gh1_71, gh1_72, gi_236, gi_237, gi_239, \
                         gi_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_0 * fi_158[k]
                   + pb_z[k] * gi_236[k];

        t_516[k] = f_5 * gh0_71[k]
                   - f_6 * gh1_71[k]
                   + pb_x[k] * gi_239[k];

        t_517[k] = pb_y[k] * gi_237[k];

        t_518[k] = f_5 * gh0_72[k]
                   - f_6 * gh1_72[k]
                   + pb_x[k] * gi_240[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pb_x, pb_z, fi_161, gh0_73, gh0_74, gh1_73, \
                         gh1_74, gi_238, gi_241, gi_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_3 * gh0_73[k]
                   - f_4 * gh1_73[k]
                   + pb_x[k] * gi_241[k];

        t_520[k] = f_0 * fi_161[k]
                   + pb_z[k] * gi_238[k];

        t_521[k] = f_3 * gh0_74[k]
                   - f_4 * gh1_74[k]
                   + pb_x[k] * gi_242[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pb_x, pb_y, gh0_75, gh0_77, \
                         gh1_75, gh1_77, gi_240, gi_243, gi_244, gi_245, \
                         gi_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_3 * gh0_75[k]
                   - f_4 * gh1_75[k]
                   + pb_x[k] * gi_243[k];

        t_523[k] = pb_y[k] * gi_240[k];

        t_524[k] = f_3 * gh0_77[k]
                   - f_4 * gh1_77[k]
                   + pb_x[k] * gi_244[k];

        t_525[k] = pb_x[k] * gi_245[k];

        t_526[k] = pb_x[k] * gi_246[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, t_532, pb_x, pb_y, gh0_73, gh1_73, \
                         gi_245, gi_247, gi_248, gi_249, gi_250, \
                         gi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = pb_x[k] * gi_247[k];

        t_528[k] = pb_x[k] * gi_248[k];

        t_529[k] = pb_x[k] * gi_249[k];

        t_530[k] = pb_x[k] * gi_250[k];

        t_531[k] = pb_x[k] * gi_251[k];

        t_532[k] = f_1 * gh0_73[k]
                   - f_2 * gh1_73[k]
                   + pb_y[k] * gi_245[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pb_y, pb_z, fi_169, gh0_74, gh0_75, gh1_74, \
                         gh1_75, gi_245, gi_247, gi_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_0 * fi_169[k]
                   + pb_z[k] * gi_245[k];

        t_534[k] = f_9 * gh0_74[k]
                   - f_10 * gh1_74[k]
                   + pb_y[k] * gi_247[k];

        t_535[k] = f_7 * gh0_75[k]
                   - f_8 * gh1_75[k]
                   + pb_y[k] * gi_248[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pb_y, pb_z, fi_175, gh0_76, gh0_77, \
                         gh1_76, gh1_77, gi_249, gi_250, gi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_5 * gh0_76[k]
                   - f_6 * gh1_76[k]
                   + pb_y[k] * gi_249[k];

        t_537[k] = f_3 * gh0_77[k]
                   - f_4 * gh1_77[k]
                   + pb_y[k] * gi_250[k];

        t_538[k] = pb_y[k] * gi_251[k];

        t_539[k] = f_0 * fi_175[k]
                   + f_1 * gh0_77[k]
                   - f_2 * gh1_77[k]
                   + pb_z[k] * gi_251[k];
    }
}

auto
compute_prim_gk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_54 = buffer.data(dk1 + 54);
    const auto *dk1_95 = buffer.data(dk1 + 95);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
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
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
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
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
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
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_71 = buffer.data(fk + 71);
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
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_128 = buffer.data(fk + 128);

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
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_41 = buffer.data(gh0 + 41);
    const auto *gh0_42 = buffer.data(gh0 + 42);
    const auto *gh0_43 = buffer.data(gh0 + 43);
    const auto *gh0_44 = buffer.data(gh0 + 44);
    const auto *gh0_45 = buffer.data(gh0 + 45);
    const auto *gh0_46 = buffer.data(gh0 + 46);
    const auto *gh0_47 = buffer.data(gh0 + 47);
    const auto *gh0_48 = buffer.data(gh0 + 48);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_74 = buffer.data(gh0 + 74);
    const auto *gh0_75 = buffer.data(gh0 + 75);
    const auto *gh0_76 = buffer.data(gh0 + 76);
    const auto *gh0_77 = buffer.data(gh0 + 77);

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
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_42 = buffer.data(gh1 + 42);
    const auto *gh1_43 = buffer.data(gh1 + 43);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);
    const auto *gh1_48 = buffer.data(gh1 + 48);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_62 = buffer.data(gh1 + 62);
    const auto *gh1_63 = buffer.data(gh1 + 63);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_74 = buffer.data(gh1 + 74);
    const auto *gh1_75 = buffer.data(gh1 + 75);
    const auto *gh1_76 = buffer.data(gh1 + 76);
    const auto *gh1_77 = buffer.data(gh1 + 77);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, gi_0, \
                         gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, gh1_2, \
                         gh1_3, gi_3, gi_4, gi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_6[k] = pb_y[k] * gi_4[k];

        t_7[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_8[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, gh0_4, gh0_5, gh1_4, gh1_5, gi_6, \
                         gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * gh0_4[k]
                 - f_4 * gh1_4[k]
                 + pb_y[k] * gi_6[k];

        t_10[k] = pb_y[k] * gi_7[k];

        t_11[k] = f_7 * gh0_4[k]
                  - f_8 * gh1_4[k]
                  + pb_z[k] * gi_7[k];

        t_12[k] = f_9 * gh0_5[k]
                  - f_10 * gh1_5[k]
                  + pb_y[k] * gi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, gh0_6, gh0_7, gh1_6, gh1_7, gi_9, \
                         gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gh0_6[k]
                  - f_6 * gh1_6[k]
                  + pb_y[k] * gi_9[k];

        t_14[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_15[k] = pb_y[k] * gi_11[k];

        t_16[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_x, pb_y, fi_12, fi_18, gh0_8, gh0_9, \
                         gh1_8, gh1_9, gi_12, gi_13, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_12[k]
                  + pb_x[k] * gi_12[k];

        t_18[k] = f_0 * fi_18[k]
                  + pb_x[k] * gi_17[k];

        t_19[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_20[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_9[k]
                  + pb_y[k] * gi_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, gh0_10, gh0_11, gh0_12, gh1_10, gh1_11, \
                         gh1_12, gi_14, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_10[k]
                  + pb_y[k] * gi_14[k];

        t_22[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_11[k]
                  + pb_y[k] * gi_15[k];

        t_23[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_12[k]
                  + pb_y[k] * gi_16[k];

        t_24[k] = pb_y[k] * gi_17[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_y, pb_z, fi_0, fi_1, fk_0, fk_3, \
                         gh0_12, gh1_12, gi_17, gi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_12[k]
                  + pb_z[k] * gi_17[k];

        t_26[k] = pa_y[k] * fk_0[k];

        t_27[k] = f_11 * fi_0[k]
                  + pb_y[k] * gi_18[k];

        t_28[k] = f_12 * fi_1[k]
                  + pa_y[k] * fk_3[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_y, fi_3, fi_5, fi_8, fk_4, \
                         fk_5, fk_7, fk_8, fk_11, fk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * fk_4[k];

        t_30[k] = f_13 * fi_3[k]
                  + pa_y[k] * fk_5[k];

        t_31[k] = pa_y[k] * fk_7[k];

        t_32[k] = f_0 * fi_5[k]
                  + pa_y[k] * fk_8[k];

        t_33[k] = pa_y[k] * fk_11[k];

        t_34[k] = f_14 * fi_8[k]
                  + pa_y[k] * fk_12[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pb_x, fi_12, fi_14, fi_15, fi_20, \
                         fk_16, fk_17, fk_18, fk_19, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_y[k] * fk_16[k];

        t_36[k] = f_13 * fi_20[k]
                  + pb_x[k] * gi_19[k];

        t_37[k] = f_15 * fi_12[k]
                  + pa_y[k] * fk_17[k];

        t_38[k] = f_14 * fi_14[k]
                  + pa_y[k] * fk_18[k];

        t_39[k] = f_0 * fi_15[k]
                  + pa_y[k] * fk_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, fi_16, fi_17, fi_18, \
                         fk_0, fk_20, fk_21, fk_23, gi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * fi_16[k]
                  + pa_y[k] * fk_20[k];

        t_41[k] = f_12 * fi_17[k]
                  + pa_y[k] * fk_21[k];

        t_42[k] = f_11 * fi_18[k]
                  + pb_y[k] * gi_20[k];

        t_43[k] = pa_y[k] * fk_23[k];

        t_44[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_z, pb_z, fi_0, fi_2, fi_4, fk_3, \
                         fk_4, fk_5, fk_7, gi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_11 * fi_0[k]
                  + pb_z[k] * gi_21[k];

        t_46[k] = pa_z[k] * fk_3[k];

        t_47[k] = f_12 * fi_2[k]
                  + pa_z[k] * fk_4[k];

        t_48[k] = pa_z[k] * fk_5[k];

        t_49[k] = f_13 * fi_4[k]
                  + pa_z[k] * fk_7[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_x, fi_7, fi_11, fi_31, fk_8, \
                         fk_11, fk_12, fk_16, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_z[k] * fk_8[k];

        t_51[k] = f_0 * fi_7[k]
                  + pa_z[k] * fk_11[k];

        t_52[k] = pa_z[k] * fk_12[k];

        t_53[k] = f_14 * fi_11[k]
                  + pa_z[k] * fk_16[k];

        t_54[k] = f_13 * fi_31[k]
                  + pb_x[k] * gi_23[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_z, pb_z, fi_12, fi_13, fi_14, fi_15, \
                         fk_17, fk_18, fk_19, fk_20, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * fk_17[k];

        t_56[k] = f_11 * fi_12[k]
                  + pb_z[k] * gi_22[k];

        t_57[k] = f_12 * fi_13[k]
                  + pa_z[k] * fk_18[k];

        t_58[k] = f_13 * fi_14[k]
                  + pa_z[k] * fk_19[k];

        t_59[k] = f_0 * fi_15[k]
                  + pa_z[k] * fk_20[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, pb_y, dk0_0, dk1_0, fi_16, fi_18, \
                         fi_19, fk_21, fk_23, fk_24, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_14 * fi_16[k]
                  + pa_z[k] * fk_21[k];

        t_61[k] = f_15 * fi_18[k]
                  + pa_z[k] * fk_23[k];

        t_62[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_y[k] * fk_24[k];

        t_63[k] = f_12 * fi_19[k]
                  + pb_y[k] * gi_24[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_z, fi_33, gh0_13, gh0_15, gh1_13, gh1_15, \
                         gi_24, gi_25, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * gi_24[k];

        t_65[k] = f_12 * fi_33[k]
                  + f_9 * gh0_15[k]
                  - f_10 * gh1_15[k]
                  + pb_x[k] * gi_26[k];

        t_66[k] = f_3 * gh0_13[k]
                  - f_4 * gh1_13[k]
                  + pb_z[k] * gi_25[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pb_z, fi_34, gh0_14, gh0_17, gh1_14, gh1_17, \
                         gi_26, gi_27, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_12 * fi_34[k]
                  + f_7 * gh0_17[k]
                  - f_8 * gh1_17[k]
                  + pb_x[k] * gi_28[k];

        t_68[k] = pb_z[k] * gi_26[k];

        t_69[k] = f_5 * gh0_14[k]
                  - f_6 * gh1_14[k]
                  + pb_z[k] * gi_27[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_x, pb_z, fi_35, gh0_15, gh0_20, gh1_15, gh1_20, \
                         gi_28, gi_29, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_12 * fi_35[k]
                  + f_5 * gh0_20[k]
                  - f_6 * gh1_20[k]
                  + pb_x[k] * gi_31[k];

        t_71[k] = pb_z[k] * gi_28[k];

        t_72[k] = f_3 * gh0_15[k]
                  - f_4 * gh1_15[k]
                  + pb_z[k] * gi_29[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, pb_z, fi_36, gh0_16, gh0_21, gh1_16, gh1_21, \
                         gi_30, gi_31, gi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * gh0_16[k]
                  - f_8 * gh1_16[k]
                  + pb_z[k] * gi_30[k];

        t_74[k] = f_12 * fi_36[k]
                  + f_3 * gh0_21[k]
                  - f_4 * gh1_21[k]
                  + pb_x[k] * gi_35[k];

        t_75[k] = pb_z[k] * gi_31[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_z, gh0_17, gh0_18, gh0_19, gh1_17, gh1_18, \
                         gh1_19, gi_32, gi_33, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * gh0_17[k]
                  - f_4 * gh1_17[k]
                  + pb_z[k] * gi_32[k];

        t_77[k] = f_5 * gh0_18[k]
                  - f_6 * gh1_18[k]
                  + pb_z[k] * gi_33[k];

        t_78[k] = f_9 * gh0_19[k]
                  - f_10 * gh1_19[k]
                  + pb_z[k] * gi_34[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_x, pb_x, pb_z, dk0_1, dk1_54, fi_37, \
                         fk_46, gh0_21, gh1_21, gi_36, gi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_12 * fi_37[k]
                  + pb_x[k] * gi_36[k];

        t_80[k] = f_16 * dk0_1[k]
                  - f_17 * dk1_54[k]
                  + pa_x[k] * fk_46[k];

        t_81[k] = pb_z[k] * gi_36[k];

        t_82[k] = f_3 * gh0_21[k]
                  - f_4 * gh1_21[k]
                  + pb_z[k] * gi_37[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_z, gh0_22, gh0_23, gh0_24, gh1_22, gh1_23, \
                         gh1_24, gi_38, gi_39, gi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_5 * gh0_22[k]
                  - f_6 * gh1_22[k]
                  + pb_z[k] * gi_38[k];

        t_84[k] = f_7 * gh0_23[k]
                  - f_8 * gh1_23[k]
                  + pb_z[k] * gi_39[k];

        t_85[k] = f_9 * gh0_24[k]
                  - f_10 * gh1_24[k]
                  + pb_z[k] * gi_40[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_y, pa_z, pb_y, pb_z, fi_21, fk_25, \
                         fk_31, fk_32, gh0_25, gh1_25, gi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_12 * fi_21[k]
                  + pb_y[k] * gi_41[k];

        t_87[k] = f_1 * gh0_25[k]
                  - f_2 * gh1_25[k]
                  + pb_z[k] * gi_41[k];

        t_88[k] = pa_y[k] * fk_31[k];

        t_89[k] = pa_z[k] * fk_25[k];

        t_90[k] = pa_y[k] * fk_32[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, t_96, t_97, pa_y, pa_z, fk_26, fk_27, \
                         fk_28, fk_29, fk_33, fk_34, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pa_z[k] * fk_26[k];

        t_92[k] = pa_y[k] * fk_33[k];

        t_93[k] = pa_z[k] * fk_27[k];

        t_94[k] = pa_y[k] * fk_34[k];

        t_95[k] = pa_z[k] * fk_28[k];

        t_96[k] = pa_y[k] * fk_35[k];

        t_97[k] = pa_z[k] * fk_29[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_y, pb_z, fi_20, fi_27, fi_28, fi_29, \
                         fk_36, fk_37, fk_38, gi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_11 * fi_20[k]
                  + pb_z[k] * gi_42[k];

        t_99[k] = f_14 * fi_27[k]
                  + pa_y[k] * fk_36[k];

        t_100[k] = f_0 * fi_28[k]
                   + pa_y[k] * fk_37[k];

        t_101[k] = f_13 * fi_29[k]
                   + pa_y[k] * fk_38[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_y, pa_z, pb_y, dk0_0, dk1_0, fi_30, \
                         fi_31, fk_30, fk_39, fk_40, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * fi_30[k]
                   + pa_y[k] * fk_39[k];

        t_103[k] = f_11 * fi_31[k]
                   + pb_y[k] * gi_43[k];

        t_104[k] = pa_y[k] * fk_40[k];

        t_105[k] = f_16 * dk0_0[k]
                   - f_17 * dk1_0[k]
                   + pa_z[k] * fk_30[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_x, pb_y, pb_z, fi_22, fi_39, gh0_26, \
                         gh0_29, gh1_26, gh1_29, gi_44, gi_45, gi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_y[k] * gi_44[k];

        t_107[k] = f_12 * fi_22[k]
                   + pb_z[k] * gi_44[k];

        t_108[k] = f_3 * gh0_26[k]
                   - f_4 * gh1_26[k]
                   + pb_y[k] * gi_45[k];

        t_109[k] = f_12 * fi_39[k]
                   + f_9 * gh0_29[k]
                   - f_10 * gh1_29[k]
                   + pb_x[k] * gi_47[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, pb_y, fi_40, gh0_27, gh0_32, gh1_27, \
                         gh1_32, gi_46, gi_47, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_5 * gh0_27[k]
                   - f_6 * gh1_27[k]
                   + pb_y[k] * gi_46[k];

        t_111[k] = pb_y[k] * gi_47[k];

        t_112[k] = f_12 * fi_40[k]
                   + f_7 * gh0_32[k]
                   - f_8 * gh1_32[k]
                   + pb_x[k] * gi_50[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, gh0_28, gh0_29, gh1_28, gh1_29, gi_48, \
                         gi_49, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_7 * gh0_28[k]
                   - f_8 * gh1_28[k]
                   + pb_y[k] * gi_48[k];

        t_114[k] = f_3 * gh0_29[k]
                   - f_4 * gh1_29[k]
                   + pb_y[k] * gi_49[k];

        t_115[k] = pb_y[k] * gi_50[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_x, pb_y, fi_41, gh0_30, gh0_31, gh0_33, \
                         gh1_30, gh1_31, gh1_33, gi_51, gi_52, gi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_12 * fi_41[k]
                   + f_5 * gh0_33[k]
                   - f_6 * gh1_33[k]
                   + pb_x[k] * gi_54[k];

        t_117[k] = f_9 * gh0_30[k]
                   - f_10 * gh1_30[k]
                   + pb_y[k] * gi_51[k];

        t_118[k] = f_5 * gh0_31[k]
                   - f_6 * gh1_31[k]
                   + pb_y[k] * gi_52[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, fi_42, fi_43, gh0_32, gh0_38, \
                         gh1_32, gh1_38, gi_53, gi_54, gi_55, gi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * gh0_32[k]
                   - f_4 * gh1_32[k]
                   + pb_y[k] * gi_53[k];

        t_120[k] = pb_y[k] * gi_54[k];

        t_121[k] = f_12 * fi_42[k]
                   + f_3 * gh0_38[k]
                   - f_4 * gh1_38[k]
                   + pb_x[k] * gi_55[k];

        t_122[k] = f_12 * fi_43[k]
                   + pb_x[k] * gi_61[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_y, pb_z, fi_26, gh0_34, gh0_35, \
                         gh0_36, gh1_34, gh1_35, gh1_36, gi_56, gi_57, \
                         gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_1 * gh0_34[k]
                   - f_2 * gh1_34[k]
                   + pb_y[k] * gi_56[k];

        t_124[k] = f_12 * fi_26[k]
                   + pb_z[k] * gi_56[k];

        t_125[k] = f_9 * gh0_35[k]
                   - f_10 * gh1_35[k]
                   + pb_y[k] * gi_57[k];

        t_126[k] = f_7 * gh0_36[k]
                   - f_8 * gh1_36[k]
                   + pb_y[k] * gi_58[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pb_y, dk0_2, dk1_95, fk_53, gh0_37, \
                         gh0_38, gh1_37, gh1_38, gi_59, gi_60, gi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_5 * gh0_37[k]
                   - f_6 * gh1_37[k]
                   + pb_y[k] * gi_59[k];

        t_128[k] = f_3 * gh0_38[k]
                   - f_4 * gh1_38[k]
                   + pb_y[k] * gi_60[k];

        t_129[k] = pb_y[k] * gi_61[k];

        t_130[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_95[k]
                   + pa_x[k] * fk_53[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_x, pb_y, fi_32, fi_44, fi_46, fi_47, \
                         fk_54, fk_55, fk_56, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_15 * fi_44[k]
                   + pa_x[k] * fk_54[k];

        t_132[k] = f_13 * fi_32[k]
                   + pb_y[k] * gi_62[k];

        t_133[k] = f_14 * fi_46[k]
                   + pa_x[k] * fk_55[k];

        t_134[k] = f_14 * fi_47[k]
                   + pa_x[k] * fk_56[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pa_x, fi_48, fi_49, fi_50, fi_52, \
                         fi_53, fk_57, fk_58, fk_59, fk_61, fk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * fi_48[k]
                   + pa_x[k] * fk_57[k];

        t_136[k] = f_0 * fi_49[k]
                   + pa_x[k] * fk_58[k];

        t_137[k] = f_13 * fi_50[k]
                   + pa_x[k] * fk_59[k];

        t_138[k] = f_13 * fi_52[k]
                   + pa_x[k] * fk_61[k];

        t_139[k] = f_12 * fi_53[k]
                   + pa_x[k] * fk_62[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, pa_x, pb_x, fi_56, fi_57, \
                         fk_65, fk_71, fk_73, fk_74, fk_75, gi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_12 * fi_56[k]
                   + pa_x[k] * fk_65[k];

        t_141[k] = f_11 * fi_57[k]
                   + pb_x[k] * gi_63[k];

        t_142[k] = pa_x[k] * fk_71[k];

        t_143[k] = pa_x[k] * fk_73[k];

        t_144[k] = pa_x[k] * fk_74[k];

        t_145[k] = pa_x[k] * fk_75[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, t_151, pa_x, pa_z, pb_z, fi_32, \
                         fk_41, fk_42, fk_76, fk_77, fk_78, gi_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * fk_76[k];

        t_147[k] = pa_x[k] * fk_77[k];

        t_148[k] = pa_x[k] * fk_78[k];

        t_149[k] = pa_z[k] * fk_41[k];

        t_150[k] = f_11 * fi_32[k]
                   + pb_z[k] * gi_64[k];

        t_151[k] = pa_z[k] * fk_42[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pa_x, pa_z, fi_64, fi_65, fi_66, \
                         fk_43, fk_44, fk_79, fk_80, fk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * fi_64[k]
                   + pa_x[k] * fk_79[k];

        t_153[k] = pa_z[k] * fk_43[k];

        t_154[k] = f_0 * fi_65[k]
                   + pa_x[k] * fk_80[k];

        t_155[k] = pa_z[k] * fk_44[k];

        t_156[k] = f_13 * fi_66[k]
                   + pa_x[k] * fk_81[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, pa_x, pa_z, fi_67, fk_45, \
                         fk_82, fk_84, fk_85, fk_86, fk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pa_z[k] * fk_45[k];

        t_158[k] = f_12 * fi_67[k]
                   + pa_x[k] * fk_82[k];

        t_159[k] = pa_x[k] * fk_84[k];

        t_160[k] = pa_x[k] * fk_85[k];

        t_161[k] = pa_x[k] * fk_86[k];

        t_162[k] = pa_x[k] * fk_87[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pa_y, fi_73, fk_47, \
                         fk_48, fk_88, fk_89, fk_90, fk_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * fk_88[k];

        t_164[k] = pa_x[k] * fk_89[k];

        t_165[k] = pa_x[k] * fk_90[k];

        t_166[k] = pa_y[k] * fk_47[k];

        t_167[k] = pa_y[k] * fk_48[k];

        t_168[k] = f_14 * fi_73[k]
                   + pa_x[k] * fk_91[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_x, pa_y, fi_74, fi_75, fk_49, \
                         fk_50, fk_51, fk_92, fk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pa_y[k] * fk_49[k];

        t_170[k] = f_0 * fi_74[k]
                   + pa_x[k] * fk_92[k];

        t_171[k] = pa_y[k] * fk_50[k];

        t_172[k] = f_13 * fi_75[k]
                   + pa_x[k] * fk_93[k];

        t_173[k] = pa_y[k] * fk_51[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, pa_x, pa_y, fi_76, fk_52, \
                         fk_94, fk_95, fk_96, fk_97, fk_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_12 * fi_76[k]
                   + pa_x[k] * fk_94[k];

        t_175[k] = pa_y[k] * fk_52[k];

        t_176[k] = pa_x[k] * fk_95[k];

        t_177[k] = pa_x[k] * fk_96[k];

        t_178[k] = pa_x[k] * fk_97[k];

        t_179[k] = pa_x[k] * fk_98[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pa_x, pb_z, fi_38, fi_83, fk_99, \
                         fk_100, fk_101, fk_103, gi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_x[k] * fk_99[k];

        t_181[k] = pa_x[k] * fk_100[k];

        t_182[k] = pa_x[k] * fk_101[k];

        t_183[k] = f_15 * fi_83[k]
                   + pa_x[k] * fk_103[k];

        t_184[k] = f_13 * fi_38[k]
                   + pb_z[k] * gi_65[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, pa_x, fi_85, fi_86, fi_87, fi_88, \
                         fi_89, fk_105, fk_106, fk_107, fk_108, \
                         fk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * fi_85[k]
                   + pa_x[k] * fk_105[k];

        t_186[k] = f_14 * fi_86[k]
                   + pa_x[k] * fk_106[k];

        t_187[k] = f_0 * fi_87[k]
                   + pa_x[k] * fk_107[k];

        t_188[k] = f_0 * fi_88[k]
                   + pa_x[k] * fk_108[k];

        t_189[k] = f_13 * fi_89[k]
                   + pa_x[k] * fk_109[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pa_x, pb_x, fi_91, fi_92, fi_95, \
                         fi_101, fk_111, fk_112, fk_115, fk_121, \
                         gi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_13 * fi_91[k]
                   + pa_x[k] * fk_111[k];

        t_191[k] = f_12 * fi_92[k]
                   + pa_x[k] * fk_112[k];

        t_192[k] = f_12 * fi_95[k]
                   + pa_x[k] * fk_115[k];

        t_193[k] = f_11 * fi_101[k]
                   + pb_x[k] * gi_66[k];

        t_194[k] = pa_x[k] * fk_121[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, t_200, pa_x, fk_122, fk_123, \
                         fk_124, fk_125, fk_126, fk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_x[k] * fk_122[k];

        t_196[k] = pa_x[k] * fk_123[k];

        t_197[k] = pa_x[k] * fk_124[k];

        t_198[k] = pa_x[k] * fk_125[k];

        t_199[k] = pa_x[k] * fk_126[k];

        t_200[k] = pa_x[k] * fk_128[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pb_x, pb_y, fi_44, gh0_39, gh0_40, \
                         gh0_41, gh1_39, gh1_40, gh1_41, gi_67, gi_68, \
                         gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_1 * gh0_39[k]
                   - f_2 * gh1_39[k]
                   + pb_x[k] * gi_67[k];

        t_202[k] = f_0 * fi_44[k]
                   + pb_y[k] * gi_67[k];

        t_203[k] = f_9 * gh0_40[k]
                   - f_10 * gh1_40[k]
                   + pb_x[k] * gi_68[k];

        t_204[k] = f_9 * gh0_41[k]
                   - f_10 * gh1_41[k]
                   + pb_x[k] * gi_69[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pb_x, gh0_42, gh0_43, gh0_44, gh1_42, gh1_43, \
                         gh1_44, gi_70, gi_71, gi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * gh0_42[k]
                   - f_8 * gh1_42[k]
                   + pb_x[k] * gi_70[k];

        t_206[k] = f_7 * gh0_43[k]
                   - f_8 * gh1_43[k]
                   + pb_x[k] * gi_71[k];

        t_207[k] = f_5 * gh0_44[k]
                   - f_6 * gh1_44[k]
                   + pb_x[k] * gi_72[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pb_x, gh0_45, gh0_46, gh0_47, gh1_45, gh1_46, \
                         gh1_47, gi_73, gi_74, gi_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_5 * gh0_45[k]
                   - f_6 * gh1_45[k]
                   + pb_x[k] * gi_73[k];

        t_209[k] = f_5 * gh0_46[k]
                   - f_6 * gh1_46[k]
                   + pb_x[k] * gi_74[k];

        t_210[k] = f_3 * gh0_47[k]
                   - f_4 * gh1_47[k]
                   + pb_x[k] * gi_75[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pb_x, gh0_49, gh0_50, gh0_51, gh1_49, \
                         gh1_50, gh1_51, gi_76, gi_77, gi_78, gi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_3 * gh0_49[k]
                   - f_4 * gh1_49[k]
                   + pb_x[k] * gi_76[k];

        t_212[k] = f_3 * gh0_50[k]
                   - f_4 * gh1_50[k]
                   + pb_x[k] * gi_77[k];

        t_213[k] = f_3 * gh0_51[k]
                   - f_4 * gh1_51[k]
                   + pb_x[k] * gi_78[k];

        t_214[k] = pb_x[k] * gi_79[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pb_x, pb_y, fi_57, gh0_47, gh1_47, \
                         gi_79, gi_81, gi_82, gi_83, gi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pb_x[k] * gi_81[k];

        t_216[k] = pb_x[k] * gi_82[k];

        t_217[k] = pb_x[k] * gi_83[k];

        t_218[k] = pb_x[k] * gi_84[k];

        t_219[k] = f_0 * fi_57[k]
                   + f_1 * gh0_47[k]
                   - f_2 * gh1_47[k]
                   + pb_y[k] * gi_79[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_z, gh0_47, gh0_48, gh0_49, gh1_47, \
                         gh1_48, gh1_49, gi_79, gi_80, gi_81, gi_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pb_z[k] * gi_79[k];

        t_221[k] = f_3 * gh0_47[k]
                   - f_4 * gh1_47[k]
                   + pb_z[k] * gi_80[k];

        t_222[k] = f_5 * gh0_48[k]
                   - f_6 * gh1_48[k]
                   + pb_z[k] * gi_81[k];

        t_223[k] = f_7 * gh0_49[k]
                   - f_8 * gh1_49[k]
                   + pb_z[k] * gi_82[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_z, pb_y, pb_z, fi_62, fk_54, gh0_50, \
                         gh0_51, gh1_50, gh1_51, gi_83, gi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_9 * gh0_50[k]
                   - f_10 * gh1_50[k]
                   + pb_z[k] * gi_83[k];

        t_225[k] = f_0 * fi_62[k]
                   + pb_y[k] * gi_84[k];

        t_226[k] = f_1 * gh0_51[k]
                   - f_2 * gh1_51[k]
                   + pb_z[k] * gi_84[k];

        t_227[k] = pa_z[k] * fk_54[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, pa_z, pb_z, fi_44, fi_45, fi_47, \
                         fk_55, fk_56, fk_57, fk_58, gi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_11 * fi_44[k]
                   + pb_z[k] * gi_85[k];

        t_229[k] = pa_z[k] * fk_55[k];

        t_230[k] = f_12 * fi_45[k]
                   + pa_z[k] * fk_56[k];

        t_231[k] = pa_z[k] * fk_57[k];

        t_232[k] = f_13 * fi_47[k]
                   + pa_z[k] * fk_58[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, t_237, pa_z, fi_49, fi_52, fk_59, fk_61, \
                         fk_62, fk_65, fk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_z[k] * fk_59[k];

        t_234[k] = f_0 * fi_49[k]
                   + pa_z[k] * fk_61[k];

        t_235[k] = pa_z[k] * fk_62[k];

        t_236[k] = f_14 * fi_52[k]
                   + pa_z[k] * fk_65[k];

        t_237[k] = pa_z[k] * fk_71[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_z, pb_z, fi_57, fi_58, fi_59, fi_60, \
                         fk_73, fk_74, fk_75, gi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_11 * fi_57[k]
                   + pb_z[k] * gi_86[k];

        t_239[k] = f_12 * fi_58[k]
                   + pa_z[k] * fk_73[k];

        t_240[k] = f_13 * fi_59[k]
                   + pa_z[k] * fk_74[k];

        t_241[k] = f_0 * fi_60[k]
                   + pa_z[k] * fk_75[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_z, pb_x, pb_y, fi_61, fi_62, fi_72, \
                         fk_76, fk_78, gh0_52, gh1_52, gi_87, gi_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_14 * fi_61[k]
                   + pa_z[k] * fk_76[k];

        t_243[k] = f_13 * fi_72[k]
                   + pb_y[k] * gi_87[k];

        t_244[k] = f_15 * fi_62[k]
                   + pa_z[k] * fk_78[k];

        t_245[k] = f_1 * gh0_52[k]
                   - f_2 * gh1_52[k]
                   + pb_x[k] * gi_88[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pb_x, pb_z, fi_63, gh0_53, gh0_54, gh1_53, \
                         gh1_54, gi_88, gi_89, gi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * fi_63[k]
                   + pb_z[k] * gi_88[k];

        t_247[k] = f_9 * gh0_53[k]
                   - f_10 * gh1_53[k]
                   + pb_x[k] * gi_89[k];

        t_248[k] = f_9 * gh0_54[k]
                   - f_10 * gh1_54[k]
                   + pb_x[k] * gi_90[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_x, gh0_55, gh0_56, gh0_57, gh1_55, gh1_56, \
                         gh1_57, gi_91, gi_92, gi_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_7 * gh0_55[k]
                   - f_8 * gh1_55[k]
                   + pb_x[k] * gi_91[k];

        t_250[k] = f_7 * gh0_56[k]
                   - f_8 * gh1_56[k]
                   + pb_x[k] * gi_92[k];

        t_251[k] = f_5 * gh0_57[k]
                   - f_6 * gh1_57[k]
                   + pb_x[k] * gi_93[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pb_x, gh0_58, gh0_59, gh0_60, gh1_58, gh1_59, \
                         gh1_60, gi_94, gi_95, gi_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_5 * gh0_58[k]
                   - f_6 * gh1_58[k]
                   + pb_x[k] * gi_94[k];

        t_253[k] = f_5 * gh0_59[k]
                   - f_6 * gh1_59[k]
                   + pb_x[k] * gi_95[k];

        t_254[k] = f_3 * gh0_60[k]
                   - f_4 * gh1_60[k]
                   + pb_x[k] * gi_96[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pb_x, gh0_61, gh0_62, gh0_64, gh1_61, \
                         gh1_62, gh1_64, gi_97, gi_98, gi_99, gi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_3 * gh0_61[k]
                   - f_4 * gh1_61[k]
                   + pb_x[k] * gi_97[k];

        t_256[k] = f_3 * gh0_62[k]
                   - f_4 * gh1_62[k]
                   + pb_x[k] * gi_98[k];

        t_257[k] = f_3 * gh0_64[k]
                   - f_4 * gh1_64[k]
                   + pb_x[k] * gi_99[k];

        t_258[k] = pb_x[k] * gi_100[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, pa_z, pb_x, dk0_1, dk1_54, fk_83, \
                         gi_101, gi_102, gi_103, gi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = pb_x[k] * gi_101[k];

        t_260[k] = pb_x[k] * gi_102[k];

        t_261[k] = pb_x[k] * gi_103[k];

        t_262[k] = pb_x[k] * gi_105[k];

        t_263[k] = f_16 * dk0_1[k]
                   - f_17 * dk1_54[k]
                   + pa_z[k] * fk_83[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pb_y, pb_z, fi_68, fi_78, fi_79, gh0_61, gh0_62, \
                         gh1_61, gh1_62, gi_100, gi_101, gi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_12 * fi_68[k]
                   + pb_z[k] * gi_100[k];

        t_265[k] = f_12 * fi_78[k]
                   + f_9 * gh0_61[k]
                   - f_10 * gh1_61[k]
                   + pb_y[k] * gi_101[k];

        t_266[k] = f_12 * fi_79[k]
                   + f_7 * gh0_62[k]
                   - f_8 * gh1_62[k]
                   + pb_y[k] * gi_102[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pb_y, fi_80, fi_81, fi_82, gh0_63, gh0_64, \
                         gh1_63, gh1_64, gi_103, gi_104, gi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * fi_80[k]
                   + f_5 * gh0_63[k]
                   - f_6 * gh1_63[k]
                   + pb_y[k] * gi_103[k];

        t_268[k] = f_12 * fi_81[k]
                   + f_3 * gh0_64[k]
                   - f_4 * gh1_64[k]
                   + pb_y[k] * gi_104[k];

        t_269[k] = f_12 * fi_82[k]
                   + pb_y[k] * gi_105[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, pa_y, dk0_2, dk1_95, fi_84, \
                         fk_102, fk_103, fk_104, fk_105, fk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_95[k]
                   + pa_y[k] * fk_102[k];

        t_271[k] = pa_y[k] * fk_103[k];

        t_272[k] = pa_y[k] * fk_104[k];

        t_273[k] = f_12 * fi_84[k]
                   + pa_y[k] * fk_105[k];

        t_274[k] = pa_y[k] * fk_106[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, t_280, pa_y, fi_85, fi_87, fi_89, \
                         fk_107, fk_108, fk_109, fk_111, fk_112, \
                         fk_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_13 * fi_85[k]
                   + pa_y[k] * fk_107[k];

        t_276[k] = pa_y[k] * fk_108[k];

        t_277[k] = f_0 * fi_87[k]
                   + pa_y[k] * fk_109[k];

        t_278[k] = pa_y[k] * fk_111[k];

        t_279[k] = f_14 * fi_89[k]
                   + pa_y[k] * fk_112[k];

        t_280[k] = pa_y[k] * fk_115[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pa_y, pb_z, fi_77, fi_96, fi_97, fi_98, \
                         fk_121, fk_123, fk_124, gi_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_15 * fi_96[k]
                   + pa_y[k] * fk_121[k];

        t_282[k] = f_13 * fi_77[k]
                   + pb_z[k] * gi_106[k];

        t_283[k] = f_14 * fi_97[k]
                   + pa_y[k] * fk_123[k];

        t_284[k] = f_0 * fi_98[k]
                   + pa_y[k] * fk_124[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pa_y, pb_y, fi_99, fi_100, fi_101, \
                         fk_125, fk_126, fk_128, gi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_13 * fi_99[k]
                   + pa_y[k] * fk_125[k];

        t_286[k] = f_12 * fi_100[k]
                   + pa_y[k] * fk_126[k];

        t_287[k] = f_11 * fi_101[k]
                   + pb_y[k] * gi_107[k];

        t_288[k] = pa_y[k] * fk_128[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pb_x, pb_z, fi_83, gh0_65, gh0_66, \
                         gh0_67, gh1_65, gh1_66, gh1_67, gi_108, gi_109, \
                         gi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_1 * gh0_65[k]
                   - f_2 * gh1_65[k]
                   + pb_x[k] * gi_108[k];

        t_290[k] = f_0 * fi_83[k]
                   + pb_z[k] * gi_108[k];

        t_291[k] = f_9 * gh0_66[k]
                   - f_10 * gh1_66[k]
                   + pb_x[k] * gi_109[k];

        t_292[k] = f_9 * gh0_67[k]
                   - f_10 * gh1_67[k]
                   + pb_x[k] * gi_110[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pb_x, gh0_68, gh0_69, gh0_70, gh1_68, gh1_69, \
                         gh1_70, gi_111, gi_112, gi_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_7 * gh0_68[k]
                   - f_8 * gh1_68[k]
                   + pb_x[k] * gi_111[k];

        t_294[k] = f_7 * gh0_69[k]
                   - f_8 * gh1_69[k]
                   + pb_x[k] * gi_112[k];

        t_295[k] = f_5 * gh0_70[k]
                   - f_6 * gh1_70[k]
                   + pb_x[k] * gi_113[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, pb_x, gh0_71, gh0_72, gh0_73, gh1_71, gh1_72, \
                         gh1_73, gi_114, gi_115, gi_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_5 * gh0_71[k]
                   - f_6 * gh1_71[k]
                   + pb_x[k] * gi_114[k];

        t_297[k] = f_5 * gh0_72[k]
                   - f_6 * gh1_72[k]
                   + pb_x[k] * gi_115[k];

        t_298[k] = f_3 * gh0_73[k]
                   - f_4 * gh1_73[k]
                   + pb_x[k] * gi_116[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pb_x, gh0_74, gh0_75, gh0_77, gh1_74, \
                         gh1_75, gh1_77, gi_117, gi_118, gi_119, \
                         gi_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_3 * gh0_74[k]
                   - f_4 * gh1_74[k]
                   + pb_x[k] * gi_117[k];

        t_300[k] = f_3 * gh0_75[k]
                   - f_4 * gh1_75[k]
                   + pb_x[k] * gi_118[k];

        t_301[k] = f_3 * gh0_77[k]
                   - f_4 * gh1_77[k]
                   + pb_x[k] * gi_119[k];

        t_302[k] = pb_x[k] * gi_120[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, pb_x, pb_y, gh0_73, gh1_73, \
                         gi_120, gi_121, gi_122, gi_123, gi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pb_x[k] * gi_121[k];

        t_304[k] = pb_x[k] * gi_122[k];

        t_305[k] = pb_x[k] * gi_123[k];

        t_306[k] = pb_x[k] * gi_125[k];

        t_307[k] = f_1 * gh0_73[k]
                   - f_2 * gh1_73[k]
                   + pb_y[k] * gi_120[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pb_y, pb_z, fi_96, gh0_74, gh0_75, gh1_74, \
                         gh1_75, gi_120, gi_121, gi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * fi_96[k]
                   + pb_z[k] * gi_120[k];

        t_309[k] = f_9 * gh0_74[k]
                   - f_10 * gh1_74[k]
                   + pb_y[k] * gi_121[k];

        t_310[k] = f_7 * gh0_75[k]
                   - f_8 * gh1_75[k]
                   + pb_y[k] * gi_122[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_y, pb_z, fi_101, gh0_76, gh0_77, \
                         gh1_76, gh1_77, gi_123, gi_124, gi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_5 * gh0_76[k]
                   - f_6 * gh1_76[k]
                   + pb_y[k] * gi_123[k];

        t_312[k] = f_3 * gh0_77[k]
                   - f_4 * gh1_77[k]
                   + pb_y[k] * gi_124[k];

        t_313[k] = pb_y[k] * gi_125[k];

        t_314[k] = f_0 * fi_101[k]
                   + f_1 * gh0_77[k]
                   - f_2 * gh1_77[k]
                   + pb_z[k] * gi_125[k];
    }
}

auto
compute_prim_gk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

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
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
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
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
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
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);

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
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
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
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_43 = buffer.data(gh0 + 43);
    const auto *gh0_44 = buffer.data(gh0 + 44);
    const auto *gh0_45 = buffer.data(gh0 + 45);
    const auto *gh0_46 = buffer.data(gh0 + 46);
    const auto *gh0_47 = buffer.data(gh0 + 47);
    const auto *gh0_48 = buffer.data(gh0 + 48);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_74 = buffer.data(gh0 + 74);
    const auto *gh0_75 = buffer.data(gh0 + 75);
    const auto *gh0_76 = buffer.data(gh0 + 76);
    const auto *gh0_77 = buffer.data(gh0 + 77);
    const auto *gh0_78 = buffer.data(gh0 + 78);
    const auto *gh0_79 = buffer.data(gh0 + 79);
    const auto *gh0_80 = buffer.data(gh0 + 80);
    const auto *gh0_81 = buffer.data(gh0 + 81);
    const auto *gh0_82 = buffer.data(gh0 + 82);
    const auto *gh0_83 = buffer.data(gh0 + 83);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_1 = buffer.data(gh1 + 1);
    const auto *gh1_2 = buffer.data(gh1 + 2);
    const auto *gh1_3 = buffer.data(gh1 + 3);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_42 = buffer.data(gh1 + 42);
    const auto *gh1_43 = buffer.data(gh1 + 43);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);
    const auto *gh1_48 = buffer.data(gh1 + 48);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_74 = buffer.data(gh1 + 74);
    const auto *gh1_75 = buffer.data(gh1 + 75);
    const auto *gh1_76 = buffer.data(gh1 + 76);
    const auto *gh1_77 = buffer.data(gh1 + 77);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_87 = buffer.data(gh1 + 87);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_89 = buffer.data(gh1 + 89);
    const auto *gh1_90 = buffer.data(gh1 + 90);
    const auto *gh1_91 = buffer.data(gh1 + 91);
    const auto *gh1_92 = buffer.data(gh1 + 92);
    const auto *gh1_93 = buffer.data(gh1 + 93);
    const auto *gh1_94 = buffer.data(gh1 + 94);
    const auto *gh1_95 = buffer.data(gh1 + 95);
    const auto *gh1_96 = buffer.data(gh1 + 96);
    const auto *gh1_97 = buffer.data(gh1 + 97);
    const auto *gh1_98 = buffer.data(gh1 + 98);
    const auto *gh1_107 = buffer.data(gh1 + 107);
    const auto *gh1_109 = buffer.data(gh1 + 109);
    const auto *gh1_110 = buffer.data(gh1 + 110);
    const auto *gh1_111 = buffer.data(gh1 + 111);
    const auto *gh1_112 = buffer.data(gh1 + 112);
    const auto *gh1_113 = buffer.data(gh1 + 113);
    const auto *gh1_114 = buffer.data(gh1 + 114);
    const auto *gh1_115 = buffer.data(gh1 + 115);
    const auto *gh1_116 = buffer.data(gh1 + 116);
    const auto *gh1_117 = buffer.data(gh1 + 117);
    const auto *gh1_118 = buffer.data(gh1 + 118);
    const auto *gh1_119 = buffer.data(gh1 + 119);
    const auto *gh1_120 = buffer.data(gh1 + 120);

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
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
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
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
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
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_102 = buffer.data(gi + 102);
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
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
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
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fi_0, gh0_0, gh0_1, gh1_0, \
                         gh1_1, gi_0, gi_1, gi_2, gi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_2[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];

        t_3[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, gh0_2, gh0_3, gh0_4, gh1_2, gh1_3, \
                         gh1_4, gi_4, gi_5, gi_6, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_5[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];

        t_6[k] = f_3 * gh0_4[k]
                 - f_4 * gh1_4[k]
                 + pb_y[k] * gi_6[k];

        t_7[k] = f_7 * gh0_4[k]
                 - f_8 * gh1_4[k]
                 + pb_z[k] * gi_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, gh0_5, gh0_6, gh0_7, gh1_5, gh1_6, \
                         gh1_7, gi_8, gi_9, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * gh0_5[k]
                 - f_10 * gh1_5[k]
                 + pb_y[k] * gi_8[k];

        t_9[k] = f_5 * gh0_6[k]
                 - f_6 * gh1_6[k]
                 + pb_y[k] * gi_9[k];

        t_10[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_11[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, fi_12, fi_17, gh0_8, gh0_9, \
                         gh1_8, gh1_10, gi_12, gi_13, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * fi_12[k]
                  + pb_x[k] * gi_12[k];

        t_13[k] = f_0 * fi_17[k]
                  + pb_x[k] * gi_17[k];

        t_14[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_15[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_10[k]
                  + pb_y[k] * gi_13[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, gh0_10, gh0_11, gh0_12, gh1_11, \
                         gh1_12, gh1_13, gi_14, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_11[k]
                  + pb_y[k] * gi_14[k];

        t_17[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_12[k]
                  + pb_y[k] * gi_15[k];

        t_18[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_13[k]
                  + pb_y[k] * gi_16[k];

        t_19[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_13[k]
                  + pb_z[k] * gi_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, fi_0, fi_1, fi_3, fi_5, \
                         fk_0, fk_1, fk_3, fk_5, gi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * fk_0[k];

        t_21[k] = f_11 * fi_0[k]
                  + pb_y[k] * gi_18[k];

        t_22[k] = f_12 * fi_1[k]
                  + pa_y[k] * fk_1[k];

        t_23[k] = f_13 * fi_3[k]
                  + pa_y[k] * fk_3[k];

        t_24[k] = f_0 * fi_5[k]
                  + pa_y[k] * fk_5[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_x, fi_8, fi_12, fi_22, fk_0, \
                         fk_8, fk_12, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_14 * fi_8[k]
                  + pa_y[k] * fk_8[k];

        t_26[k] = f_13 * fi_22[k]
                  + pb_x[k] * gi_22[k];

        t_27[k] = f_15 * fi_12[k]
                  + pa_y[k] * fk_12[k];

        t_28[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_z, pb_z, fi_0, fi_2, fi_4, fi_6, fk_2, \
                         fk_4, fk_6, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * fi_0[k]
                  + pb_z[k] * gi_23[k];

        t_30[k] = f_12 * fi_2[k]
                  + pa_z[k] * fk_2[k];

        t_31[k] = f_13 * fi_4[k]
                  + pa_z[k] * fk_4[k];

        t_32[k] = f_12 * fi_6[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, fi_7, fi_9, fi_10, fi_11, fk_7, fk_9, \
                         fk_10, fk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * fi_7[k]
                  + pa_z[k] * fk_7[k];

        t_34[k] = f_12 * fi_9[k]
                  + pa_z[k] * fk_9[k];

        t_35[k] = f_13 * fi_10[k]
                  + pa_z[k] * fk_10[k];

        t_36[k] = f_14 * fi_11[k]
                  + pa_z[k] * fk_11[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, fi_13, fi_14, fi_15, fi_31, \
                         fk_13, fk_14, fk_15, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * fi_31[k]
                  + pb_x[k] * gi_31[k];

        t_38[k] = f_12 * fi_13[k]
                  + pa_z[k] * fk_13[k];

        t_39[k] = f_13 * fi_14[k]
                  + pa_z[k] * fk_14[k];

        t_40[k] = f_0 * fi_15[k]
                  + pa_z[k] * fk_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, dk0_0, dk1_0, fi_16, fi_17, \
                         fi_18, fk_16, fk_17, fk_18, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_14 * fi_16[k]
                  + pa_z[k] * fk_16[k];

        t_42[k] = f_15 * fi_17[k]
                  + pa_z[k] * fk_17[k];

        t_43[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_y[k] * fk_18[k];

        t_44[k] = f_12 * fi_18[k]
                  + pb_y[k] * gi_32[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, pb_z, fi_33, fi_34, gh0_15, gh0_17, gh0_19, \
                         gh1_25, gh1_27, gh1_29, gi_33, gi_34, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * fi_33[k]
                  + f_9 * gh0_17[k]
                  - f_10 * gh1_27[k]
                  + pb_x[k] * gi_34[k];

        t_46[k] = f_3 * gh0_15[k]
                  - f_4 * gh1_25[k]
                  + pb_z[k] * gi_33[k];

        t_47[k] = f_12 * fi_34[k]
                  + f_7 * gh0_19[k]
                  - f_8 * gh1_29[k]
                  + pb_x[k] * gi_36[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_z, fi_35, gh0_16, gh0_17, gh0_22, gh1_26, \
                         gh1_27, gh1_32, gi_35, gi_37, gi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * gh0_16[k]
                  - f_6 * gh1_26[k]
                  + pb_z[k] * gi_35[k];

        t_49[k] = f_12 * fi_35[k]
                  + f_5 * gh0_22[k]
                  - f_6 * gh1_32[k]
                  + pb_x[k] * gi_39[k];

        t_50[k] = f_3 * gh0_17[k]
                  - f_4 * gh1_27[k]
                  + pb_z[k] * gi_37[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, pb_z, fi_36, gh0_18, gh0_19, gh0_23, gh1_28, \
                         gh1_29, gh1_33, gi_38, gi_40, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_7 * gh0_18[k]
                  - f_8 * gh1_28[k]
                  + pb_z[k] * gi_38[k];

        t_52[k] = f_12 * fi_36[k]
                  + f_3 * gh0_23[k]
                  - f_4 * gh1_33[k]
                  + pb_x[k] * gi_43[k];

        t_53[k] = f_3 * gh0_19[k]
                  - f_4 * gh1_29[k]
                  + pb_z[k] * gi_40[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_z, fi_37, gh0_20, gh0_21, gh1_30, gh1_31, \
                         gi_41, gi_42, gi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * gh0_20[k]
                  - f_6 * gh1_30[k]
                  + pb_z[k] * gi_41[k];

        t_55[k] = f_9 * gh0_21[k]
                  - f_10 * gh1_31[k]
                  + pb_z[k] * gi_42[k];

        t_56[k] = f_12 * fi_37[k]
                  + pb_x[k] * gi_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pb_z, dk0_1, dk1_1, fk_20, gh0_23, gh0_24, \
                         gh1_33, gh1_34, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_16 * dk0_1[k]
                  - f_17 * dk1_1[k]
                  + pa_x[k] * fk_20[k];

        t_58[k] = f_3 * gh0_23[k]
                  - f_4 * gh1_33[k]
                  + pb_z[k] * gi_45[k];

        t_59[k] = f_5 * gh0_24[k]
                  - f_6 * gh1_34[k]
                  + pb_z[k] * gi_46[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_z, gh0_25, gh0_26, gh0_27, gh1_35, gh1_36, \
                         gh1_37, gi_47, gi_48, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_7 * gh0_25[k]
                  - f_8 * gh1_35[k]
                  + pb_z[k] * gi_47[k];

        t_61[k] = f_9 * gh0_26[k]
                  - f_10 * gh1_36[k]
                  + pb_z[k] * gi_48[k];

        t_62[k] = f_1 * gh0_27[k]
                  - f_2 * gh1_37[k]
                  + pb_z[k] * gi_49[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_z, pb_y, pb_z, dk0_0, dk1_0, fi_23, fk_19, \
                         gh0_28, gh1_40, gi_50, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_z[k] * fk_19[k];

        t_64[k] = f_12 * fi_23[k]
                  + pb_z[k] * gi_50[k];

        t_65[k] = f_3 * gh0_28[k]
                  - f_4 * gh1_40[k]
                  + pb_y[k] * gi_51[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_y, fi_40, fi_41, gh0_29, gh0_31, gh0_34, \
                         gh1_41, gh1_43, gh1_46, gi_53, gi_54, gi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * fi_40[k]
                  + f_9 * gh0_31[k]
                  - f_10 * gh1_43[k]
                  + pb_x[k] * gi_54[k];

        t_67[k] = f_5 * gh0_29[k]
                  - f_6 * gh1_41[k]
                  + pb_y[k] * gi_53[k];

        t_68[k] = f_12 * fi_41[k]
                  + f_7 * gh0_34[k]
                  - f_8 * gh1_46[k]
                  + pb_x[k] * gi_57[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_y, fi_42, gh0_30, gh0_31, gh0_35, gh1_42, \
                         gh1_43, gh1_47, gi_55, gi_56, gi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_7 * gh0_30[k]
                  - f_8 * gh1_42[k]
                  + pb_y[k] * gi_55[k];

        t_70[k] = f_3 * gh0_31[k]
                  - f_4 * gh1_43[k]
                  + pb_y[k] * gi_56[k];

        t_71[k] = f_12 * fi_42[k]
                  + f_5 * gh0_35[k]
                  - f_6 * gh1_47[k]
                  + pb_x[k] * gi_61[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, gh0_32, gh0_33, gh0_34, gh1_44, gh1_45, \
                         gh1_46, gi_58, gi_59, gi_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * gh0_32[k]
                  - f_10 * gh1_44[k]
                  + pb_y[k] * gi_58[k];

        t_73[k] = f_5 * gh0_33[k]
                  - f_6 * gh1_45[k]
                  + pb_y[k] * gi_59[k];

        t_74[k] = f_3 * gh0_34[k]
                  - f_4 * gh1_46[k]
                  + pb_y[k] * gi_60[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, pb_y, fi_43, fi_44, gh0_36, gh0_40, gh1_48, \
                         gh1_52, gi_62, gi_63, gi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_12 * fi_43[k]
                  + f_3 * gh0_40[k]
                  - f_4 * gh1_52[k]
                  + pb_x[k] * gi_62[k];

        t_76[k] = f_12 * fi_44[k]
                  + pb_x[k] * gi_68[k];

        t_77[k] = f_1 * gh0_36[k]
                  - f_2 * gh1_48[k]
                  + pb_y[k] * gi_63[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_y, gh0_37, gh0_38, gh0_39, gh1_49, gh1_50, \
                         gh1_51, gi_64, gi_65, gi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * gh0_37[k]
                  - f_10 * gh1_49[k]
                  + pb_y[k] * gi_64[k];

        t_79[k] = f_7 * gh0_38[k]
                  - f_8 * gh1_50[k]
                  + pb_y[k] * gi_65[k];

        t_80[k] = f_5 * gh0_39[k]
                  - f_6 * gh1_51[k]
                  + pb_y[k] * gi_66[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_x, pb_y, dk0_2, dk1_2, fi_32, fi_45, \
                         fk_21, fk_22, gh0_40, gh1_52, gi_67, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * gh0_40[k]
                  - f_4 * gh1_52[k]
                  + pb_y[k] * gi_67[k];

        t_82[k] = f_16 * dk0_2[k]
                  - f_17 * dk1_2[k]
                  + pa_x[k] * fk_21[k];

        t_83[k] = f_15 * fi_45[k]
                  + pa_x[k] * fk_22[k];

        t_84[k] = f_13 * fi_32[k]
                  + pb_y[k] * gi_69[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_x, fi_47, fi_49, fi_52, fi_56, fk_23, \
                         fk_25, fk_27, fk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_14 * fi_47[k]
                  + pa_x[k] * fk_23[k];

        t_86[k] = f_0 * fi_49[k]
                  + pa_x[k] * fk_25[k];

        t_87[k] = f_13 * fi_52[k]
                  + pa_x[k] * fk_27[k];

        t_88[k] = f_12 * fi_56[k]
                  + pa_x[k] * fk_30[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pb_x, pb_z, fi_38, fi_57, fi_84, fk_34, \
                         fk_42, gi_73, gi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_11 * fi_57[k]
                  + pb_x[k] * gi_73[k];

        t_90[k] = pa_x[k] * fk_34[k];

        t_91[k] = f_15 * fi_84[k]
                  + pa_x[k] * fk_42[k];

        t_92[k] = f_13 * fi_38[k]
                  + pb_z[k] * gi_74[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_x, fi_88, fi_91, fi_95, fi_96, fk_44, \
                         fk_46, fk_49, fk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_14 * fi_88[k]
                  + pa_x[k] * fk_44[k];

        t_94[k] = f_0 * fi_91[k]
                  + pa_x[k] * fk_46[k];

        t_95[k] = f_13 * fi_95[k]
                  + pa_x[k] * fk_49[k];

        t_96[k] = f_12 * fi_96[k]
                  + pa_x[k] * fk_53[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pb_x, pb_y, fi_45, fi_103, fk_59, \
                         gh0_43, gh1_64, gi_79, gi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_11 * fi_103[k]
                  + pb_x[k] * gi_79[k];

        t_98[k] = pa_x[k] * fk_59[k];

        t_99[k] = f_1 * gh0_43[k]
                  - f_2 * gh1_64[k]
                  + pb_x[k] * gi_80[k];

        t_100[k] = f_0 * fi_45[k]
                   + pb_y[k] * gi_80[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pb_x, gh0_44, gh0_45, gh0_46, gh1_66, gh1_67, \
                         gh1_68, gi_81, gi_82, gi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * gh0_44[k]
                   - f_10 * gh1_66[k]
                   + pb_x[k] * gi_81[k];

        t_102[k] = f_9 * gh0_45[k]
                   - f_10 * gh1_67[k]
                   + pb_x[k] * gi_82[k];

        t_103[k] = f_7 * gh0_46[k]
                   - f_8 * gh1_68[k]
                   + pb_x[k] * gi_83[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_x, gh0_47, gh0_48, gh0_49, gh1_69, gh1_70, \
                         gh1_71, gi_84, gi_85, gi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_7 * gh0_47[k]
                   - f_8 * gh1_69[k]
                   + pb_x[k] * gi_84[k];

        t_105[k] = f_5 * gh0_48[k]
                   - f_6 * gh1_70[k]
                   + pb_x[k] * gi_85[k];

        t_106[k] = f_5 * gh0_49[k]
                   - f_6 * gh1_71[k]
                   + pb_x[k] * gi_86[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, gh0_50, gh0_51, gh0_53, gh1_72, gh1_73, \
                         gh1_75, gi_87, gi_88, gi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_5 * gh0_50[k]
                   - f_6 * gh1_72[k]
                   + pb_x[k] * gi_87[k];

        t_108[k] = f_3 * gh0_51[k]
                   - f_4 * gh1_73[k]
                   + pb_x[k] * gi_88[k];

        t_109[k] = f_3 * gh0_53[k]
                   - f_4 * gh1_75[k]
                   + pb_x[k] * gi_89[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, pb_y, fi_57, gh0_51, gh0_54, gh0_55, \
                         gh1_73, gh1_76, gh1_77, gi_90, gi_91, gi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * gh0_54[k]
                   - f_4 * gh1_76[k]
                   + pb_x[k] * gi_90[k];

        t_111[k] = f_3 * gh0_55[k]
                   - f_4 * gh1_77[k]
                   + pb_x[k] * gi_91[k];

        t_112[k] = f_0 * fi_57[k]
                   + f_1 * gh0_51[k]
                   - f_2 * gh1_73[k]
                   + pb_y[k] * gi_92[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_z, gh0_51, gh0_52, gh0_53, gh1_73, gh1_74, \
                         gh1_75, gi_93, gi_94, gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * gh0_51[k]
                   - f_4 * gh1_73[k]
                   + pb_z[k] * gi_93[k];

        t_114[k] = f_5 * gh0_52[k]
                   - f_6 * gh1_74[k]
                   + pb_z[k] * gi_94[k];

        t_115[k] = f_7 * gh0_53[k]
                   - f_8 * gh1_75[k]
                   + pb_z[k] * gi_95[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_z, pb_y, pb_z, fi_46, fi_63, fk_24, \
                         gh0_54, gh0_55, gh1_76, gh1_77, gi_96, gi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_9 * gh0_54[k]
                   - f_10 * gh1_76[k]
                   + pb_z[k] * gi_96[k];

        t_117[k] = f_0 * fi_63[k]
                   + pb_y[k] * gi_98[k];

        t_118[k] = f_1 * gh0_55[k]
                   - f_2 * gh1_77[k]
                   + pb_z[k] * gi_98[k];

        t_119[k] = f_12 * fi_46[k]
                   + pa_z[k] * fk_24[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pa_z, fi_48, fi_50, fi_51, fi_53, \
                         fi_54, fk_26, fk_28, fk_29, fk_31, fk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_13 * fi_48[k]
                   + pa_z[k] * fk_26[k];

        t_121[k] = f_12 * fi_50[k]
                   + pa_z[k] * fk_28[k];

        t_122[k] = f_0 * fi_51[k]
                   + pa_z[k] * fk_29[k];

        t_123[k] = f_12 * fi_53[k]
                   + pa_z[k] * fk_31[k];

        t_124[k] = f_13 * fi_54[k]
                   + pa_z[k] * fk_32[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_z, pb_z, fi_55, fi_57, fi_58, \
                         fi_59, fk_33, fk_34, fk_35, fk_36, gi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_14 * fi_55[k]
                   + pa_z[k] * fk_33[k];

        t_126[k] = pa_z[k] * fk_34[k];

        t_127[k] = f_11 * fi_57[k]
                   + pb_z[k] * gi_102[k];

        t_128[k] = f_12 * fi_58[k]
                   + pa_z[k] * fk_35[k];

        t_129[k] = f_13 * fi_59[k]
                   + pa_z[k] * fk_36[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_z, pb_y, fi_60, fi_61, fi_63, fi_73, \
                         fk_37, fk_38, fk_39, gi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * fi_60[k]
                   + pa_z[k] * fk_37[k];

        t_131[k] = f_14 * fi_61[k]
                   + pa_z[k] * fk_38[k];

        t_132[k] = f_13 * fi_73[k]
                   + pb_y[k] * gi_108[k];

        t_133[k] = f_15 * fi_63[k]
                   + pa_z[k] * fk_39[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pb_x, gh0_57, gh0_58, gh0_59, gh1_86, gh1_87, \
                         gh1_88, gi_109, gi_110, gi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_1 * gh0_57[k]
                   - f_2 * gh1_86[k]
                   + pb_x[k] * gi_109[k];

        t_135[k] = f_9 * gh0_58[k]
                   - f_10 * gh1_87[k]
                   + pb_x[k] * gi_110[k];

        t_136[k] = f_9 * gh0_59[k]
                   - f_10 * gh1_88[k]
                   + pb_x[k] * gi_111[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, gh0_60, gh0_61, gh0_62, gh1_89, gh1_90, \
                         gh1_91, gi_112, gi_113, gi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_7 * gh0_60[k]
                   - f_8 * gh1_89[k]
                   + pb_x[k] * gi_112[k];

        t_138[k] = f_7 * gh0_61[k]
                   - f_8 * gh1_90[k]
                   + pb_x[k] * gi_113[k];

        t_139[k] = f_5 * gh0_62[k]
                   - f_6 * gh1_91[k]
                   + pb_x[k] * gi_114[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, gh0_63, gh0_64, gh0_65, gh1_92, gh1_93, \
                         gh1_94, gi_115, gi_116, gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_5 * gh0_63[k]
                   - f_6 * gh1_92[k]
                   + pb_x[k] * gi_115[k];

        t_141[k] = f_5 * gh0_64[k]
                   - f_6 * gh1_93[k]
                   + pb_x[k] * gi_116[k];

        t_142[k] = f_3 * gh0_65[k]
                   - f_4 * gh1_94[k]
                   + pb_x[k] * gi_117[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, gh0_66, gh0_67, gh0_69, gh1_95, gh1_96, \
                         gh1_98, gi_118, gi_119, gi_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * gh0_66[k]
                   - f_4 * gh1_95[k]
                   + pb_x[k] * gi_118[k];

        t_144[k] = f_3 * gh0_67[k]
                   - f_4 * gh1_96[k]
                   + pb_x[k] * gi_119[k];

        t_145[k] = f_3 * gh0_69[k]
                   - f_4 * gh1_98[k]
                   + pb_x[k] * gi_120[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_z, pb_y, pb_z, dk0_1, dk1_1, fi_67, fi_79, \
                         fk_40, gh0_66, gh1_95, gi_121, gi_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_16 * dk0_1[k]
                   - f_17 * dk1_1[k]
                   + pa_z[k] * fk_40[k];

        t_147[k] = f_12 * fi_67[k]
                   + pb_z[k] * gi_121[k];

        t_148[k] = f_12 * fi_79[k]
                   + f_9 * gh0_66[k]
                   - f_10 * gh1_95[k]
                   + pb_y[k] * gi_123[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_y, fi_80, fi_81, fi_82, gh0_67, gh0_68, \
                         gh0_69, gh1_96, gh1_97, gh1_98, gi_124, gi_125, \
                         gi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_12 * fi_80[k]
                   + f_7 * gh0_67[k]
                   - f_8 * gh1_96[k]
                   + pb_y[k] * gi_124[k];

        t_150[k] = f_12 * fi_81[k]
                   + f_5 * gh0_68[k]
                   - f_6 * gh1_97[k]
                   + pb_y[k] * gi_125[k];

        t_151[k] = f_12 * fi_82[k]
                   + f_3 * gh0_69[k]
                   - f_4 * gh1_98[k]
                   + pb_y[k] * gi_126[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_y, pb_y, dk0_2, dk1_2, fi_83, fi_85, \
                         fi_87, fk_41, fk_43, fk_45, gi_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_12 * fi_83[k]
                   + pb_y[k] * gi_127[k];

        t_153[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_2[k]
                   + pa_y[k] * fk_41[k];

        t_154[k] = f_12 * fi_85[k]
                   + pa_y[k] * fk_43[k];

        t_155[k] = f_13 * fi_87[k]
                   + pa_y[k] * fk_45[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pa_y, fi_89, fi_90, fi_92, fi_93, \
                         fi_94, fk_47, fk_48, fk_50, fk_51, fk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * fi_89[k]
                   + pa_y[k] * fk_47[k];

        t_157[k] = f_12 * fi_90[k]
                   + pa_y[k] * fk_48[k];

        t_158[k] = f_14 * fi_92[k]
                   + pa_y[k] * fk_50[k];

        t_159[k] = f_13 * fi_93[k]
                   + pa_y[k] * fk_51[k];

        t_160[k] = f_12 * fi_94[k]
                   + pa_y[k] * fk_52[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pa_y, pb_z, fi_77, fi_97, fi_99, fi_100, \
                         fk_54, fk_55, fk_56, gi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_15 * fi_97[k]
                   + pa_y[k] * fk_54[k];

        t_162[k] = f_13 * fi_77[k]
                   + pb_z[k] * gi_131[k];

        t_163[k] = f_14 * fi_99[k]
                   + pa_y[k] * fk_55[k];

        t_164[k] = f_0 * fi_100[k]
                   + pa_y[k] * fk_56[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pa_y, pb_y, fi_101, fi_102, fi_103, \
                         fk_57, fk_58, fk_59, gi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_13 * fi_101[k]
                   + pa_y[k] * fk_57[k];

        t_166[k] = f_12 * fi_102[k]
                   + pa_y[k] * fk_58[k];

        t_167[k] = f_11 * fi_103[k]
                   + pb_y[k] * gi_137[k];

        t_168[k] = pa_y[k] * fk_59[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pb_x, pb_z, fi_84, gh0_71, gh0_72, \
                         gh0_73, gh1_107, gh1_109, gh1_110, gi_138, gi_140, \
                         gi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_1 * gh0_71[k]
                   - f_2 * gh1_107[k]
                   + pb_x[k] * gi_138[k];

        t_170[k] = f_0 * fi_84[k]
                   + pb_z[k] * gi_138[k];

        t_171[k] = f_9 * gh0_72[k]
                   - f_10 * gh1_109[k]
                   + pb_x[k] * gi_140[k];

        t_172[k] = f_9 * gh0_73[k]
                   - f_10 * gh1_110[k]
                   + pb_x[k] * gi_141[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pb_x, gh0_74, gh0_75, gh0_76, gh1_111, gh1_112, \
                         gh1_113, gi_142, gi_143, gi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_7 * gh0_74[k]
                   - f_8 * gh1_111[k]
                   + pb_x[k] * gi_142[k];

        t_174[k] = f_7 * gh0_75[k]
                   - f_8 * gh1_112[k]
                   + pb_x[k] * gi_143[k];

        t_175[k] = f_5 * gh0_76[k]
                   - f_6 * gh1_113[k]
                   + pb_x[k] * gi_144[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pb_x, gh0_77, gh0_78, gh0_79, gh1_114, gh1_115, \
                         gh1_116, gi_145, gi_146, gi_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_5 * gh0_77[k]
                   - f_6 * gh1_114[k]
                   + pb_x[k] * gi_145[k];

        t_177[k] = f_5 * gh0_78[k]
                   - f_6 * gh1_115[k]
                   + pb_x[k] * gi_146[k];

        t_178[k] = f_3 * gh0_79[k]
                   - f_4 * gh1_116[k]
                   + pb_x[k] * gi_147[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pb_x, gh0_80, gh0_81, gh0_83, gh1_117, gh1_118, \
                         gh1_120, gi_148, gi_149, gi_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_3 * gh0_80[k]
                   - f_4 * gh1_117[k]
                   + pb_x[k] * gi_148[k];

        t_180[k] = f_3 * gh0_81[k]
                   - f_4 * gh1_118[k]
                   + pb_x[k] * gi_149[k];

        t_181[k] = f_3 * gh0_83[k]
                   - f_4 * gh1_120[k]
                   + pb_x[k] * gi_150[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_y, pb_z, fi_97, gh0_79, gh0_80, \
                         gh0_81, gh1_116, gh1_117, gh1_118, gi_151, gi_153, \
                         gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_1 * gh0_79[k]
                   - f_2 * gh1_116[k]
                   + pb_y[k] * gi_151[k];

        t_183[k] = f_0 * fi_97[k]
                   + pb_z[k] * gi_151[k];

        t_184[k] = f_9 * gh0_80[k]
                   - f_10 * gh1_117[k]
                   + pb_y[k] * gi_153[k];

        t_185[k] = f_7 * gh0_81[k]
                   - f_8 * gh1_118[k]
                   + pb_y[k] * gi_154[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pb_z, fi_103, gh0_82, gh0_83, gh1_119, \
                         gh1_120, gi_155, gi_156, gi_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_5 * gh0_82[k]
                   - f_6 * gh1_119[k]
                   + pb_y[k] * gi_155[k];

        t_187[k] = f_3 * gh0_83[k]
                   - f_4 * gh1_120[k]
                   + pb_y[k] * gi_156[k];

        t_188[k] = f_0 * fi_103[k]
                   + f_1 * gh0_83[k]
                   - f_2 * gh1_120[k]
                   + pb_z[k] * gi_157[k];
    }
}

auto
compute_prim_gk_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_11 = 1.0 / p;
    const auto f_12 = 1.5 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
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
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
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
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_73 = buffer.data(fk + 73);
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
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_130 = buffer.data(fk + 130);

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
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_47 = buffer.data(gh0 + 47);
    const auto *gh0_48 = buffer.data(gh0 + 48);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_77 = buffer.data(gh0 + 77);
    const auto *gh0_78 = buffer.data(gh0 + 78);
    const auto *gh0_79 = buffer.data(gh0 + 79);
    const auto *gh0_80 = buffer.data(gh0 + 80);
    const auto *gh0_81 = buffer.data(gh0 + 81);
    const auto *gh0_82 = buffer.data(gh0 + 82);
    const auto *gh0_83 = buffer.data(gh0 + 83);
    const auto *gh0_84 = buffer.data(gh0 + 84);
    const auto *gh0_85 = buffer.data(gh0 + 85);
    const auto *gh0_86 = buffer.data(gh0 + 86);
    const auto *gh0_87 = buffer.data(gh0 + 87);
    const auto *gh0_88 = buffer.data(gh0 + 88);
    const auto *gh0_89 = buffer.data(gh0 + 89);

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
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_63 = buffer.data(gh1 + 63);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_74 = buffer.data(gh1 + 74);
    const auto *gh1_75 = buffer.data(gh1 + 75);
    const auto *gh1_80 = buffer.data(gh1 + 80);
    const auto *gh1_81 = buffer.data(gh1 + 81);
    const auto *gh1_82 = buffer.data(gh1 + 82);
    const auto *gh1_83 = buffer.data(gh1 + 83);
    const auto *gh1_84 = buffer.data(gh1 + 84);
    const auto *gh1_85 = buffer.data(gh1 + 85);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_87 = buffer.data(gh1 + 87);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_89 = buffer.data(gh1 + 89);
    const auto *gh1_90 = buffer.data(gh1 + 90);
    const auto *gh1_91 = buffer.data(gh1 + 91);
    const auto *gh1_92 = buffer.data(gh1 + 92);

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
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
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
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
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
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
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
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, gi_0, \
                         gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, \
                         gh1_2, gh1_3, gi_3, gi_4, gi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_6[k] = pb_z[k] * gi_3[k];

        t_7[k] = pb_y[k] * gi_4[k];

        t_8[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_9[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pb_y, pb_z, gh0_4, gh0_5, gh1_4, \
                         gh1_5, gi_5, gi_6, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * gi_5[k];

        t_11[k] = f_3 * gh0_4[k]
                  - f_4 * gh1_4[k]
                  + pb_y[k] * gi_6[k];

        t_12[k] = pb_y[k] * gi_7[k];

        t_13[k] = f_7 * gh0_4[k]
                  - f_8 * gh1_4[k]
                  + pb_z[k] * gi_7[k];

        t_14[k] = f_9 * gh0_5[k]
                  - f_10 * gh1_5[k]
                  + pb_y[k] * gi_8[k];

        t_15[k] = pb_z[k] * gi_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, gh0_6, gh0_7, gh1_6, gh1_7, gi_9, \
                         gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * gh0_6[k]
                  - f_6 * gh1_6[k]
                  + pb_y[k] * gi_9[k];

        t_17[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_18[k] = pb_y[k] * gi_11[k];

        t_19[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, gh0_8, gh0_9, gh0_10, gh1_8, \
                         gh1_9, gh1_10, gi_12, gi_14, gi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_21[k] = pb_z[k] * gi_12[k];

        t_22[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_9[k]
                  + pb_y[k] * gi_14[k];

        t_23[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_10[k]
                  + pb_y[k] * gi_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_y, pb_z, fk_0, gh0_11, gh0_12, \
                         gh1_11, gh1_12, gi_16, gi_17, gi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_11[k]
                  + pb_y[k] * gi_16[k];

        t_25[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_12[k]
                  + pb_y[k] * gi_17[k];

        t_26[k] = pb_y[k] * gi_18[k];

        t_27[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_12[k]
                  + pb_z[k] * gi_18[k];

        t_28[k] = pa_y[k] * fk_0[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_y, fi_1, fi_3, fi_5, fk_3, \
                         fk_4, fk_5, fk_7, fk_8, fk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * fi_1[k]
                  + pa_y[k] * fk_3[k];

        t_30[k] = pa_y[k] * fk_4[k];

        t_31[k] = f_12 * fi_3[k]
                  + pa_y[k] * fk_5[k];

        t_32[k] = pa_y[k] * fk_7[k];

        t_33[k] = f_0 * fi_5[k]
                  + pa_y[k] * fk_8[k];

        t_34[k] = pa_y[k] * fk_11[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, fi_9, fi_14, fi_16, fi_17, fk_12, \
                         fk_16, fk_17, fk_19, fk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * fi_9[k]
                  + pa_y[k] * fk_12[k];

        t_36[k] = pa_y[k] * fk_16[k];

        t_37[k] = f_14 * fi_14[k]
                  + pa_y[k] * fk_17[k];

        t_38[k] = f_13 * fi_16[k]
                  + pa_y[k] * fk_19[k];

        t_39[k] = f_0 * fi_17[k]
                  + pa_y[k] * fk_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, fi_18, fi_19, fi_20, \
                         fk_0, fk_21, fk_22, fk_23, gi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_12 * fi_18[k]
                  + pa_y[k] * fk_21[k];

        t_41[k] = f_11 * fi_19[k]
                  + pa_y[k] * fk_22[k];

        t_42[k] = f_15 * fi_20[k]
                  + pb_y[k] * gi_21[k];

        t_43[k] = pa_y[k] * fk_23[k];

        t_44[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_z, pb_y, pb_z, fi_0, fi_2, fk_3, \
                         fk_4, fk_5, gi_22, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * fi_0[k]
                  + pb_z[k] * gi_22[k];

        t_46[k] = pa_z[k] * fk_3[k];

        t_47[k] = f_11 * fi_2[k]
                  + pa_z[k] * fk_4[k];

        t_48[k] = pa_z[k] * fk_5[k];

        t_49[k] = pb_y[k] * gi_23[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_y, fi_4, fi_6, fi_8, fk_7, \
                         fk_8, fk_10, fk_11, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_12 * fi_4[k]
                  + pa_z[k] * fk_7[k];

        t_51[k] = pa_z[k] * fk_8[k];

        t_52[k] = f_11 * fi_6[k]
                  + pa_z[k] * fk_10[k];

        t_53[k] = pb_y[k] * gi_24[k];

        t_54[k] = f_0 * fi_8[k]
                  + pa_z[k] * fk_11[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_z, pb_y, fi_10, fi_11, fi_13, fk_12, \
                         fk_14, fk_15, fk_16, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * fk_12[k];

        t_56[k] = f_11 * fi_10[k]
                  + pa_z[k] * fk_14[k];

        t_57[k] = f_12 * fi_11[k]
                  + pa_z[k] * fk_15[k];

        t_58[k] = pb_y[k] * gi_25[k];

        t_59[k] = f_13 * fi_13[k]
                  + pa_z[k] * fk_16[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_z, pb_z, fi_14, fi_15, fi_16, fi_17, \
                         fk_17, fk_19, fk_20, fk_21, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_z[k] * fk_17[k];

        t_61[k] = f_15 * fi_14[k]
                  + pb_z[k] * gi_26[k];

        t_62[k] = f_11 * fi_15[k]
                  + pa_z[k] * fk_19[k];

        t_63[k] = f_12 * fi_16[k]
                  + pa_z[k] * fk_20[k];

        t_64[k] = f_0 * fi_17[k]
                  + pa_z[k] * fk_21[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_y, pa_z, pb_y, dk0_0, dk1_0, fi_18, fi_20, \
                         fk_22, fk_23, fk_24, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_13 * fi_18[k]
                  + pa_z[k] * fk_22[k];

        t_66[k] = pb_y[k] * gi_31[k];

        t_67[k] = f_14 * fi_20[k]
                  + pa_z[k] * fk_23[k];

        t_68[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_y[k] * fk_24[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_z, fi_32, gh0_13, gh0_15, gh1_15, gh1_17, \
                         gi_32, gi_33, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_z[k] * gi_32[k];

        t_70[k] = f_11 * fi_32[k]
                  + f_9 * gh0_15[k]
                  - f_10 * gh1_17[k]
                  + pb_x[k] * gi_34[k];

        t_71[k] = f_3 * gh0_13[k]
                  - f_4 * gh1_15[k]
                  + pb_z[k] * gi_33[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, pb_z, fi_33, gh0_14, gh0_17, gh1_16, gh1_19, \
                         gi_34, gi_35, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * fi_33[k]
                  + f_7 * gh0_17[k]
                  - f_8 * gh1_19[k]
                  + pb_x[k] * gi_36[k];

        t_73[k] = pb_z[k] * gi_34[k];

        t_74[k] = f_5 * gh0_14[k]
                  - f_6 * gh1_16[k]
                  + pb_z[k] * gi_35[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, pb_z, fi_34, gh0_15, gh0_20, gh1_17, gh1_22, \
                         gi_36, gi_37, gi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_11 * fi_34[k]
                  + f_5 * gh0_20[k]
                  - f_6 * gh1_22[k]
                  + pb_x[k] * gi_39[k];

        t_76[k] = pb_z[k] * gi_36[k];

        t_77[k] = f_3 * gh0_15[k]
                  - f_4 * gh1_17[k]
                  + pb_z[k] * gi_37[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_x, pb_z, fi_35, gh0_16, gh0_21, gh1_18, gh1_23, \
                         gi_38, gi_39, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_7 * gh0_16[k]
                  - f_8 * gh1_18[k]
                  + pb_z[k] * gi_38[k];

        t_79[k] = f_11 * fi_35[k]
                  + f_3 * gh0_21[k]
                  - f_4 * gh1_23[k]
                  + pb_x[k] * gi_43[k];

        t_80[k] = pb_z[k] * gi_39[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_z, gh0_17, gh0_18, gh0_19, gh1_19, gh1_20, \
                         gh1_21, gi_40, gi_41, gi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * gh0_17[k]
                  - f_4 * gh1_19[k]
                  + pb_z[k] * gi_40[k];

        t_82[k] = f_5 * gh0_18[k]
                  - f_6 * gh1_20[k]
                  + pb_z[k] * gi_41[k];

        t_83[k] = f_9 * gh0_19[k]
                  - f_10 * gh1_21[k]
                  + pb_z[k] * gi_42[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pb_x, pb_z, dk0_1, dk1_1, fi_36, fk_47, \
                         gh0_21, gh1_23, gi_44, gi_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * fi_36[k]
                  + pb_x[k] * gi_44[k];

        t_85[k] = f_16 * dk0_1[k]
                  - f_17 * dk1_1[k]
                  + pa_x[k] * fk_47[k];

        t_86[k] = pb_z[k] * gi_44[k];

        t_87[k] = f_3 * gh0_21[k]
                  - f_4 * gh1_23[k]
                  + pb_z[k] * gi_45[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pb_z, gh0_22, gh0_23, gh0_24, gh1_24, gh1_25, \
                         gh1_26, gi_46, gi_47, gi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_5 * gh0_22[k]
                  - f_6 * gh1_24[k]
                  + pb_z[k] * gi_46[k];

        t_89[k] = f_7 * gh0_23[k]
                  - f_8 * gh1_25[k]
                  + pb_z[k] * gi_47[k];

        t_90[k] = f_9 * gh0_24[k]
                  - f_10 * gh1_26[k]
                  + pb_z[k] * gi_48[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pa_y, pa_z, pb_y, pb_z, fi_23, fk_25, \
                         fk_31, fk_32, gh0_25, gh1_27, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_11 * fi_23[k]
                  + pb_y[k] * gi_49[k];

        t_92[k] = f_1 * gh0_25[k]
                  - f_2 * gh1_27[k]
                  + pb_z[k] * gi_49[k];

        t_93[k] = pa_y[k] * fk_31[k];

        t_94[k] = pa_z[k] * fk_25[k];

        t_95[k] = pa_y[k] * fk_32[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, t_102, pa_y, pa_z, fk_26, \
                         fk_27, fk_28, fk_29, fk_33, fk_34, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * fk_26[k];

        t_97[k] = pa_y[k] * fk_33[k];

        t_98[k] = pa_z[k] * fk_27[k];

        t_99[k] = pa_y[k] * fk_34[k];

        t_100[k] = pa_z[k] * fk_28[k];

        t_101[k] = pa_y[k] * fk_35[k];

        t_102[k] = pa_z[k] * fk_29[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_y, pb_z, fi_22, fi_26, fi_27, fi_28, \
                         fk_36, fk_37, fk_38, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_15 * fi_22[k]
                   + pb_z[k] * gi_50[k];

        t_104[k] = f_13 * fi_26[k]
                   + pa_y[k] * fk_36[k];

        t_105[k] = f_0 * fi_27[k]
                   + pa_y[k] * fk_37[k];

        t_106[k] = f_12 * fi_28[k]
                   + pa_y[k] * fk_38[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pa_z, pb_y, dk0_0, dk1_0, fi_29, \
                         fi_30, fk_30, fk_39, fk_40, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_11 * fi_29[k]
                   + pa_y[k] * fk_39[k];

        t_108[k] = f_15 * fi_30[k]
                   + pb_y[k] * gi_51[k];

        t_109[k] = pa_y[k] * fk_40[k];

        t_110[k] = f_16 * dk0_0[k]
                   - f_17 * dk1_0[k]
                   + pa_z[k] * fk_30[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_y, pb_z, fi_24, fi_38, gh0_26, \
                         gh0_29, gh1_28, gh1_31, gi_52, gi_53, gi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pb_y[k] * gi_52[k];

        t_112[k] = f_11 * fi_24[k]
                   + pb_z[k] * gi_52[k];

        t_113[k] = f_3 * gh0_26[k]
                   - f_4 * gh1_28[k]
                   + pb_y[k] * gi_53[k];

        t_114[k] = f_11 * fi_38[k]
                   + f_9 * gh0_29[k]
                   - f_10 * gh1_31[k]
                   + pb_x[k] * gi_55[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_x, pb_y, fi_39, gh0_27, gh0_32, gh1_29, \
                         gh1_34, gi_54, gi_55, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_5 * gh0_27[k]
                   - f_6 * gh1_29[k]
                   + pb_y[k] * gi_54[k];

        t_116[k] = pb_y[k] * gi_55[k];

        t_117[k] = f_11 * fi_39[k]
                   + f_7 * gh0_32[k]
                   - f_8 * gh1_34[k]
                   + pb_x[k] * gi_58[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_y, gh0_28, gh0_29, gh1_30, gh1_31, gi_56, \
                         gi_57, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_7 * gh0_28[k]
                   - f_8 * gh1_30[k]
                   + pb_y[k] * gi_56[k];

        t_119[k] = f_3 * gh0_29[k]
                   - f_4 * gh1_31[k]
                   + pb_y[k] * gi_57[k];

        t_120[k] = pb_y[k] * gi_58[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_x, pb_y, fi_40, gh0_30, gh0_31, gh0_33, \
                         gh1_32, gh1_33, gh1_35, gi_59, gi_60, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_11 * fi_40[k]
                   + f_5 * gh0_33[k]
                   - f_6 * gh1_35[k]
                   + pb_x[k] * gi_62[k];

        t_122[k] = f_9 * gh0_30[k]
                   - f_10 * gh1_32[k]
                   + pb_y[k] * gi_59[k];

        t_123[k] = f_5 * gh0_31[k]
                   - f_6 * gh1_33[k]
                   + pb_y[k] * gi_60[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_x, pb_y, fi_41, fi_42, gh0_32, gh0_38, \
                         gh1_34, gh1_40, gi_61, gi_62, gi_63, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_3 * gh0_32[k]
                   - f_4 * gh1_34[k]
                   + pb_y[k] * gi_61[k];

        t_125[k] = pb_y[k] * gi_62[k];

        t_126[k] = f_11 * fi_41[k]
                   + f_3 * gh0_38[k]
                   - f_4 * gh1_40[k]
                   + pb_x[k] * gi_63[k];

        t_127[k] = f_11 * fi_42[k]
                   + pb_x[k] * gi_69[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_y, pb_z, fi_25, gh0_34, gh0_35, \
                         gh0_36, gh1_36, gh1_37, gh1_38, gi_64, gi_65, \
                         gi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * gh0_34[k]
                   - f_2 * gh1_36[k]
                   + pb_y[k] * gi_64[k];

        t_129[k] = f_11 * fi_25[k]
                   + pb_z[k] * gi_64[k];

        t_130[k] = f_9 * gh0_35[k]
                   - f_10 * gh1_37[k]
                   + pb_y[k] * gi_65[k];

        t_131[k] = f_7 * gh0_36[k]
                   - f_8 * gh1_38[k]
                   + pb_y[k] * gi_66[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pb_y, dk0_2, dk1_2, fk_55, gh0_37, \
                         gh0_38, gh1_39, gh1_40, gi_67, gi_68, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_5 * gh0_37[k]
                   - f_6 * gh1_39[k]
                   + pb_y[k] * gi_67[k];

        t_133[k] = f_3 * gh0_38[k]
                   - f_4 * gh1_40[k]
                   + pb_y[k] * gi_68[k];

        t_134[k] = pb_y[k] * gi_69[k];

        t_135[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_2[k]
                   + pa_x[k] * fk_55[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pa_x, fi_43, fi_45, fi_46, fi_47, \
                         fi_49, fk_56, fk_58, fk_59, fk_60, fk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * fi_43[k]
                   + pa_x[k] * fk_56[k];

        t_137[k] = f_13 * fi_45[k]
                   + pa_x[k] * fk_58[k];

        t_138[k] = f_13 * fi_46[k]
                   + pa_x[k] * fk_59[k];

        t_139[k] = f_0 * fi_47[k]
                   + pa_x[k] * fk_60[k];

        t_140[k] = f_0 * fi_49[k]
                   + pa_x[k] * fk_62[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, fi_50, fi_53, fi_54, fi_57, fk_63, \
                         fk_66, fk_67, fk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_12 * fi_50[k]
                   + pa_x[k] * fk_63[k];

        t_142[k] = f_12 * fi_53[k]
                   + pa_x[k] * fk_66[k];

        t_143[k] = f_11 * fi_54[k]
                   + pa_x[k] * fk_67[k];

        t_144[k] = f_11 * fi_57[k]
                   + pa_x[k] * fk_71[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_x, pb_x, fi_58, fk_73, \
                         fk_75, fk_76, fk_77, fk_78, gi_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_15 * fi_58[k]
                   + pb_x[k] * gi_75[k];

        t_146[k] = pa_x[k] * fk_73[k];

        t_147[k] = pa_x[k] * fk_75[k];

        t_148[k] = pa_x[k] * fk_76[k];

        t_149[k] = pa_x[k] * fk_77[k];

        t_150[k] = pa_x[k] * fk_78[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pa_x, pa_z, pb_z, fi_31, fk_41, \
                         fk_42, fk_79, fk_80, gi_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = pa_x[k] * fk_79[k];

        t_152[k] = pa_x[k] * fk_80[k];

        t_153[k] = pa_z[k] * fk_41[k];

        t_154[k] = f_15 * fi_31[k]
                   + pb_z[k] * gi_76[k];

        t_155[k] = pa_z[k] * fk_42[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pa_x, pa_z, fi_65, fi_66, fi_67, \
                         fk_43, fk_44, fk_81, fk_82, fk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_13 * fi_65[k]
                   + pa_x[k] * fk_81[k];

        t_157[k] = pa_z[k] * fk_43[k];

        t_158[k] = f_0 * fi_66[k]
                   + pa_x[k] * fk_82[k];

        t_159[k] = pa_z[k] * fk_44[k];

        t_160[k] = f_12 * fi_67[k]
                   + pa_x[k] * fk_83[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pa_x, pa_z, fi_68, fk_45, \
                         fk_84, fk_86, fk_87, fk_88, fk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pa_z[k] * fk_45[k];

        t_162[k] = f_11 * fi_68[k]
                   + pa_x[k] * fk_84[k];

        t_163[k] = pa_x[k] * fk_86[k];

        t_164[k] = pa_x[k] * fk_87[k];

        t_165[k] = pa_x[k] * fk_88[k];

        t_166[k] = pa_x[k] * fk_89[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, t_172, pa_x, pa_y, fi_71, fk_48, \
                         fk_49, fk_90, fk_91, fk_92, fk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pa_x[k] * fk_90[k];

        t_168[k] = pa_x[k] * fk_91[k];

        t_169[k] = pa_x[k] * fk_92[k];

        t_170[k] = pa_y[k] * fk_48[k];

        t_171[k] = pa_y[k] * fk_49[k];

        t_172[k] = f_13 * fi_71[k]
                   + pa_x[k] * fk_93[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pa_x, pa_y, fi_72, fi_73, fk_50, \
                         fk_51, fk_52, fk_94, fk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pa_y[k] * fk_50[k];

        t_174[k] = f_0 * fi_72[k]
                   + pa_x[k] * fk_94[k];

        t_175[k] = pa_y[k] * fk_51[k];

        t_176[k] = f_12 * fi_73[k]
                   + pa_x[k] * fk_95[k];

        t_177[k] = pa_y[k] * fk_52[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, t_183, pa_x, pa_y, fi_74, fk_53, \
                         fk_96, fk_97, fk_98, fk_99, fk_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_11 * fi_74[k]
                   + pa_x[k] * fk_96[k];

        t_179[k] = pa_y[k] * fk_53[k];

        t_180[k] = pa_x[k] * fk_97[k];

        t_181[k] = pa_x[k] * fk_98[k];

        t_182[k] = pa_x[k] * fk_99[k];

        t_183[k] = pa_x[k] * fk_100[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pa_x, pb_z, fi_37, fi_81, fk_101, \
                         fk_102, fk_103, fk_105, gi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pa_x[k] * fk_101[k];

        t_185[k] = pa_x[k] * fk_102[k];

        t_186[k] = pa_x[k] * fk_103[k];

        t_187[k] = f_14 * fi_81[k]
                   + pa_x[k] * fk_105[k];

        t_188[k] = f_12 * fi_37[k]
                   + pb_z[k] * gi_77[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pa_x, fi_83, fi_84, fi_85, fi_87, \
                         fi_88, fk_108, fk_109, fk_110, fk_112, \
                         fk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * fi_83[k]
                   + pa_x[k] * fk_108[k];

        t_190[k] = f_13 * fi_84[k]
                   + pa_x[k] * fk_109[k];

        t_191[k] = f_0 * fi_85[k]
                   + pa_x[k] * fk_110[k];

        t_192[k] = f_0 * fi_87[k]
                   + pa_x[k] * fk_112[k];

        t_193[k] = f_12 * fi_88[k]
                   + pa_x[k] * fk_113[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, pa_x, pb_x, fi_91, fi_92, fi_95, \
                         fi_101, fk_116, fk_117, fk_121, fk_123, \
                         gi_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_12 * fi_91[k]
                   + pa_x[k] * fk_116[k];

        t_195[k] = f_11 * fi_92[k]
                   + pa_x[k] * fk_117[k];

        t_196[k] = f_11 * fi_95[k]
                   + pa_x[k] * fk_121[k];

        t_197[k] = f_15 * fi_101[k]
                   + pb_x[k] * gi_82[k];

        t_198[k] = pa_x[k] * fk_123[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, pa_x, fk_124, fk_125, \
                         fk_126, fk_127, fk_128, fk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pa_x[k] * fk_124[k];

        t_200[k] = pa_x[k] * fk_125[k];

        t_201[k] = pa_x[k] * fk_126[k];

        t_202[k] = pa_x[k] * fk_127[k];

        t_203[k] = pa_x[k] * fk_128[k];

        t_204[k] = pa_x[k] * fk_130[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pb_x, pb_z, gh0_47, gh0_48, gh0_49, \
                         gh1_49, gh1_50, gh1_51, gi_83, gi_85, gi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_1 * gh0_47[k]
                   - f_2 * gh1_49[k]
                   + pb_x[k] * gi_83[k];

        t_206[k] = pb_z[k] * gi_83[k];

        t_207[k] = f_9 * gh0_48[k]
                   - f_10 * gh1_50[k]
                   + pb_x[k] * gi_85[k];

        t_208[k] = f_9 * gh0_49[k]
                   - f_10 * gh1_51[k]
                   + pb_x[k] * gi_86[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_x, pb_z, gh0_50, gh0_51, gh0_52, \
                         gh1_52, gh1_53, gh1_54, gi_85, gi_87, gi_88, \
                         gi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_7 * gh0_50[k]
                   - f_8 * gh1_52[k]
                   + pb_x[k] * gi_87[k];

        t_210[k] = pb_z[k] * gi_85[k];

        t_211[k] = f_7 * gh0_51[k]
                   - f_8 * gh1_53[k]
                   + pb_x[k] * gi_88[k];

        t_212[k] = f_5 * gh0_52[k]
                   - f_6 * gh1_54[k]
                   + pb_x[k] * gi_89[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_x, pb_z, gh0_53, gh0_54, gh0_55, \
                         gh1_55, gh1_56, gh1_57, gi_87, gi_90, gi_91, \
                         gi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = pb_z[k] * gi_87[k];

        t_214[k] = f_5 * gh0_53[k]
                   - f_6 * gh1_55[k]
                   + pb_x[k] * gi_90[k];

        t_215[k] = f_5 * gh0_54[k]
                   - f_6 * gh1_56[k]
                   + pb_x[k] * gi_91[k];

        t_216[k] = f_3 * gh0_55[k]
                   - f_4 * gh1_57[k]
                   + pb_x[k] * gi_92[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_x, pb_z, gh0_57, gh0_58, gh0_59, \
                         gh1_59, gh1_60, gh1_61, gi_89, gi_93, gi_94, \
                         gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_z[k] * gi_89[k];

        t_218[k] = f_3 * gh0_57[k]
                   - f_4 * gh1_59[k]
                   + pb_x[k] * gi_93[k];

        t_219[k] = f_3 * gh0_58[k]
                   - f_4 * gh1_60[k]
                   + pb_x[k] * gi_94[k];

        t_220[k] = f_3 * gh0_59[k]
                   - f_4 * gh1_61[k]
                   + pb_x[k] * gi_95[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, t_226, pb_x, pb_y, fi_58, gh0_55, \
                         gh1_57, gi_96, gi_98, gi_99, gi_100, gi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pb_x[k] * gi_96[k];

        t_222[k] = pb_x[k] * gi_98[k];

        t_223[k] = pb_x[k] * gi_99[k];

        t_224[k] = pb_x[k] * gi_100[k];

        t_225[k] = pb_x[k] * gi_101[k];

        t_226[k] = f_0 * fi_58[k]
                   + f_1 * gh0_55[k]
                   - f_2 * gh1_57[k]
                   + pb_y[k] * gi_96[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pb_z, gh0_55, gh0_56, gh0_57, gh1_57, \
                         gh1_58, gh1_59, gi_96, gi_97, gi_98, gi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_z[k] * gi_96[k];

        t_228[k] = f_3 * gh0_55[k]
                   - f_4 * gh1_57[k]
                   + pb_z[k] * gi_97[k];

        t_229[k] = f_5 * gh0_56[k]
                   - f_6 * gh1_58[k]
                   + pb_z[k] * gi_98[k];

        t_230[k] = f_7 * gh0_57[k]
                   - f_8 * gh1_59[k]
                   + pb_z[k] * gi_99[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pa_z, pb_y, pb_z, fi_63, fk_56, gh0_58, \
                         gh0_59, gh1_60, gh1_61, gi_100, gi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_9 * gh0_58[k]
                   - f_10 * gh1_60[k]
                   + pb_z[k] * gi_100[k];

        t_232[k] = f_0 * fi_63[k]
                   + pb_y[k] * gi_101[k];

        t_233[k] = f_1 * gh0_59[k]
                   - f_2 * gh1_61[k]
                   + pb_z[k] * gi_101[k];

        t_234[k] = pa_z[k] * fk_56[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pa_z, pb_z, fi_43, fi_44, fi_46, \
                         fk_58, fk_59, fk_60, fk_62, gi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * fi_43[k]
                   + pb_z[k] * gi_102[k];

        t_236[k] = pa_z[k] * fk_58[k];

        t_237[k] = f_11 * fi_44[k]
                   + pa_z[k] * fk_59[k];

        t_238[k] = pa_z[k] * fk_60[k];

        t_239[k] = f_12 * fi_46[k]
                   + pa_z[k] * fk_62[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pa_z, fi_48, fi_49, fi_51, fk_63, \
                         fk_65, fk_66, fk_67, fk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_z[k] * fk_63[k];

        t_241[k] = f_11 * fi_48[k]
                   + pa_z[k] * fk_65[k];

        t_242[k] = f_0 * fi_49[k]
                   + pa_z[k] * fk_66[k];

        t_243[k] = pa_z[k] * fk_67[k];

        t_244[k] = f_11 * fi_51[k]
                   + pa_z[k] * fk_69[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, t_250, pa_z, pb_x, fi_52, fi_53, \
                         fk_70, fk_71, gi_108, gi_109, gi_110, gi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_12 * fi_52[k]
                   + pa_z[k] * fk_70[k];

        t_246[k] = f_13 * fi_53[k]
                   + pa_z[k] * fk_71[k];

        t_247[k] = pb_x[k] * gi_108[k];

        t_248[k] = pb_x[k] * gi_109[k];

        t_249[k] = pb_x[k] * gi_110[k];

        t_250[k] = pb_x[k] * gi_111[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, pa_z, pb_z, fi_58, fi_59, fi_60, \
                         fi_61, fk_73, fk_75, fk_76, fk_77, gi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_z[k] * fk_73[k];

        t_252[k] = f_15 * fi_58[k]
                   + pb_z[k] * gi_107[k];

        t_253[k] = f_11 * fi_59[k]
                   + pa_z[k] * fk_75[k];

        t_254[k] = f_12 * fi_60[k]
                   + pa_z[k] * fk_76[k];

        t_255[k] = f_0 * fi_61[k]
                   + pa_z[k] * fk_77[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_z, pb_x, pb_y, fi_62, fi_63, fi_70, \
                         fk_78, fk_80, gh0_60, gh1_63, gi_111, gi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_13 * fi_62[k]
                   + pa_z[k] * fk_78[k];

        t_257[k] = f_12 * fi_70[k]
                   + pb_y[k] * gi_111[k];

        t_258[k] = f_14 * fi_63[k]
                   + pa_z[k] * fk_80[k];

        t_259[k] = f_1 * gh0_60[k]
                   - f_2 * gh1_63[k]
                   + pb_x[k] * gi_112[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_x, pb_z, fi_64, gh0_61, gh0_62, gh1_64, \
                         gh1_65, gi_112, gi_113, gi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_11 * fi_64[k]
                   + pb_z[k] * gi_112[k];

        t_261[k] = f_9 * gh0_61[k]
                   - f_10 * gh1_64[k]
                   + pb_x[k] * gi_113[k];

        t_262[k] = f_9 * gh0_62[k]
                   - f_10 * gh1_65[k]
                   + pb_x[k] * gi_114[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_x, gh0_63, gh0_64, gh0_65, gh1_66, gh1_67, \
                         gh1_68, gi_115, gi_116, gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_7 * gh0_63[k]
                   - f_8 * gh1_66[k]
                   + pb_x[k] * gi_115[k];

        t_264[k] = f_7 * gh0_64[k]
                   - f_8 * gh1_67[k]
                   + pb_x[k] * gi_116[k];

        t_265[k] = f_5 * gh0_65[k]
                   - f_6 * gh1_68[k]
                   + pb_x[k] * gi_117[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_x, gh0_66, gh0_67, gh0_68, gh1_69, gh1_70, \
                         gh1_71, gi_118, gi_119, gi_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_5 * gh0_66[k]
                   - f_6 * gh1_69[k]
                   + pb_x[k] * gi_118[k];

        t_267[k] = f_5 * gh0_67[k]
                   - f_6 * gh1_70[k]
                   + pb_x[k] * gi_119[k];

        t_268[k] = f_3 * gh0_68[k]
                   - f_4 * gh1_71[k]
                   + pb_x[k] * gi_120[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pb_x, gh0_69, gh0_70, gh0_72, gh1_72, \
                         gh1_73, gh1_75, gi_121, gi_122, gi_123, \
                         gi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_3 * gh0_69[k]
                   - f_4 * gh1_72[k]
                   + pb_x[k] * gi_121[k];

        t_270[k] = f_3 * gh0_70[k]
                   - f_4 * gh1_73[k]
                   + pb_x[k] * gi_122[k];

        t_271[k] = f_3 * gh0_72[k]
                   - f_4 * gh1_75[k]
                   + pb_x[k] * gi_123[k];

        t_272[k] = pb_x[k] * gi_124[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, pa_z, pb_x, dk0_1, dk1_1, fk_85, \
                         gi_125, gi_126, gi_127, gi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = pb_x[k] * gi_125[k];

        t_274[k] = pb_x[k] * gi_126[k];

        t_275[k] = pb_x[k] * gi_127[k];

        t_276[k] = pb_x[k] * gi_129[k];

        t_277[k] = f_16 * dk0_1[k]
                   - f_17 * dk1_1[k]
                   + pa_z[k] * fk_85[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pb_y, pb_z, fi_69, fi_76, fi_77, gh0_69, gh0_70, \
                         gh1_72, gh1_73, gi_124, gi_125, gi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_11 * fi_69[k]
                   + pb_z[k] * gi_124[k];

        t_279[k] = f_11 * fi_76[k]
                   + f_9 * gh0_69[k]
                   - f_10 * gh1_72[k]
                   + pb_y[k] * gi_125[k];

        t_280[k] = f_11 * fi_77[k]
                   + f_7 * gh0_70[k]
                   - f_8 * gh1_73[k]
                   + pb_y[k] * gi_126[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pb_y, fi_78, fi_79, fi_80, gh0_71, gh0_72, \
                         gh1_74, gh1_75, gi_127, gi_128, gi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_11 * fi_78[k]
                   + f_5 * gh0_71[k]
                   - f_6 * gh1_74[k]
                   + pb_y[k] * gi_127[k];

        t_282[k] = f_11 * fi_79[k]
                   + f_3 * gh0_72[k]
                   - f_4 * gh1_75[k]
                   + pb_y[k] * gi_128[k];

        t_283[k] = f_11 * fi_80[k]
                   + pb_y[k] * gi_129[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_y, dk0_2, dk1_2, fi_82, fk_104, \
                         fk_105, fk_107, fk_108, fk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_2[k]
                   + pa_y[k] * fk_104[k];

        t_285[k] = pa_y[k] * fk_105[k];

        t_286[k] = pa_y[k] * fk_107[k];

        t_287[k] = f_11 * fi_82[k]
                   + pa_y[k] * fk_108[k];

        t_288[k] = pa_y[k] * fk_109[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, pa_y, fi_83, fi_85, fi_86, fk_110, \
                         fk_112, fk_113, fk_114, fk_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_12 * fi_83[k]
                   + pa_y[k] * fk_110[k];

        t_290[k] = pa_y[k] * fk_112[k];

        t_291[k] = f_0 * fi_85[k]
                   + pa_y[k] * fk_113[k];

        t_292[k] = f_11 * fi_86[k]
                   + pa_y[k] * fk_114[k];

        t_293[k] = pa_y[k] * fk_116[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pa_y, pb_x, fi_88, fi_89, fi_90, \
                         fk_117, fk_118, fk_119, fk_121, gi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_13 * fi_88[k]
                   + pa_y[k] * fk_117[k];

        t_295[k] = f_12 * fi_89[k]
                   + pa_y[k] * fk_118[k];

        t_296[k] = f_11 * fi_90[k]
                   + pa_y[k] * fk_119[k];

        t_297[k] = pa_y[k] * fk_121[k];

        t_298[k] = pb_x[k] * gi_134[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pa_y, pb_x, pb_z, fi_75, fi_96, \
                         fk_123, gi_134, gi_135, gi_136, gi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pb_x[k] * gi_135[k];

        t_300[k] = pb_x[k] * gi_136[k];

        t_301[k] = pb_x[k] * gi_137[k];

        t_302[k] = f_14 * fi_96[k]
                   + pa_y[k] * fk_123[k];

        t_303[k] = f_12 * fi_75[k]
                   + pb_z[k] * gi_134[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_y, fi_97, fi_98, fi_99, fi_100, \
                         fk_125, fk_126, fk_127, fk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_13 * fi_97[k]
                   + pa_y[k] * fk_125[k];

        t_305[k] = f_0 * fi_98[k]
                   + pa_y[k] * fk_126[k];

        t_306[k] = f_12 * fi_99[k]
                   + pa_y[k] * fk_127[k];

        t_307[k] = f_11 * fi_100[k]
                   + pa_y[k] * fk_128[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, pa_y, pb_x, pb_y, pb_z, fi_81, \
                         fi_101, fk_130, gh0_77, gh1_80, gi_139, \
                         gi_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_15 * fi_101[k]
                   + pb_y[k] * gi_139[k];

        t_309[k] = pa_y[k] * fk_130[k];

        t_310[k] = f_1 * gh0_77[k]
                   - f_2 * gh1_80[k]
                   + pb_x[k] * gi_140[k];

        t_311[k] = pb_y[k] * gi_140[k];

        t_312[k] = f_0 * fi_81[k]
                   + pb_z[k] * gi_140[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pb_x, pb_y, gh0_78, gh0_79, gh0_80, \
                         gh1_81, gh1_82, gh1_83, gi_142, gi_143, \
                         gi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_9 * gh0_78[k]
                   - f_10 * gh1_81[k]
                   + pb_x[k] * gi_142[k];

        t_314[k] = f_9 * gh0_79[k]
                   - f_10 * gh1_82[k]
                   + pb_x[k] * gi_143[k];

        t_315[k] = f_7 * gh0_80[k]
                   - f_8 * gh1_83[k]
                   + pb_x[k] * gi_144[k];

        t_316[k] = pb_y[k] * gi_143[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_x, pb_y, gh0_81, gh0_82, gh0_83, \
                         gh1_84, gh1_85, gh1_86, gi_145, gi_146, \
                         gi_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_7 * gh0_81[k]
                   - f_8 * gh1_84[k]
                   + pb_x[k] * gi_145[k];

        t_318[k] = f_5 * gh0_82[k]
                   - f_6 * gh1_85[k]
                   + pb_x[k] * gi_146[k];

        t_319[k] = f_5 * gh0_83[k]
                   - f_6 * gh1_86[k]
                   + pb_x[k] * gi_147[k];

        t_320[k] = pb_y[k] * gi_145[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_x, gh0_84, gh0_85, gh0_86, gh1_87, gh1_88, \
                         gh1_89, gi_148, gi_149, gi_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_5 * gh0_84[k]
                   - f_6 * gh1_87[k]
                   + pb_x[k] * gi_148[k];

        t_322[k] = f_3 * gh0_85[k]
                   - f_4 * gh1_88[k]
                   + pb_x[k] * gi_149[k];

        t_323[k] = f_3 * gh0_86[k]
                   - f_4 * gh1_89[k]
                   + pb_x[k] * gi_150[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pb_x, pb_y, gh0_87, gh0_89, \
                         gh1_90, gh1_92, gi_148, gi_151, gi_152, gi_153, \
                         gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_3 * gh0_87[k]
                   - f_4 * gh1_90[k]
                   + pb_x[k] * gi_151[k];

        t_325[k] = pb_y[k] * gi_148[k];

        t_326[k] = f_3 * gh0_89[k]
                   - f_4 * gh1_92[k]
                   + pb_x[k] * gi_152[k];

        t_327[k] = pb_x[k] * gi_153[k];

        t_328[k] = pb_x[k] * gi_154[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, pb_x, pb_y, pb_z, fi_96, gh0_85, \
                         gh1_88, gi_153, gi_155, gi_156, gi_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pb_x[k] * gi_155[k];

        t_330[k] = pb_x[k] * gi_156[k];

        t_331[k] = pb_x[k] * gi_158[k];

        t_332[k] = f_1 * gh0_85[k]
                   - f_2 * gh1_88[k]
                   + pb_y[k] * gi_153[k];

        t_333[k] = f_0 * fi_96[k]
                   + pb_z[k] * gi_153[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pb_y, gh0_86, gh0_87, gh0_88, gh1_89, gh1_90, \
                         gh1_91, gi_154, gi_155, gi_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_9 * gh0_86[k]
                   - f_10 * gh1_89[k]
                   + pb_y[k] * gi_154[k];

        t_335[k] = f_7 * gh0_87[k]
                   - f_8 * gh1_90[k]
                   + pb_y[k] * gi_155[k];

        t_336[k] = f_5 * gh0_88[k]
                   - f_6 * gh1_91[k]
                   + pb_y[k] * gi_156[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pb_y, pb_z, fi_101, gh0_89, gh1_92, gi_157, \
                         gi_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * gh0_89[k]
                   - f_4 * gh1_92[k]
                   + pb_y[k] * gi_157[k];

        t_338[k] = pb_y[k] * gi_158[k];

        t_339[k] = f_0 * fi_101[k]
                   + f_1 * gh0_89[k]
                   - f_2 * gh1_92[k]
                   + pb_z[k] * gi_158[k];
    }
}

auto
compute_prim_gk_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_11 = 1.0 / p;
    const auto f_12 = 1.5 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_24 = buffer.data(dk1 + 24);
    const auto *dk1_44 = buffer.data(dk1 + 44);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_21 = buffer.data(fi + 21);
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
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
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
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_45 = buffer.data(fk + 45);
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
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_77 = buffer.data(fk + 77);

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
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_74 = buffer.data(gh0 + 74);
    const auto *gh0_75 = buffer.data(gh0 + 75);
    const auto *gh0_80 = buffer.data(gh0 + 80);
    const auto *gh0_81 = buffer.data(gh0 + 81);
    const auto *gh0_82 = buffer.data(gh0 + 82);
    const auto *gh0_83 = buffer.data(gh0 + 83);
    const auto *gh0_84 = buffer.data(gh0 + 84);
    const auto *gh0_85 = buffer.data(gh0 + 85);
    const auto *gh0_86 = buffer.data(gh0 + 86);
    const auto *gh0_87 = buffer.data(gh0 + 87);
    const auto *gh0_88 = buffer.data(gh0 + 88);
    const auto *gh0_89 = buffer.data(gh0 + 89);
    const auto *gh0_90 = buffer.data(gh0 + 90);
    const auto *gh0_91 = buffer.data(gh0 + 91);
    const auto *gh0_92 = buffer.data(gh0 + 92);

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
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_62 = buffer.data(gh1 + 62);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_74 = buffer.data(gh1 + 74);
    const auto *gh1_75 = buffer.data(gh1 + 75);
    const auto *gh1_76 = buffer.data(gh1 + 76);
    const auto *gh1_77 = buffer.data(gh1 + 77);
    const auto *gh1_83 = buffer.data(gh1 + 83);
    const auto *gh1_84 = buffer.data(gh1 + 84);
    const auto *gh1_85 = buffer.data(gh1 + 85);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_87 = buffer.data(gh1 + 87);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_89 = buffer.data(gh1 + 89);
    const auto *gh1_90 = buffer.data(gh1 + 90);
    const auto *gh1_91 = buffer.data(gh1 + 91);
    const auto *gh1_92 = buffer.data(gh1 + 92);
    const auto *gh1_93 = buffer.data(gh1 + 93);
    const auto *gh1_94 = buffer.data(gh1 + 94);
    const auto *gh1_95 = buffer.data(gh1 + 95);

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
    const auto *gi_20 = buffer.data(gi + 20);
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
    const auto *gi_60 = buffer.data(gi + 60);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, gi_0, \
                         gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, gh1_2, \
                         gh1_3, gi_3, gi_4, gi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_6[k] = pb_y[k] * gi_4[k];

        t_7[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_8[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, gh0_4, gh0_5, gh1_4, gh1_5, gi_6, \
                         gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * gh0_4[k]
                 - f_4 * gh1_4[k]
                 + pb_y[k] * gi_6[k];

        t_10[k] = pb_y[k] * gi_7[k];

        t_11[k] = f_7 * gh0_4[k]
                  - f_8 * gh1_4[k]
                  + pb_z[k] * gi_7[k];

        t_12[k] = f_9 * gh0_5[k]
                  - f_10 * gh1_5[k]
                  + pb_y[k] * gi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, gh0_6, gh0_7, gh1_6, gh1_7, gi_9, \
                         gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gh0_6[k]
                  - f_6 * gh1_6[k]
                  + pb_y[k] * gi_9[k];

        t_14[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_15[k] = pb_y[k] * gi_11[k];

        t_16[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, gh0_8, gh0_9, gh0_10, gh1_8, gh1_9, gh1_10, \
                         gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_18[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_9[k]
                  + pb_y[k] * gi_13[k];

        t_19[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_10[k]
                  + pb_y[k] * gi_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, pb_z, fk_0, gh0_11, gh0_12, \
                         gh1_11, gh1_12, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_11[k]
                  + pb_y[k] * gi_15[k];

        t_21[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_12[k]
                  + pb_y[k] * gi_16[k];

        t_22[k] = pb_y[k] * gi_17[k];

        t_23[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_12[k]
                  + pb_z[k] * gi_17[k];

        t_24[k] = pa_y[k] * fk_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, fi_1, fi_3, fi_5, fi_8, fi_12, \
                         fk_3, fk_5, fk_8, fk_12, fk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * fi_1[k]
                  + pa_y[k] * fk_3[k];

        t_26[k] = f_12 * fi_3[k]
                  + pa_y[k] * fk_5[k];

        t_27[k] = f_0 * fi_5[k]
                  + pa_y[k] * fk_8[k];

        t_28[k] = f_13 * fi_8[k]
                  + pa_y[k] * fk_12[k];

        t_29[k] = f_14 * fi_12[k]
                  + pa_y[k] * fk_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_z, pb_z, fi_0, fi_2, fi_4, fi_7, \
                         fk_0, fk_4, fk_7, fk_11, gi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * fk_0[k];

        t_31[k] = f_15 * fi_0[k]
                  + pb_z[k] * gi_20[k];

        t_32[k] = f_11 * fi_2[k]
                  + pa_z[k] * fk_4[k];

        t_33[k] = f_12 * fi_4[k]
                  + pa_z[k] * fk_7[k];

        t_34[k] = f_0 * fi_7[k]
                  + pa_z[k] * fk_11[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, fi_11, fi_13, fi_14, fi_15, \
                         fi_16, fk_16, fk_18, fk_19, fk_20, fk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * fi_11[k]
                  + pa_z[k] * fk_16[k];

        t_36[k] = f_11 * fi_13[k]
                  + pa_z[k] * fk_18[k];

        t_37[k] = f_12 * fi_14[k]
                  + pa_z[k] * fk_19[k];

        t_38[k] = f_0 * fi_15[k]
                  + pa_z[k] * fk_20[k];

        t_39[k] = f_13 * fi_16[k]
                  + pa_z[k] * fk_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_z, dk0_0, dk1_0, fi_18, fk_23, \
                         fk_24, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_14 * fi_18[k]
                  + pa_z[k] * fk_23[k];

        t_41[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_y[k] * fk_24[k];

        t_42[k] = pb_z[k] * gi_22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_z, fi_24, fi_25, gh0_15, gh0_17, gh0_19, \
                         gh1_16, gh1_18, gh1_20, gi_23, gi_24, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * fi_24[k]
                  + f_9 * gh0_17[k]
                  - f_10 * gh1_18[k]
                  + pb_x[k] * gi_24[k];

        t_44[k] = f_3 * gh0_15[k]
                  - f_4 * gh1_16[k]
                  + pb_z[k] * gi_23[k];

        t_45[k] = f_11 * fi_25[k]
                  + f_7 * gh0_19[k]
                  - f_8 * gh1_20[k]
                  + pb_x[k] * gi_26[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pb_z, fi_26, gh0_16, gh0_22, gh1_17, \
                         gh1_23, gi_24, gi_25, gi_26, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * gi_24[k];

        t_47[k] = f_5 * gh0_16[k]
                  - f_6 * gh1_17[k]
                  + pb_z[k] * gi_25[k];

        t_48[k] = f_11 * fi_26[k]
                  + f_5 * gh0_22[k]
                  - f_6 * gh1_23[k]
                  + pb_x[k] * gi_29[k];

        t_49[k] = pb_z[k] * gi_26[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, pb_z, fi_27, gh0_17, gh0_18, gh0_23, gh1_18, \
                         gh1_19, gh1_24, gi_27, gi_28, gi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * gh0_17[k]
                  - f_4 * gh1_18[k]
                  + pb_z[k] * gi_27[k];

        t_51[k] = f_7 * gh0_18[k]
                  - f_8 * gh1_19[k]
                  + pb_z[k] * gi_28[k];

        t_52[k] = f_11 * fi_27[k]
                  + f_3 * gh0_23[k]
                  - f_4 * gh1_24[k]
                  + pb_x[k] * gi_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_z, gh0_19, gh0_20, gh0_21, gh1_20, gh1_21, \
                         gh1_22, gi_29, gi_30, gi_31, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_z[k] * gi_29[k];

        t_54[k] = f_3 * gh0_19[k]
                  - f_4 * gh1_20[k]
                  + pb_z[k] * gi_30[k];

        t_55[k] = f_5 * gh0_20[k]
                  - f_6 * gh1_21[k]
                  + pb_z[k] * gi_31[k];

        t_56[k] = f_9 * gh0_21[k]
                  - f_10 * gh1_22[k]
                  + pb_z[k] * gi_32[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_x, pb_z, dk0_1, dk1_24, fi_28, \
                         fk_26, gh0_23, gh1_24, gi_34, gi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_11 * fi_28[k]
                  + pb_x[k] * gi_34[k];

        t_58[k] = f_16 * dk0_1[k]
                  - f_17 * dk1_24[k]
                  + pa_x[k] * fk_26[k];

        t_59[k] = pb_z[k] * gi_34[k];

        t_60[k] = f_3 * gh0_23[k]
                  - f_4 * gh1_24[k]
                  + pb_z[k] * gi_35[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_z, gh0_24, gh0_25, gh0_26, gh1_25, gh1_26, \
                         gh1_27, gi_36, gi_37, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * gh0_24[k]
                  - f_6 * gh1_25[k]
                  + pb_z[k] * gi_36[k];

        t_62[k] = f_7 * gh0_25[k]
                  - f_8 * gh1_26[k]
                  + pb_z[k] * gi_37[k];

        t_63[k] = f_9 * gh0_26[k]
                  - f_10 * gh1_27[k]
                  + pb_z[k] * gi_38[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_y, pb_z, dk0_0, dk1_0, fi_21, fk_25, \
                         gh0_27, gh1_28, gi_39, gi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * gh0_27[k]
                  - f_2 * gh1_28[k]
                  + pb_z[k] * gi_39[k];

        t_65[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_z[k] * fk_25[k];

        t_66[k] = pb_y[k] * gi_40[k];

        t_67[k] = f_11 * fi_21[k]
                  + pb_z[k] * gi_40[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, fi_30, gh0_28, gh0_29, gh0_31, \
                         gh1_29, gh1_30, gh1_32, gi_41, gi_42, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * gh0_28[k]
                  - f_4 * gh1_29[k]
                  + pb_y[k] * gi_41[k];

        t_69[k] = f_11 * fi_30[k]
                  + f_9 * gh0_31[k]
                  - f_10 * gh1_32[k]
                  + pb_x[k] * gi_43[k];

        t_70[k] = f_5 * gh0_29[k]
                  - f_6 * gh1_30[k]
                  + pb_y[k] * gi_42[k];

        t_71[k] = pb_y[k] * gi_43[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, pb_y, fi_31, gh0_30, gh0_31, gh0_34, \
                         gh1_31, gh1_32, gh1_35, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * fi_31[k]
                  + f_7 * gh0_34[k]
                  - f_8 * gh1_35[k]
                  + pb_x[k] * gi_46[k];

        t_73[k] = f_7 * gh0_30[k]
                  - f_8 * gh1_31[k]
                  + pb_y[k] * gi_44[k];

        t_74[k] = f_3 * gh0_31[k]
                  - f_4 * gh1_32[k]
                  + pb_y[k] * gi_45[k];

        t_75[k] = pb_y[k] * gi_46[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_x, pb_y, fi_32, gh0_32, gh0_33, gh0_35, gh1_33, \
                         gh1_34, gh1_36, gi_47, gi_48, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_11 * fi_32[k]
                  + f_5 * gh0_35[k]
                  - f_6 * gh1_36[k]
                  + pb_x[k] * gi_50[k];

        t_77[k] = f_9 * gh0_32[k]
                  - f_10 * gh1_33[k]
                  + pb_y[k] * gi_47[k];

        t_78[k] = f_5 * gh0_33[k]
                  - f_6 * gh1_34[k]
                  + pb_y[k] * gi_48[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pb_y, fi_33, fi_34, gh0_34, gh0_40, \
                         gh1_35, gh1_41, gi_49, gi_50, gi_51, gi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * gh0_34[k]
                  - f_4 * gh1_35[k]
                  + pb_y[k] * gi_49[k];

        t_80[k] = pb_y[k] * gi_50[k];

        t_81[k] = f_11 * fi_33[k]
                  + f_3 * gh0_40[k]
                  - f_4 * gh1_41[k]
                  + pb_x[k] * gi_51[k];

        t_82[k] = f_11 * fi_34[k]
                  + pb_x[k] * gi_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, gh0_36, gh0_37, gh0_38, gh1_37, gh1_38, \
                         gh1_39, gi_52, gi_53, gi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * gh0_36[k]
                  - f_2 * gh1_37[k]
                  + pb_y[k] * gi_52[k];

        t_84[k] = f_9 * gh0_37[k]
                  - f_10 * gh1_38[k]
                  + pb_y[k] * gi_53[k];

        t_85[k] = f_7 * gh0_38[k]
                  - f_8 * gh1_39[k]
                  + pb_y[k] * gi_54[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pb_y, dk0_2, dk1_44, fk_27, gh0_39, \
                         gh0_40, gh1_40, gh1_41, gi_55, gi_56, gi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_5 * gh0_39[k]
                  - f_6 * gh1_40[k]
                  + pb_y[k] * gi_55[k];

        t_87[k] = f_3 * gh0_40[k]
                  - f_4 * gh1_41[k]
                  + pb_y[k] * gi_56[k];

        t_88[k] = pb_y[k] * gi_57[k];

        t_89[k] = f_16 * dk0_2[k]
                  - f_17 * dk1_44[k]
                  + pa_x[k] * fk_27[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_x, fi_35, fi_37, fi_39, fi_41, \
                         fi_44, fk_28, fk_29, fk_31, fk_33, fk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_14 * fi_35[k]
                  + pa_x[k] * fk_28[k];

        t_91[k] = f_13 * fi_37[k]
                  + pa_x[k] * fk_29[k];

        t_92[k] = f_0 * fi_39[k]
                  + pa_x[k] * fk_31[k];

        t_93[k] = f_12 * fi_41[k]
                  + pa_x[k] * fk_33[k];

        t_94[k] = f_11 * fi_44[k]
                  + pa_x[k] * fk_36[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pa_x, pb_z, fi_29, fi_62, fi_65, fi_67, \
                         fk_45, fk_54, fk_56, fk_58, gi_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pa_x[k] * fk_45[k];

        t_96[k] = f_14 * fi_62[k]
                  + pa_x[k] * fk_54[k];

        t_97[k] = f_12 * fi_29[k]
                  + pb_z[k] * gi_60[k];

        t_98[k] = f_13 * fi_65[k]
                  + pa_x[k] * fk_56[k];

        t_99[k] = f_0 * fi_67[k]
                  + pa_x[k] * fk_58[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, fi_70, fi_74, fk_61, fk_65, \
                         fk_77, gh0_49, gh1_50, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_12 * fi_70[k]
                   + pa_x[k] * fk_61[k];

        t_101[k] = f_11 * fi_74[k]
                   + pa_x[k] * fk_65[k];

        t_102[k] = pa_x[k] * fk_77[k];

        t_103[k] = f_1 * gh0_49[k]
                   - f_2 * gh1_50[k]
                   + pb_x[k] * gi_62[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_x, gh0_50, gh0_51, gh0_52, gh1_51, gh1_52, \
                         gh1_53, gi_63, gi_64, gi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_9 * gh0_50[k]
                   - f_10 * gh1_51[k]
                   + pb_x[k] * gi_63[k];

        t_105[k] = f_9 * gh0_51[k]
                   - f_10 * gh1_52[k]
                   + pb_x[k] * gi_64[k];

        t_106[k] = f_7 * gh0_52[k]
                   - f_8 * gh1_53[k]
                   + pb_x[k] * gi_65[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, gh0_53, gh0_54, gh0_55, gh1_54, gh1_55, \
                         gh1_56, gi_66, gi_67, gi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_7 * gh0_53[k]
                   - f_8 * gh1_54[k]
                   + pb_x[k] * gi_66[k];

        t_108[k] = f_5 * gh0_54[k]
                   - f_6 * gh1_55[k]
                   + pb_x[k] * gi_67[k];

        t_109[k] = f_5 * gh0_55[k]
                   - f_6 * gh1_56[k]
                   + pb_x[k] * gi_68[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, gh0_56, gh0_57, gh0_59, gh1_57, gh1_58, \
                         gh1_60, gi_69, gi_70, gi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_5 * gh0_56[k]
                   - f_6 * gh1_57[k]
                   + pb_x[k] * gi_69[k];

        t_111[k] = f_3 * gh0_57[k]
                   - f_4 * gh1_58[k]
                   + pb_x[k] * gi_70[k];

        t_112[k] = f_3 * gh0_59[k]
                   - f_4 * gh1_60[k]
                   + pb_x[k] * gi_71[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pb_x, gh0_60, gh0_61, gh1_61, \
                         gh1_62, gi_72, gi_73, gi_74, gi_76, gi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * gh0_60[k]
                   - f_4 * gh1_61[k]
                   + pb_x[k] * gi_72[k];

        t_114[k] = f_3 * gh0_61[k]
                   - f_4 * gh1_62[k]
                   + pb_x[k] * gi_73[k];

        t_115[k] = pb_x[k] * gi_74[k];

        t_116[k] = pb_x[k] * gi_76[k];

        t_117[k] = pb_x[k] * gi_77[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pb_x, pb_y, pb_z, fi_48, gh0_57, \
                         gh1_58, gi_74, gi_75, gi_78, gi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * gi_78[k];

        t_119[k] = pb_x[k] * gi_79[k];

        t_120[k] = f_0 * fi_48[k]
                   + f_1 * gh0_57[k]
                   - f_2 * gh1_58[k]
                   + pb_y[k] * gi_74[k];

        t_121[k] = pb_z[k] * gi_74[k];

        t_122[k] = f_3 * gh0_57[k]
                   - f_4 * gh1_58[k]
                   + pb_z[k] * gi_75[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_z, gh0_58, gh0_59, gh0_60, gh1_59, gh1_60, \
                         gh1_61, gi_76, gi_77, gi_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_5 * gh0_58[k]
                   - f_6 * gh1_59[k]
                   + pb_z[k] * gi_76[k];

        t_124[k] = f_7 * gh0_59[k]
                   - f_8 * gh1_60[k]
                   + pb_z[k] * gi_77[k];

        t_125[k] = f_9 * gh0_60[k]
                   - f_10 * gh1_61[k]
                   + pb_z[k] * gi_78[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_z, pb_y, pb_z, fi_36, fi_38, fi_53, \
                         fk_30, fk_32, gh0_61, gh1_62, gi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_0 * fi_53[k]
                   + pb_y[k] * gi_79[k];

        t_127[k] = f_1 * gh0_61[k]
                   - f_2 * gh1_62[k]
                   + pb_z[k] * gi_79[k];

        t_128[k] = f_11 * fi_36[k]
                   + pa_z[k] * fk_30[k];

        t_129[k] = f_12 * fi_38[k]
                   + pa_z[k] * fk_32[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_z, pb_z, fi_40, fi_43, fi_48, \
                         fi_49, fk_35, fk_39, fk_45, fk_47, gi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * fi_40[k]
                   + pa_z[k] * fk_35[k];

        t_131[k] = f_13 * fi_43[k]
                   + pa_z[k] * fk_39[k];

        t_132[k] = pa_z[k] * fk_45[k];

        t_133[k] = f_15 * fi_48[k]
                   + pb_z[k] * gi_80[k];

        t_134[k] = f_11 * fi_49[k]
                   + pa_z[k] * fk_47[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, fi_50, fi_51, fi_52, fi_55, \
                         fk_48, fk_49, fk_50, gi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_12 * fi_50[k]
                   + pa_z[k] * fk_48[k];

        t_136[k] = f_0 * fi_51[k]
                   + pa_z[k] * fk_49[k];

        t_137[k] = f_13 * fi_52[k]
                   + pa_z[k] * fk_50[k];

        t_138[k] = f_12 * fi_55[k]
                   + pb_y[k] * gi_81[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_z, pb_x, fi_53, fk_51, gh0_63, gh0_64, \
                         gh1_65, gh1_66, gi_82, gi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * fi_53[k]
                   + pa_z[k] * fk_51[k];

        t_140[k] = f_1 * gh0_63[k]
                   - f_2 * gh1_65[k]
                   + pb_x[k] * gi_82[k];

        t_141[k] = f_9 * gh0_64[k]
                   - f_10 * gh1_66[k]
                   + pb_x[k] * gi_83[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_x, gh0_65, gh0_66, gh0_67, gh1_67, gh1_68, \
                         gh1_69, gi_84, gi_85, gi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_9 * gh0_65[k]
                   - f_10 * gh1_67[k]
                   + pb_x[k] * gi_84[k];

        t_143[k] = f_7 * gh0_66[k]
                   - f_8 * gh1_68[k]
                   + pb_x[k] * gi_85[k];

        t_144[k] = f_7 * gh0_67[k]
                   - f_8 * gh1_69[k]
                   + pb_x[k] * gi_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, gh0_68, gh0_69, gh0_70, gh1_70, gh1_71, \
                         gh1_72, gi_87, gi_88, gi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_5 * gh0_68[k]
                   - f_6 * gh1_70[k]
                   + pb_x[k] * gi_87[k];

        t_146[k] = f_5 * gh0_69[k]
                   - f_6 * gh1_71[k]
                   + pb_x[k] * gi_88[k];

        t_147[k] = f_5 * gh0_70[k]
                   - f_6 * gh1_72[k]
                   + pb_x[k] * gi_89[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pb_x, gh0_71, gh0_72, gh0_73, gh1_73, gh1_74, \
                         gh1_75, gi_90, gi_91, gi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_3 * gh0_71[k]
                   - f_4 * gh1_73[k]
                   + pb_x[k] * gi_90[k];

        t_149[k] = f_3 * gh0_72[k]
                   - f_4 * gh1_74[k]
                   + pb_x[k] * gi_91[k];

        t_150[k] = f_3 * gh0_73[k]
                   - f_4 * gh1_75[k]
                   + pb_x[k] * gi_92[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, t_156, pb_x, gh0_75, gh1_77, \
                         gi_93, gi_94, gi_95, gi_96, gi_97, gi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_3 * gh0_75[k]
                   - f_4 * gh1_77[k]
                   + pb_x[k] * gi_93[k];

        t_152[k] = pb_x[k] * gi_94[k];

        t_153[k] = pb_x[k] * gi_95[k];

        t_154[k] = pb_x[k] * gi_96[k];

        t_155[k] = pb_x[k] * gi_97[k];

        t_156[k] = pb_x[k] * gi_99[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pa_z, pb_y, pb_z, dk0_1, dk1_24, fi_54, fi_57, \
                         fk_52, gh0_72, gh1_74, gi_94, gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_16 * dk0_1[k]
                   - f_17 * dk1_24[k]
                   + pa_z[k] * fk_52[k];

        t_158[k] = f_11 * fi_54[k]
                   + pb_z[k] * gi_94[k];

        t_159[k] = f_11 * fi_57[k]
                   + f_9 * gh0_72[k]
                   - f_10 * gh1_74[k]
                   + pb_y[k] * gi_95[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pb_y, fi_58, fi_59, fi_60, gh0_73, gh0_74, \
                         gh0_75, gh1_75, gh1_76, gh1_77, gi_96, gi_97, \
                         gi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_11 * fi_58[k]
                   + f_7 * gh0_73[k]
                   - f_8 * gh1_75[k]
                   + pb_y[k] * gi_96[k];

        t_161[k] = f_11 * fi_59[k]
                   + f_5 * gh0_74[k]
                   - f_6 * gh1_76[k]
                   + pb_y[k] * gi_97[k];

        t_162[k] = f_11 * fi_60[k]
                   + f_3 * gh0_75[k]
                   - f_4 * gh1_77[k]
                   + pb_y[k] * gi_98[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pb_y, dk0_2, dk1_44, fi_61, fi_63, \
                         fi_64, fk_53, fk_55, fk_57, gi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * fi_61[k]
                   + pb_y[k] * gi_99[k];

        t_164[k] = f_16 * dk0_2[k]
                   - f_17 * dk1_44[k]
                   + pa_y[k] * fk_53[k];

        t_165[k] = f_11 * fi_63[k]
                   + pa_y[k] * fk_55[k];

        t_166[k] = f_12 * fi_64[k]
                   + pa_y[k] * fk_57[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pb_z, fi_56, fi_66, fi_68, fi_75, \
                         fk_59, fk_62, fk_71, gi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * fi_66[k]
                   + pa_y[k] * fk_59[k];

        t_168[k] = f_13 * fi_68[k]
                   + pa_y[k] * fk_62[k];

        t_169[k] = f_14 * fi_75[k]
                   + pa_y[k] * fk_71[k];

        t_170[k] = f_12 * fi_56[k]
                   + pb_z[k] * gi_100[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, fi_76, fi_77, fi_78, fi_79, fk_72, \
                         fk_73, fk_74, fk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_13 * fi_76[k]
                   + pa_y[k] * fk_72[k];

        t_172[k] = f_0 * fi_77[k]
                   + pa_y[k] * fk_73[k];

        t_173[k] = f_12 * fi_78[k]
                   + pa_y[k] * fk_74[k];

        t_174[k] = f_11 * fi_79[k]
                   + pa_y[k] * fk_75[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pb_x, pb_y, pb_z, fi_62, fi_80, \
                         fk_77, gh0_80, gh1_83, gi_101, gi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_15 * fi_80[k]
                   + pb_y[k] * gi_101[k];

        t_176[k] = pa_y[k] * fk_77[k];

        t_177[k] = f_1 * gh0_80[k]
                   - f_2 * gh1_83[k]
                   + pb_x[k] * gi_102[k];

        t_178[k] = f_0 * fi_62[k]
                   + pb_z[k] * gi_102[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pb_x, gh0_81, gh0_82, gh0_83, gh1_84, gh1_85, \
                         gh1_86, gi_103, gi_104, gi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_9 * gh0_81[k]
                   - f_10 * gh1_84[k]
                   + pb_x[k] * gi_103[k];

        t_180[k] = f_9 * gh0_82[k]
                   - f_10 * gh1_85[k]
                   + pb_x[k] * gi_104[k];

        t_181[k] = f_7 * gh0_83[k]
                   - f_8 * gh1_86[k]
                   + pb_x[k] * gi_105[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pb_x, gh0_84, gh0_85, gh0_86, gh1_87, gh1_88, \
                         gh1_89, gi_106, gi_107, gi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_7 * gh0_84[k]
                   - f_8 * gh1_87[k]
                   + pb_x[k] * gi_106[k];

        t_183[k] = f_5 * gh0_85[k]
                   - f_6 * gh1_88[k]
                   + pb_x[k] * gi_107[k];

        t_184[k] = f_5 * gh0_86[k]
                   - f_6 * gh1_89[k]
                   + pb_x[k] * gi_108[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pb_x, gh0_87, gh0_88, gh0_89, gh1_90, gh1_91, \
                         gh1_92, gi_109, gi_110, gi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_5 * gh0_87[k]
                   - f_6 * gh1_90[k]
                   + pb_x[k] * gi_109[k];

        t_186[k] = f_3 * gh0_88[k]
                   - f_4 * gh1_91[k]
                   + pb_x[k] * gi_110[k];

        t_187[k] = f_3 * gh0_89[k]
                   - f_4 * gh1_92[k]
                   + pb_x[k] * gi_111[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pb_x, gh0_90, gh0_92, gh1_93, \
                         gh1_95, gi_112, gi_113, gi_114, gi_115, \
                         gi_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_3 * gh0_90[k]
                   - f_4 * gh1_93[k]
                   + pb_x[k] * gi_112[k];

        t_189[k] = f_3 * gh0_92[k]
                   - f_4 * gh1_95[k]
                   + pb_x[k] * gi_113[k];

        t_190[k] = pb_x[k] * gi_114[k];

        t_191[k] = pb_x[k] * gi_115[k];

        t_192[k] = pb_x[k] * gi_116[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_x, pb_y, pb_z, fi_75, gh0_88, gh1_91, \
                         gi_114, gi_117, gi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * gi_117[k];

        t_194[k] = pb_x[k] * gi_119[k];

        t_195[k] = f_1 * gh0_88[k]
                   - f_2 * gh1_91[k]
                   + pb_y[k] * gi_114[k];

        t_196[k] = f_0 * fi_75[k]
                   + pb_z[k] * gi_114[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, gh0_89, gh0_90, gh0_91, gh1_92, gh1_93, \
                         gh1_94, gi_115, gi_116, gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_9 * gh0_89[k]
                   - f_10 * gh1_92[k]
                   + pb_y[k] * gi_115[k];

        t_198[k] = f_7 * gh0_90[k]
                   - f_8 * gh1_93[k]
                   + pb_y[k] * gi_116[k];

        t_199[k] = f_5 * gh0_91[k]
                   - f_6 * gh1_94[k]
                   + pb_y[k] * gi_117[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_y, pb_z, fi_80, gh0_92, gh1_95, gi_118, \
                         gi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_3 * gh0_92[k]
                   - f_4 * gh1_95[k]
                   + pb_y[k] * gi_118[k];

        t_201[k] = pb_y[k] * gi_119[k];

        t_202[k] = f_0 * fi_80[k]
                   + f_1 * gh0_92[k]
                   - f_2 * gh1_95[k]
                   + pb_z[k] * gi_119[k];
    }
}

auto
compute_prim_gk_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_20 = buffer.data(fi + 20);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_23 = buffer.data(gh1 + 23);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_26 = buffer.data(gi + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_7, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_9, gh1_10, gh1_11, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_12, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_13, fk_5, fk_6, \
                         gh0_15, gh1_15, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_13[k]
                  + f_1 * gh0_15[k]
                  - f_2 * gh1_15[k]
                  + pb_y[k] * gi_17[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_15, fi_16, fi_17, gh0_18, gh0_19, gh0_20, \
                         gh1_18, gh1_19, gh1_20, gi_20, gi_21, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_15[k]
                  + f_6 * gh0_18[k]
                  - f_7 * gh1_18[k]
                  + pb_y[k] * gi_20[k];

        t_21[k] = f_5 * fi_16[k]
                  + f_8 * gh0_19[k]
                  - f_9 * gh1_19[k]
                  + pb_y[k] * gi_21[k];

        t_22[k] = f_5 * fi_17[k]
                  + f_10 * gh0_20[k]
                  - f_11 * gh1_20[k]
                  + pb_y[k] * gi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_18, fk_7, fk_8, \
                         gh0_21, gh1_21, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_18[k]
                  + f_12 * gh0_21[k]
                  - f_13 * gh1_21[k]
                  + pb_y[k] * gi_23[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_20, gh0_23, gh1_23, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_20[k]
                  + f_1 * gh0_23[k]
                  - f_2 * gh1_23[k]
                  + pb_z[k] * gi_26[k];
    }
}

auto
compute_prim_gk_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_20 = buffer.data(fi + 20);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_43 = buffer.data(gh1 + 43);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_91 = buffer.data(gh1 + 91);
    const auto *gh1_92 = buffer.data(gh1 + 92);
    const auto *gh1_93 = buffer.data(gh1 + 93);
    const auto *gh1_94 = buffer.data(gh1 + 94);
    const auto *gh1_117 = buffer.data(gh1 + 117);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_145 = buffer.data(gi + 145);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_26, \
                         gh1_28, gh1_30, gi_33, gi_35, gi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_26[k]
                 + pb_x[k] * gi_33[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_28[k]
                 + pb_x[k] * gi_35[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_30[k]
                 + pb_x[k] * gi_37[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_31, gi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_31[k]
                 + pb_x[k] * gi_39[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_41, gh1_43, gh1_44, gi_48, gi_50, gi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_41[k]
                  + pb_x[k] * gi_48[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_43[k]
                  + pb_x[k] * gi_50[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_44[k]
                  + pb_x[k] * gi_52[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_49, gi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_49[k]
                  + pb_x[k] * gi_53[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_13, fk_5, fk_6, \
                         gh0_15, gh1_71, gi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_13[k]
                  + f_1 * gh0_15[k]
                  - f_2 * gh1_71[k]
                  + pb_y[k] * gi_84[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_15, fi_16, fi_17, gh0_18, gh0_19, gh0_20, \
                         gh1_91, gh1_92, gh1_93, gi_111, gi_112, \
                         gi_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_15[k]
                  + f_6 * gh0_18[k]
                  - f_7 * gh1_91[k]
                  + pb_y[k] * gi_111[k];

        t_21[k] = f_5 * fi_16[k]
                  + f_8 * gh0_19[k]
                  - f_9 * gh1_92[k]
                  + pb_y[k] * gi_112[k];

        t_22[k] = f_5 * fi_17[k]
                  + f_10 * gh0_20[k]
                  - f_11 * gh1_93[k]
                  + pb_y[k] * gi_113[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_18, fk_7, fk_8, \
                         gh0_21, gh1_94, gi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_18[k]
                  + f_12 * gh0_21[k]
                  - f_13 * gh1_94[k]
                  + pb_y[k] * gi_114[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_20, gh0_23, gh1_117, gi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_20[k]
                  + f_1 * gh0_23[k]
                  - f_2 * gh1_117[k]
                  + pb_z[k] * gi_145[k];
    }
}

auto
compute_prim_gk_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);
    const auto f_13 = 1.0 / p;

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

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_53 = buffer.data(fi + 53);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_1 = buffer.data(gh0 + 1);
    const auto *gh0_2 = buffer.data(gh0 + 2);
    const auto *gh0_3 = buffer.data(gh0 + 3);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_14 = buffer.data(gh0 + 14);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_41 = buffer.data(gh0 + 41);
    const auto *gh0_43 = buffer.data(gh0 + 43);
    const auto *gh0_44 = buffer.data(gh0 + 44);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_74 = buffer.data(gh0 + 74);
    const auto *gh0_75 = buffer.data(gh0 + 75);
    const auto *gh0_91 = buffer.data(gh0 + 91);
    const auto *gh0_92 = buffer.data(gh0 + 92);
    const auto *gh0_93 = buffer.data(gh0 + 93);
    const auto *gh0_94 = buffer.data(gh0 + 94);
    const auto *gh0_103 = buffer.data(gh0 + 103);
    const auto *gh0_105 = buffer.data(gh0 + 105);
    const auto *gh0_106 = buffer.data(gh0 + 106);
    const auto *gh0_107 = buffer.data(gh0 + 107);
    const auto *gh0_109 = buffer.data(gh0 + 109);
    const auto *gh0_110 = buffer.data(gh0 + 110);
    const auto *gh0_111 = buffer.data(gh0 + 111);
    const auto *gh0_112 = buffer.data(gh0 + 112);
    const auto *gh0_113 = buffer.data(gh0 + 113);
    const auto *gh0_114 = buffer.data(gh0 + 114);
    const auto *gh0_115 = buffer.data(gh0 + 115);
    const auto *gh0_116 = buffer.data(gh0 + 116);
    const auto *gh0_117 = buffer.data(gh0 + 117);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_1 = buffer.data(gh1 + 1);
    const auto *gh1_2 = buffer.data(gh1 + 2);
    const auto *gh1_3 = buffer.data(gh1 + 3);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_62 = buffer.data(gh1 + 62);
    const auto *gh1_63 = buffer.data(gh1 + 63);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_77 = buffer.data(gh1 + 77);
    const auto *gh1_78 = buffer.data(gh1 + 78);
    const auto *gh1_79 = buffer.data(gh1 + 79);
    const auto *gh1_80 = buffer.data(gh1 + 80);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_89 = buffer.data(gh1 + 89);
    const auto *gh1_90 = buffer.data(gh1 + 90);
    const auto *gh1_91 = buffer.data(gh1 + 91);
    const auto *gh1_92 = buffer.data(gh1 + 92);
    const auto *gh1_93 = buffer.data(gh1 + 93);
    const auto *gh1_94 = buffer.data(gh1 + 94);
    const auto *gh1_95 = buffer.data(gh1 + 95);
    const auto *gh1_96 = buffer.data(gh1 + 96);
    const auto *gh1_97 = buffer.data(gh1 + 97);
    const auto *gh1_98 = buffer.data(gh1 + 98);
    const auto *gh1_99 = buffer.data(gh1 + 99);

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
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
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
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fi_0, gh0_0, gh0_1, gh1_0, \
                         gh1_1, gi_0, gi_1, gi_2, gi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_2[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];

        t_3[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, gh0_2, gh0_3, gh0_4, gh1_2, gh1_3, \
                         gh1_4, gi_4, gi_5, gi_6, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_5[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];

        t_6[k] = f_3 * gh0_4[k]
                 - f_4 * gh1_4[k]
                 + pb_y[k] * gi_6[k];

        t_7[k] = f_7 * gh0_4[k]
                 - f_8 * gh1_4[k]
                 + pb_z[k] * gi_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, gh0_5, gh0_7, gh0_8, gh1_5, gh1_6, \
                         gh1_7, gi_8, gi_9, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * gh0_5[k]
                 - f_10 * gh1_5[k]
                 + pb_y[k] * gi_8[k];

        t_9[k] = f_5 * gh0_7[k]
                 - f_6 * gh1_6[k]
                 + pb_y[k] * gi_9[k];

        t_10[k] = f_3 * gh0_8[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_11[k] = f_9 * gh0_8[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, gh0_9, gh0_11, gh0_12, gh1_8, gh1_10, gh1_11, \
                         gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * gh0_9[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_13[k] = f_9 * gh0_11[k]
                  - f_10 * gh1_10[k]
                  + pb_y[k] * gi_13[k];

        t_14[k] = f_7 * gh0_12[k]
                  - f_8 * gh1_11[k]
                  + pb_y[k] * gi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pb_y, pb_z, fk_0, gh0_13, gh0_14, \
                         gh1_12, gh1_13, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gh0_13[k]
                  - f_6 * gh1_12[k]
                  + pb_y[k] * gi_15[k];

        t_16[k] = f_3 * gh0_14[k]
                  - f_4 * gh1_13[k]
                  + pb_y[k] * gi_16[k];

        t_17[k] = f_1 * gh0_14[k]
                  - f_2 * gh1_13[k]
                  + pb_z[k] * gi_17[k];

        t_18[k] = pa_y[k] * fk_0[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_14, fk_0, fk_1, \
                         gh0_26, gh1_20, gi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * fk_0[k];

        t_20[k] = f_11 * dk0_0[k]
                  - f_12 * dk1_0[k]
                  + pa_y[k] * fk_1[k];

        t_21[k] = f_13 * fi_14[k]
                  + f_9 * gh0_26[k]
                  - f_10 * gh1_20[k]
                  + pb_x[k] * gi_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, fi_15, fi_16, fi_17, gh0_28, gh0_30, gh0_31, \
                         gh1_22, gh1_24, gh1_25, gi_22, gi_23, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_13 * fi_15[k]
                  + f_7 * gh0_28[k]
                  - f_8 * gh1_22[k]
                  + pb_x[k] * gi_22[k];

        t_23[k] = f_13 * fi_16[k]
                  + f_5 * gh0_30[k]
                  - f_6 * gh1_24[k]
                  + pb_x[k] * gi_23[k];

        t_24[k] = f_13 * fi_17[k]
                  + f_3 * gh0_31[k]
                  - f_4 * gh1_25[k]
                  + pb_x[k] * gi_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, \
                         fi_19, fk_2, fk_3, gh0_41, gh1_33, gi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * dk0_1[k]
                  - f_12 * dk1_1[k]
                  + pa_x[k] * fk_3[k];

        t_26[k] = f_11 * dk0_0[k]
                  - f_12 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_27[k] = f_13 * fi_19[k]
                  + f_9 * gh0_41[k]
                  - f_10 * gh1_33[k]
                  + pb_x[k] * gi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, fi_20, fi_21, fi_22, gh0_43, gh0_44, gh0_49, \
                         gh1_35, gh1_36, gh1_41, gi_28, gi_29, gi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_13 * fi_20[k]
                  + f_7 * gh0_43[k]
                  - f_8 * gh1_35[k]
                  + pb_x[k] * gi_28[k];

        t_29[k] = f_13 * fi_21[k]
                  + f_5 * gh0_44[k]
                  - f_6 * gh1_36[k]
                  + pb_x[k] * gi_29[k];

        t_30[k] = f_13 * fi_22[k]
                  + f_3 * gh0_49[k]
                  - f_4 * gh1_41[k]
                  + pb_x[k] * gi_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_x, dk0_2, dk1_2, fk_4, fk_5, fk_8, \
                         gh0_61, gh1_52, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_11 * dk0_2[k]
                  - f_12 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_32[k] = pa_x[k] * fk_5[k];

        t_33[k] = pa_x[k] * fk_8[k];

        t_34[k] = f_1 * gh0_61[k]
                  - f_2 * gh1_52[k]
                  + pb_x[k] * gi_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, gh0_63, gh0_64, gh0_65, gh1_54, gh1_55, \
                         gh1_56, gi_35, gi_36, gi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_9 * gh0_63[k]
                  - f_10 * gh1_54[k]
                  + pb_x[k] * gi_35[k];

        t_36[k] = f_9 * gh0_64[k]
                  - f_10 * gh1_55[k]
                  + pb_x[k] * gi_36[k];

        t_37[k] = f_7 * gh0_65[k]
                  - f_8 * gh1_56[k]
                  + pb_x[k] * gi_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, gh0_67, gh0_68, gh0_69, gh1_57, gh1_58, \
                         gh1_59, gi_38, gi_39, gi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_7 * gh0_67[k]
                  - f_8 * gh1_57[k]
                  + pb_x[k] * gi_38[k];

        t_39[k] = f_5 * gh0_68[k]
                  - f_6 * gh1_58[k]
                  + pb_x[k] * gi_39[k];

        t_40[k] = f_5 * gh0_69[k]
                  - f_6 * gh1_59[k]
                  + pb_x[k] * gi_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_x, gh0_70, gh0_71, gh0_73, gh1_60, gh1_61, \
                         gh1_63, gi_41, gi_42, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * gh0_70[k]
                  - f_6 * gh1_60[k]
                  + pb_x[k] * gi_41[k];

        t_42[k] = f_3 * gh0_71[k]
                  - f_4 * gh1_61[k]
                  + pb_x[k] * gi_42[k];

        t_43[k] = f_3 * gh0_73[k]
                  - f_4 * gh1_63[k]
                  + pb_x[k] * gi_43[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_y, fi_31, gh0_71, gh0_74, gh0_75, gh1_61, \
                         gh1_64, gh1_65, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * gh0_74[k]
                  - f_4 * gh1_64[k]
                  + pb_x[k] * gi_44[k];

        t_45[k] = f_3 * gh0_75[k]
                  - f_4 * gh1_65[k]
                  + pb_x[k] * gi_45[k];

        t_46[k] = f_0 * fi_31[k]
                  + f_1 * gh0_71[k]
                  - f_2 * gh1_61[k]
                  + pb_y[k] * gi_46[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_z, gh0_71, gh0_72, gh0_73, gh1_61, gh1_62, \
                         gh1_63, gi_47, gi_48, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * gh0_71[k]
                  - f_4 * gh1_61[k]
                  + pb_z[k] * gi_47[k];

        t_48[k] = f_5 * gh0_72[k]
                  - f_6 * gh1_62[k]
                  + pb_z[k] * gi_48[k];

        t_49[k] = f_7 * gh0_73[k]
                  - f_8 * gh1_63[k]
                  + pb_z[k] * gi_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_z, dk0_1, dk1_1, fk_5, fk_6, gh0_74, \
                         gh0_75, gh1_64, gh1_65, gi_50, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * gh0_74[k]
                  - f_10 * gh1_64[k]
                  + pb_z[k] * gi_50[k];

        t_51[k] = f_1 * gh0_75[k]
                  - f_2 * gh1_65[k]
                  + pb_z[k] * gi_51[k];

        t_52[k] = pa_z[k] * fk_5[k];

        t_53[k] = f_11 * dk0_1[k]
                  - f_12 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, fi_37, fi_38, fi_39, gh0_91, gh0_92, gh0_93, \
                         gh1_77, gh1_78, gh1_79, gi_54, gi_55, gi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_13 * fi_37[k]
                  + f_9 * gh0_91[k]
                  - f_10 * gh1_77[k]
                  + pb_y[k] * gi_54[k];

        t_55[k] = f_13 * fi_38[k]
                  + f_7 * gh0_92[k]
                  - f_8 * gh1_78[k]
                  + pb_y[k] * gi_55[k];

        t_56[k] = f_13 * fi_39[k]
                  + f_5 * gh0_93[k]
                  - f_6 * gh1_79[k]
                  + pb_y[k] * gi_56[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pb_y, dk0_2, dk1_2, fi_40, fk_7, fk_8, \
                         gh0_94, gh1_80, gi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_13 * fi_40[k]
                  + f_3 * gh0_94[k]
                  - f_4 * gh1_80[k]
                  + pb_y[k] * gi_57[k];

        t_58[k] = f_11 * dk0_2[k]
                  - f_12 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_59[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_x, gh0_103, gh0_105, gh0_106, gh1_86, gh1_88, \
                         gh1_89, gi_60, gi_61, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * gh0_103[k]
                  - f_2 * gh1_86[k]
                  + pb_x[k] * gi_60[k];

        t_61[k] = f_9 * gh0_105[k]
                  - f_10 * gh1_88[k]
                  + pb_x[k] * gi_61[k];

        t_62[k] = f_9 * gh0_106[k]
                  - f_10 * gh1_89[k]
                  + pb_x[k] * gi_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, gh0_107, gh0_109, gh0_110, gh1_90, gh1_91, \
                         gh1_92, gi_63, gi_64, gi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * gh0_107[k]
                  - f_8 * gh1_90[k]
                  + pb_x[k] * gi_63[k];

        t_64[k] = f_7 * gh0_109[k]
                  - f_8 * gh1_91[k]
                  + pb_x[k] * gi_64[k];

        t_65[k] = f_5 * gh0_110[k]
                  - f_6 * gh1_92[k]
                  + pb_x[k] * gi_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, gh0_111, gh0_112, gh0_113, gh1_93, gh1_94, \
                         gh1_95, gi_66, gi_67, gi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * gh0_111[k]
                  - f_6 * gh1_93[k]
                  + pb_x[k] * gi_66[k];

        t_67[k] = f_5 * gh0_112[k]
                  - f_6 * gh1_94[k]
                  + pb_x[k] * gi_67[k];

        t_68[k] = f_3 * gh0_113[k]
                  - f_4 * gh1_95[k]
                  + pb_x[k] * gi_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, gh0_114, gh0_115, gh0_117, gh1_96, gh1_97, \
                         gh1_99, gi_69, gi_70, gi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * gh0_114[k]
                  - f_4 * gh1_96[k]
                  + pb_x[k] * gi_69[k];

        t_70[k] = f_3 * gh0_115[k]
                  - f_4 * gh1_97[k]
                  + pb_x[k] * gi_70[k];

        t_71[k] = f_3 * gh0_117[k]
                  - f_4 * gh1_99[k]
                  + pb_x[k] * gi_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, gh0_113, gh0_114, gh0_115, gh1_95, gh1_96, \
                         gh1_97, gi_72, gi_73, gi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * gh0_113[k]
                  - f_2 * gh1_95[k]
                  + pb_y[k] * gi_72[k];

        t_73[k] = f_9 * gh0_114[k]
                  - f_10 * gh1_96[k]
                  + pb_y[k] * gi_73[k];

        t_74[k] = f_7 * gh0_115[k]
                  - f_8 * gh1_97[k]
                  + pb_y[k] * gi_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_y, pb_z, fi_53, gh0_116, gh0_117, gh1_98, \
                         gh1_99, gi_75, gi_76, gi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * gh0_116[k]
                  - f_6 * gh1_98[k]
                  + pb_y[k] * gi_75[k];

        t_76[k] = f_3 * gh0_117[k]
                  - f_4 * gh1_99[k]
                  + pb_y[k] * gi_76[k];

        t_77[k] = f_0 * fi_53[k]
                  + f_1 * gh0_117[k]
                  - f_2 * gh1_99[k]
                  + pb_z[k] * gi_77[k];
    }
}

auto
compute_prim_gk_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_90 = buffer.data(gh1 + 90);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_152 = buffer.data(gi + 152);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_17, \
                         gh1_19, gh1_21, gi_33, gi_35, gi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_17[k]
                 + pb_x[k] * gi_33[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_19[k]
                 + pb_x[k] * gi_35[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_21[k]
                 + pb_x[k] * gi_37[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_1, dk1_5, fi_6, fi_7, fk_8, gh0_7, \
                         gh1_22, gi_39, gi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_22[k]
                 + pb_x[k] * gi_39[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_40[k];

        t_9[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_9, \
                         gh0_10, gh1_30, gh1_32, gi_51, gi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_30[k]
                  + pb_x[k] * gi_51[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_32[k]
                  + pb_x[k] * gi_53[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_11, gh0_12, gh1_33, \
                         gh1_38, gi_55, gi_56, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_33[k]
                  + pb_x[k] * gi_55[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_38[k]
                  + pb_x[k] * gi_56[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_62[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dk0_2, dk1_14, fi_13, fk_14, \
                         fk_15, fk_23, gh0_15, gh1_55, gi_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = pa_x[k] * fk_15[k];

        t_18[k] = pa_x[k] * fk_23[k];

        t_19[k] = f_0 * fi_13[k]
                  + f_1 * gh0_15[k]
                  - f_2 * gh1_55[k]
                  + pb_y[k] * gi_91[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dk0_1, dk1_5, fi_15, fk_15, fk_16, \
                         gh0_18, gh1_70, gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fk_15[k];

        t_21[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_16[k];

        t_22[k] = f_5 * fi_15[k]
                  + f_6 * gh0_18[k]
                  - f_7 * gh1_70[k]
                  + pb_y[k] * gi_117[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fi_16, fi_17, fi_18, gh0_19, gh0_20, gh0_21, \
                         gh1_71, gh1_72, gh1_73, gi_118, gi_119, \
                         gi_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_16[k]
                  + f_8 * gh0_19[k]
                  - f_9 * gh1_71[k]
                  + pb_y[k] * gi_118[k];

        t_24[k] = f_5 * fi_17[k]
                  + f_10 * gh0_20[k]
                  - f_11 * gh1_72[k]
                  + pb_y[k] * gi_119[k];

        t_25[k] = f_5 * fi_18[k]
                  + f_12 * gh0_21[k]
                  - f_13 * gh1_73[k]
                  + pb_y[k] * gi_120[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dk0_2, dk1_14, fi_19, fk_22, fk_23, \
                         gi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fi_19[k]
                  + pb_y[k] * gi_121[k];

        t_27[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_22[k];

        t_28[k] = pa_y[k] * fk_23[k];
    }

#pragma omp simd aligned(t_29, pb_z, fi_20, gh0_23, gh1_90, gi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fi_20[k]
                  + f_1 * gh0_23[k]
                  - f_2 * gh1_90[k]
                  + pb_z[k] * gi_152[k];
    }
}

auto
compute_prim_gk_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_11 = 1.0 / p;
    const auto f_12 = 1.5 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

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
    const auto *fi_17 = buffer.data(fi + 17);
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
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
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
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);

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
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_47 = buffer.data(gh0 + 47);
    const auto *gh0_48 = buffer.data(gh0 + 48);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_78 = buffer.data(gh0 + 78);
    const auto *gh0_79 = buffer.data(gh0 + 79);
    const auto *gh0_80 = buffer.data(gh0 + 80);
    const auto *gh0_81 = buffer.data(gh0 + 81);
    const auto *gh0_82 = buffer.data(gh0 + 82);
    const auto *gh0_83 = buffer.data(gh0 + 83);
    const auto *gh0_84 = buffer.data(gh0 + 84);
    const auto *gh0_85 = buffer.data(gh0 + 85);
    const auto *gh0_86 = buffer.data(gh0 + 86);
    const auto *gh0_87 = buffer.data(gh0 + 87);
    const auto *gh0_88 = buffer.data(gh0 + 88);
    const auto *gh0_89 = buffer.data(gh0 + 89);
    const auto *gh0_90 = buffer.data(gh0 + 90);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_1 = buffer.data(gh1 + 1);
    const auto *gh1_2 = buffer.data(gh1 + 2);
    const auto *gh1_3 = buffer.data(gh1 + 3);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_62 = buffer.data(gh1 + 62);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_74 = buffer.data(gh1 + 74);
    const auto *gh1_75 = buffer.data(gh1 + 75);
    const auto *gh1_76 = buffer.data(gh1 + 76);
    const auto *gh1_77 = buffer.data(gh1 + 77);
    const auto *gh1_83 = buffer.data(gh1 + 83);
    const auto *gh1_85 = buffer.data(gh1 + 85);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_87 = buffer.data(gh1 + 87);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_89 = buffer.data(gh1 + 89);
    const auto *gh1_90 = buffer.data(gh1 + 90);
    const auto *gh1_91 = buffer.data(gh1 + 91);
    const auto *gh1_92 = buffer.data(gh1 + 92);
    const auto *gh1_93 = buffer.data(gh1 + 93);
    const auto *gh1_94 = buffer.data(gh1 + 94);
    const auto *gh1_95 = buffer.data(gh1 + 95);
    const auto *gh1_96 = buffer.data(gh1 + 96);

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
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_21 = buffer.data(gi + 21);
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
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
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
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, gi_0, \
                         gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, \
                         gh1_2, gh1_3, gi_3, gi_4, gi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_6[k] = pb_z[k] * gi_3[k];

        t_7[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_8[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];

        t_9[k] = pb_z[k] * gi_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, gh0_4, gh0_5, gh1_4, gh1_5, gi_6, \
                         gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * gh0_4[k]
                  - f_4 * gh1_4[k]
                  + pb_y[k] * gi_6[k];

        t_11[k] = f_7 * gh0_4[k]
                  - f_8 * gh1_4[k]
                  + pb_z[k] * gi_7[k];

        t_12[k] = f_9 * gh0_5[k]
                  - f_10 * gh1_5[k]
                  + pb_y[k] * gi_8[k];

        t_13[k] = pb_z[k] * gi_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, pb_z, gh0_6, gh0_7, gh0_8, gh1_6, \
                         gh1_7, gh1_8, gi_9, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * gh0_6[k]
                  - f_6 * gh1_6[k]
                  + pb_y[k] * gi_9[k];

        t_15[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_16[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];

        t_17[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, gh0_9, gh0_10, gh0_11, gh1_10, \
                         gh1_11, gh1_12, gi_12, gi_14, gi_15, gi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_z[k] * gi_12[k];

        t_19[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_10[k]
                  + pb_y[k] * gi_14[k];

        t_20[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_11[k]
                  + pb_y[k] * gi_15[k];

        t_21[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_12[k]
                  + pb_y[k] * gi_16[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pb_y, pb_z, fi_1, fk_0, fk_1, gh0_12, \
                         gh1_13, gi_17, gi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_13[k]
                  + pb_y[k] * gi_17[k];

        t_23[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_13[k]
                  + pb_z[k] * gi_18[k];

        t_24[k] = pa_y[k] * fk_0[k];

        t_25[k] = f_11 * fi_1[k]
                  + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, fi_3, fi_5, fi_7, fi_9, \
                         fk_0, fk_3, fk_5, fk_7, fk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_12 * fi_3[k]
                  + pa_y[k] * fk_3[k];

        t_27[k] = f_0 * fi_5[k]
                  + pa_y[k] * fk_5[k];

        t_28[k] = f_13 * fi_7[k]
                  + pa_y[k] * fk_7[k];

        t_29[k] = f_14 * fi_9[k]
                  + pa_y[k] * fk_9[k];

        t_30[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_z, pb_z, fi_0, fi_2, fi_4, fi_6, fk_2, \
                         fk_4, fk_6, gi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_15 * fi_0[k]
                  + pb_z[k] * gi_21[k];

        t_32[k] = f_11 * fi_2[k]
                  + pa_z[k] * fk_2[k];

        t_33[k] = f_12 * fi_4[k]
                  + pa_z[k] * fk_4[k];

        t_34[k] = f_0 * fi_6[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, fi_8, fi_10, fi_11, fi_12, fi_13, \
                         fk_8, fk_10, fk_11, fk_12, fk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * fi_8[k]
                  + pa_z[k] * fk_8[k];

        t_36[k] = f_11 * fi_10[k]
                  + pa_z[k] * fk_10[k];

        t_37[k] = f_12 * fi_11[k]
                  + pa_z[k] * fk_11[k];

        t_38[k] = f_0 * fi_12[k]
                  + pa_z[k] * fk_12[k];

        t_39[k] = f_13 * fi_13[k]
                  + pa_z[k] * fk_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_14, fi_20, \
                         fk_14, fk_15, gh0_17, gh1_19, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_14 * fi_14[k]
                  + pa_z[k] * fk_14[k];

        t_41[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_y[k] * fk_15[k];

        t_42[k] = f_11 * fi_20[k]
                  + f_9 * gh0_17[k]
                  - f_10 * gh1_19[k]
                  + pb_x[k] * gi_25[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_z, fi_21, gh0_15, gh0_16, gh0_19, gh1_17, \
                         gh1_18, gh1_21, gi_24, gi_26, gi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gh0_15[k]
                  - f_4 * gh1_17[k]
                  + pb_z[k] * gi_24[k];

        t_44[k] = f_11 * fi_21[k]
                  + f_7 * gh0_19[k]
                  - f_8 * gh1_21[k]
                  + pb_x[k] * gi_27[k];

        t_45[k] = f_5 * gh0_16[k]
                  - f_6 * gh1_18[k]
                  + pb_z[k] * gi_26[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_z, fi_22, fi_23, gh0_18, gh0_21, gh0_22, \
                         gh1_20, gh1_23, gh1_24, gi_28, gi_29, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_11 * fi_22[k]
                  + f_5 * gh0_21[k]
                  - f_6 * gh1_23[k]
                  + pb_x[k] * gi_29[k];

        t_47[k] = f_7 * gh0_18[k]
                  - f_8 * gh1_20[k]
                  + pb_z[k] * gi_28[k];

        t_48[k] = f_11 * fi_23[k]
                  + f_3 * gh0_22[k]
                  - f_4 * gh1_24[k]
                  + pb_x[k] * gi_31[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_x, pb_x, pb_z, dk0_5, dk1_5, fi_24, fk_22, \
                         gh0_20, gh1_22, gi_30, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * gh0_20[k]
                  - f_10 * gh1_22[k]
                  + pb_z[k] * gi_30[k];

        t_50[k] = f_11 * fi_24[k]
                  + pb_x[k] * gi_32[k];

        t_51[k] = f_16 * dk0_5[k]
                  - f_17 * dk1_5[k]
                  + pa_x[k] * fk_22[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_z, gh0_22, gh0_23, gh0_24, gh1_24, gh1_25, \
                         gh1_26, gi_33, gi_34, gi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * gh0_22[k]
                  - f_4 * gh1_24[k]
                  + pb_z[k] * gi_33[k];

        t_53[k] = f_5 * gh0_23[k]
                  - f_6 * gh1_25[k]
                  + pb_z[k] * gi_34[k];

        t_54[k] = f_7 * gh0_24[k]
                  - f_8 * gh1_26[k]
                  + pb_z[k] * gi_35[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_z, dk0_0, dk1_0, fk_16, gh0_25, gh0_26, \
                         gh1_27, gh1_28, gi_36, gi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_9 * gh0_25[k]
                  - f_10 * gh1_27[k]
                  + pb_z[k] * gi_36[k];

        t_56[k] = f_1 * gh0_26[k]
                  - f_2 * gh1_28[k]
                  + pb_z[k] * gi_37[k];

        t_57[k] = f_16 * dk0_0[k]
                  - f_17 * dk1_0[k]
                  + pa_z[k] * fk_16[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, pb_z, fi_17, fi_26, gh0_27, gh0_30, \
                         gh1_29, gh1_32, gi_38, gi_39, gi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * fi_17[k]
                  + pb_z[k] * gi_38[k];

        t_59[k] = f_3 * gh0_27[k]
                  - f_4 * gh1_29[k]
                  + pb_y[k] * gi_39[k];

        t_60[k] = f_11 * fi_26[k]
                  + f_9 * gh0_30[k]
                  - f_10 * gh1_32[k]
                  + pb_x[k] * gi_41[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, pb_y, fi_27, gh0_28, gh0_29, gh0_32, gh1_30, \
                         gh1_31, gh1_34, gi_40, gi_42, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * gh0_28[k]
                  - f_6 * gh1_30[k]
                  + pb_y[k] * gi_40[k];

        t_62[k] = f_11 * fi_27[k]
                  + f_7 * gh0_32[k]
                  - f_8 * gh1_34[k]
                  + pb_x[k] * gi_43[k];

        t_63[k] = f_7 * gh0_29[k]
                  - f_8 * gh1_31[k]
                  + pb_y[k] * gi_42[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, fi_28, fi_29, gh0_31, gh0_33, gh0_38, \
                         gh1_33, gh1_35, gh1_40, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * fi_28[k]
                  + f_5 * gh0_33[k]
                  - f_6 * gh1_35[k]
                  + pb_x[k] * gi_45[k];

        t_65[k] = f_9 * gh0_31[k]
                  - f_10 * gh1_33[k]
                  + pb_y[k] * gi_44[k];

        t_66[k] = f_11 * fi_29[k]
                  + f_3 * gh0_38[k]
                  - f_4 * gh1_40[k]
                  + pb_x[k] * gi_46[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pb_y, fi_30, gh0_34, gh0_35, gh1_36, gh1_37, \
                         gi_47, gi_48, gi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * fi_30[k]
                  + pb_x[k] * gi_52[k];

        t_68[k] = f_1 * gh0_34[k]
                  - f_2 * gh1_36[k]
                  + pb_y[k] * gi_47[k];

        t_69[k] = f_9 * gh0_35[k]
                  - f_10 * gh1_37[k]
                  + pb_y[k] * gi_48[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_y, gh0_36, gh0_37, gh0_38, gh1_38, gh1_39, \
                         gh1_40, gi_49, gi_50, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_7 * gh0_36[k]
                  - f_8 * gh1_38[k]
                  + pb_y[k] * gi_49[k];

        t_71[k] = f_5 * gh0_37[k]
                  - f_6 * gh1_39[k]
                  + pb_y[k] * gi_50[k];

        t_72[k] = f_3 * gh0_38[k]
                  - f_4 * gh1_40[k]
                  + pb_y[k] * gi_51[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, dk0_14, dk1_14, fi_31, fi_33, fi_35, \
                         fk_28, fk_29, fk_30, fk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_16 * dk0_14[k]
                  - f_17 * dk1_14[k]
                  + pa_x[k] * fk_28[k];

        t_74[k] = f_14 * fi_31[k]
                  + pa_x[k] * fk_29[k];

        t_75[k] = f_13 * fi_33[k]
                  + pa_x[k] * fk_30[k];

        t_76[k] = f_0 * fi_35[k]
                  + pa_x[k] * fk_32[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_x, pb_x, fi_37, fi_39, fi_40, fi_54, \
                         fk_34, fk_36, fk_38, fk_51, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_12 * fi_37[k]
                  + pa_x[k] * fk_34[k];

        t_78[k] = f_11 * fi_39[k]
                  + pa_x[k] * fk_36[k];

        t_79[k] = f_15 * fi_40[k]
                  + pb_x[k] * gi_58[k];

        t_80[k] = pa_x[k] * fk_38[k];

        t_81[k] = f_14 * fi_54[k]
                  + pa_x[k] * fk_51[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pb_z, fi_25, fi_57, fi_59, fi_61, \
                         fk_53, fk_55, fk_57, gi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_12 * fi_25[k]
                  + pb_z[k] * gi_59[k];

        t_83[k] = f_13 * fi_57[k]
                  + pa_x[k] * fk_53[k];

        t_84[k] = f_0 * fi_59[k]
                  + pa_x[k] * fk_55[k];

        t_85[k] = f_12 * fi_61[k]
                  + pa_x[k] * fk_57[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, fi_62, fi_68, fk_59, \
                         fk_65, gh0_47, gh1_49, gi_64, gi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * fi_62[k]
                  + pa_x[k] * fk_59[k];

        t_87[k] = f_15 * fi_68[k]
                  + pb_x[k] * gi_64[k];

        t_88[k] = pa_x[k] * fk_65[k];

        t_89[k] = f_1 * gh0_47[k]
                  - f_2 * gh1_49[k]
                  + pb_x[k] * gi_65[k];

        t_90[k] = pb_z[k] * gi_65[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pb_x, pb_z, gh0_48, gh0_49, gh0_50, gh1_51, \
                         gh1_52, gh1_53, gi_67, gi_68, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * gh0_48[k]
                  - f_10 * gh1_51[k]
                  + pb_x[k] * gi_67[k];

        t_92[k] = f_9 * gh0_49[k]
                  - f_10 * gh1_52[k]
                  + pb_x[k] * gi_68[k];

        t_93[k] = f_7 * gh0_50[k]
                  - f_8 * gh1_53[k]
                  + pb_x[k] * gi_69[k];

        t_94[k] = pb_z[k] * gi_67[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_x, pb_z, gh0_51, gh0_52, gh0_53, gh1_54, \
                         gh1_55, gh1_56, gi_69, gi_70, gi_71, gi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_7 * gh0_51[k]
                  - f_8 * gh1_54[k]
                  + pb_x[k] * gi_70[k];

        t_96[k] = f_5 * gh0_52[k]
                  - f_6 * gh1_55[k]
                  + pb_x[k] * gi_71[k];

        t_97[k] = pb_z[k] * gi_69[k];

        t_98[k] = f_5 * gh0_53[k]
                  - f_6 * gh1_56[k]
                  + pb_x[k] * gi_72[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pb_x, pb_z, gh0_54, gh0_55, gh0_57, \
                         gh1_57, gh1_58, gh1_60, gi_71, gi_73, gi_74, \
                         gi_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_5 * gh0_54[k]
                  - f_6 * gh1_57[k]
                  + pb_x[k] * gi_73[k];

        t_100[k] = f_3 * gh0_55[k]
                   - f_4 * gh1_58[k]
                   + pb_x[k] * gi_74[k];

        t_101[k] = pb_z[k] * gi_71[k];

        t_102[k] = f_3 * gh0_57[k]
                   - f_4 * gh1_60[k]
                   + pb_x[k] * gi_75[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pb_x, pb_y, fi_40, gh0_55, gh0_58, \
                         gh0_59, gh1_58, gh1_61, gh1_62, gi_76, gi_77, \
                         gi_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_3 * gh0_58[k]
                   - f_4 * gh1_61[k]
                   + pb_x[k] * gi_76[k];

        t_104[k] = f_3 * gh0_59[k]
                   - f_4 * gh1_62[k]
                   + pb_x[k] * gi_77[k];

        t_105[k] = pb_x[k] * gi_78[k];

        t_106[k] = f_0 * fi_40[k]
                   + f_1 * gh0_55[k]
                   - f_2 * gh1_58[k]
                   + pb_y[k] * gi_78[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, gh0_55, gh0_56, gh0_57, gh1_58, \
                         gh1_59, gh1_60, gi_78, gi_79, gi_80, gi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pb_z[k] * gi_78[k];

        t_108[k] = f_3 * gh0_55[k]
                   - f_4 * gh1_58[k]
                   + pb_z[k] * gi_79[k];

        t_109[k] = f_5 * gh0_56[k]
                   - f_6 * gh1_59[k]
                   + pb_z[k] * gi_80[k];

        t_110[k] = f_7 * gh0_57[k]
                   - f_8 * gh1_60[k]
                   + pb_z[k] * gi_81[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_z, pb_y, pb_z, fi_32, fi_45, fk_31, \
                         gh0_58, gh0_59, gh1_61, gh1_62, gi_82, gi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_9 * gh0_58[k]
                   - f_10 * gh1_61[k]
                   + pb_z[k] * gi_82[k];

        t_112[k] = f_0 * fi_45[k]
                   + pb_y[k] * gi_83[k];

        t_113[k] = f_1 * gh0_59[k]
                   - f_2 * gh1_62[k]
                   + pb_z[k] * gi_83[k];

        t_114[k] = f_11 * fi_32[k]
                   + pa_z[k] * fk_31[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_z, pb_z, fi_34, fi_36, fi_38, \
                         fi_40, fk_33, fk_35, fk_37, fk_38, gi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_12 * fi_34[k]
                   + pa_z[k] * fk_33[k];

        t_116[k] = f_0 * fi_36[k]
                   + pa_z[k] * fk_35[k];

        t_117[k] = f_13 * fi_38[k]
                   + pa_z[k] * fk_37[k];

        t_118[k] = pa_z[k] * fk_38[k];

        t_119[k] = f_15 * fi_40[k]
                   + pb_z[k] * gi_84[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, fi_41, fi_42, fi_43, fi_44, fk_39, \
                         fk_40, fk_41, fk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * fi_41[k]
                   + pa_z[k] * fk_39[k];

        t_121[k] = f_12 * fi_42[k]
                   + pa_z[k] * fk_40[k];

        t_122[k] = f_0 * fi_43[k]
                   + pa_z[k] * fk_41[k];

        t_123[k] = f_13 * fi_44[k]
                   + pa_z[k] * fk_42[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_z, pb_x, pb_y, fi_45, fi_47, fk_43, gh0_61, \
                         gh1_65, gi_85, gi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_12 * fi_47[k]
                   + pb_y[k] * gi_85[k];

        t_125[k] = f_14 * fi_45[k]
                   + pa_z[k] * fk_43[k];

        t_126[k] = f_1 * gh0_61[k]
                   - f_2 * gh1_65[k]
                   + pb_x[k] * gi_86[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pb_x, gh0_62, gh0_63, gh0_64, gh1_66, gh1_67, \
                         gh1_68, gi_87, gi_88, gi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_9 * gh0_62[k]
                   - f_10 * gh1_66[k]
                   + pb_x[k] * gi_87[k];

        t_128[k] = f_9 * gh0_63[k]
                   - f_10 * gh1_67[k]
                   + pb_x[k] * gi_88[k];

        t_129[k] = f_7 * gh0_64[k]
                   - f_8 * gh1_68[k]
                   + pb_x[k] * gi_89[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pb_x, gh0_65, gh0_66, gh0_67, gh1_69, gh1_70, \
                         gh1_71, gi_90, gi_91, gi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * gh0_65[k]
                   - f_8 * gh1_69[k]
                   + pb_x[k] * gi_90[k];

        t_131[k] = f_5 * gh0_66[k]
                   - f_6 * gh1_70[k]
                   + pb_x[k] * gi_91[k];

        t_132[k] = f_5 * gh0_67[k]
                   - f_6 * gh1_71[k]
                   + pb_x[k] * gi_92[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_x, gh0_68, gh0_69, gh0_70, gh1_72, gh1_73, \
                         gh1_74, gi_93, gi_94, gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * gh0_68[k]
                   - f_6 * gh1_72[k]
                   + pb_x[k] * gi_93[k];

        t_134[k] = f_3 * gh0_69[k]
                   - f_4 * gh1_73[k]
                   + pb_x[k] * gi_94[k];

        t_135[k] = f_3 * gh0_70[k]
                   - f_4 * gh1_74[k]
                   + pb_x[k] * gi_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_z, pb_x, dk0_5, dk1_5, fk_44, gh0_71, gh0_73, \
                         gh1_75, gh1_77, gi_96, gi_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_3 * gh0_71[k]
                   - f_4 * gh1_75[k]
                   + pb_x[k] * gi_96[k];

        t_137[k] = f_3 * gh0_73[k]
                   - f_4 * gh1_77[k]
                   + pb_x[k] * gi_97[k];

        t_138[k] = f_16 * dk0_5[k]
                   - f_17 * dk1_5[k]
                   + pa_z[k] * fk_44[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pb_y, pb_z, fi_46, fi_49, fi_50, gh0_70, gh0_71, \
                         gh1_74, gh1_75, gi_98, gi_99, gi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_11 * fi_46[k]
                   + pb_z[k] * gi_98[k];

        t_140[k] = f_11 * fi_49[k]
                   + f_9 * gh0_70[k]
                   - f_10 * gh1_74[k]
                   + pb_y[k] * gi_99[k];

        t_141[k] = f_11 * fi_50[k]
                   + f_7 * gh0_71[k]
                   - f_8 * gh1_75[k]
                   + pb_y[k] * gi_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_y, fi_51, fi_52, fi_53, gh0_72, gh0_73, \
                         gh1_76, gh1_77, gi_101, gi_102, gi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_11 * fi_51[k]
                   + f_5 * gh0_72[k]
                   - f_6 * gh1_76[k]
                   + pb_y[k] * gi_101[k];

        t_143[k] = f_11 * fi_52[k]
                   + f_3 * gh0_73[k]
                   - f_4 * gh1_77[k]
                   + pb_y[k] * gi_102[k];

        t_144[k] = f_11 * fi_53[k]
                   + pb_y[k] * gi_103[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_y, dk0_14, dk1_14, fi_55, fi_56, \
                         fi_58, fk_50, fk_52, fk_54, fk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_16 * dk0_14[k]
                   - f_17 * dk1_14[k]
                   + pa_y[k] * fk_50[k];

        t_146[k] = f_11 * fi_55[k]
                   + pa_y[k] * fk_52[k];

        t_147[k] = f_12 * fi_56[k]
                   + pa_y[k] * fk_54[k];

        t_148[k] = f_0 * fi_58[k]
                   + pa_y[k] * fk_56[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_y, pb_z, fi_48, fi_60, fi_63, fi_64, \
                         fk_58, fk_60, fk_61, gi_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_13 * fi_60[k]
                   + pa_y[k] * fk_58[k];

        t_150[k] = f_14 * fi_63[k]
                   + pa_y[k] * fk_60[k];

        t_151[k] = f_12 * fi_48[k]
                   + pb_z[k] * gi_104[k];

        t_152[k] = f_13 * fi_64[k]
                   + pa_y[k] * fk_61[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pa_y, pb_y, fi_65, fi_66, fi_67, \
                         fi_68, fk_62, fk_63, fk_64, fk_65, gi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * fi_65[k]
                   + pa_y[k] * fk_62[k];

        t_154[k] = f_12 * fi_66[k]
                   + pa_y[k] * fk_63[k];

        t_155[k] = f_11 * fi_67[k]
                   + pa_y[k] * fk_64[k];

        t_156[k] = f_15 * fi_68[k]
                   + pb_y[k] * gi_109[k];

        t_157[k] = pa_y[k] * fk_65[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_x, pb_y, pb_z, fi_54, gh0_78, gh0_79, \
                         gh1_83, gh1_85, gi_110, gi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_1 * gh0_78[k]
                   - f_2 * gh1_83[k]
                   + pb_x[k] * gi_110[k];

        t_159[k] = pb_y[k] * gi_110[k];

        t_160[k] = f_0 * fi_54[k]
                   + pb_z[k] * gi_110[k];

        t_161[k] = f_9 * gh0_79[k]
                   - f_10 * gh1_85[k]
                   + pb_x[k] * gi_112[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_x, pb_y, gh0_80, gh0_81, gh0_82, \
                         gh1_86, gh1_87, gh1_88, gi_113, gi_114, \
                         gi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_9 * gh0_80[k]
                   - f_10 * gh1_86[k]
                   + pb_x[k] * gi_113[k];

        t_163[k] = f_7 * gh0_81[k]
                   - f_8 * gh1_87[k]
                   + pb_x[k] * gi_114[k];

        t_164[k] = pb_y[k] * gi_113[k];

        t_165[k] = f_7 * gh0_82[k]
                   - f_8 * gh1_88[k]
                   + pb_x[k] * gi_115[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_x, pb_y, gh0_83, gh0_84, gh0_85, \
                         gh1_89, gh1_90, gh1_91, gi_115, gi_116, gi_117, \
                         gi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_5 * gh0_83[k]
                   - f_6 * gh1_89[k]
                   + pb_x[k] * gi_116[k];

        t_167[k] = f_5 * gh0_84[k]
                   - f_6 * gh1_90[k]
                   + pb_x[k] * gi_117[k];

        t_168[k] = pb_y[k] * gi_115[k];

        t_169[k] = f_5 * gh0_85[k]
                   - f_6 * gh1_91[k]
                   + pb_x[k] * gi_118[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_x, pb_y, gh0_86, gh0_87, gh0_88, \
                         gh1_92, gh1_93, gh1_94, gi_118, gi_119, gi_120, \
                         gi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_3 * gh0_86[k]
                   - f_4 * gh1_92[k]
                   + pb_x[k] * gi_119[k];

        t_171[k] = f_3 * gh0_87[k]
                   - f_4 * gh1_93[k]
                   + pb_x[k] * gi_120[k];

        t_172[k] = f_3 * gh0_88[k]
                   - f_4 * gh1_94[k]
                   + pb_x[k] * gi_121[k];

        t_173[k] = pb_y[k] * gi_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_x, pb_y, pb_z, fi_63, gh0_86, gh0_90, \
                         gh1_92, gh1_96, gi_122, gi_123, gi_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_3 * gh0_90[k]
                   - f_4 * gh1_96[k]
                   + pb_x[k] * gi_122[k];

        t_175[k] = pb_x[k] * gi_128[k];

        t_176[k] = f_1 * gh0_86[k]
                   - f_2 * gh1_92[k]
                   + pb_y[k] * gi_123[k];

        t_177[k] = f_0 * fi_63[k]
                   + pb_z[k] * gi_123[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pb_y, gh0_87, gh0_88, gh0_89, gh1_93, gh1_94, \
                         gh1_95, gi_124, gi_125, gi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_9 * gh0_87[k]
                   - f_10 * gh1_93[k]
                   + pb_y[k] * gi_124[k];

        t_179[k] = f_7 * gh0_88[k]
                   - f_8 * gh1_94[k]
                   + pb_y[k] * gi_125[k];

        t_180[k] = f_5 * gh0_89[k]
                   - f_6 * gh1_95[k]
                   + pb_y[k] * gi_126[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, pb_y, pb_z, fi_68, gh0_90, gh1_96, gi_127, \
                         gi_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_3 * gh0_90[k]
                   - f_4 * gh1_96[k]
                   + pb_y[k] * gi_127[k];

        t_182[k] = pb_y[k] * gi_128[k];

        t_183[k] = f_0 * fi_68[k]
                   + f_1 * gh0_90[k]
                   - f_2 * gh1_96[k]
                   + pb_z[k] * gi_128[k];
    }
}

auto
compute_prim_gk_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);
    const auto f_13 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_62 = buffer.data(fi + 62);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_1 = buffer.data(gh0 + 1);
    const auto *gh0_2 = buffer.data(gh0 + 2);
    const auto *gh0_3 = buffer.data(gh0 + 3);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_49 = buffer.data(gh0 + 49);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_58 = buffer.data(gh0 + 58);
    const auto *gh0_59 = buffer.data(gh0 + 59);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_61 = buffer.data(gh0 + 61);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_74 = buffer.data(gh0 + 74);
    const auto *gh0_75 = buffer.data(gh0 + 75);
    const auto *gh0_76 = buffer.data(gh0 + 76);
    const auto *gh0_77 = buffer.data(gh0 + 77);
    const auto *gh0_83 = buffer.data(gh0 + 83);
    const auto *gh0_85 = buffer.data(gh0 + 85);
    const auto *gh0_86 = buffer.data(gh0 + 86);
    const auto *gh0_87 = buffer.data(gh0 + 87);
    const auto *gh0_88 = buffer.data(gh0 + 88);
    const auto *gh0_89 = buffer.data(gh0 + 89);
    const auto *gh0_90 = buffer.data(gh0 + 90);
    const auto *gh0_91 = buffer.data(gh0 + 91);
    const auto *gh0_92 = buffer.data(gh0 + 92);
    const auto *gh0_93 = buffer.data(gh0 + 93);
    const auto *gh0_94 = buffer.data(gh0 + 94);
    const auto *gh0_95 = buffer.data(gh0 + 95);
    const auto *gh0_96 = buffer.data(gh0 + 96);

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
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_47 = buffer.data(gh1 + 47);
    const auto *gh1_48 = buffer.data(gh1 + 48);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_71 = buffer.data(gh1 + 71);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_78 = buffer.data(gh1 + 78);
    const auto *gh1_79 = buffer.data(gh1 + 79);
    const auto *gh1_80 = buffer.data(gh1 + 80);
    const auto *gh1_81 = buffer.data(gh1 + 81);
    const auto *gh1_82 = buffer.data(gh1 + 82);
    const auto *gh1_83 = buffer.data(gh1 + 83);
    const auto *gh1_84 = buffer.data(gh1 + 84);
    const auto *gh1_85 = buffer.data(gh1 + 85);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_87 = buffer.data(gh1 + 87);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_89 = buffer.data(gh1 + 89);
    const auto *gh1_90 = buffer.data(gh1 + 90);

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
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
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
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, gi_0, \
                         gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, gh1_2, \
                         gh1_3, gi_3, gi_4, gi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_6[k] = pb_y[k] * gi_4[k];

        t_7[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_8[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, gh0_4, gh0_5, gh1_4, gh1_5, gi_6, \
                         gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * gh0_4[k]
                 - f_4 * gh1_4[k]
                 + pb_y[k] * gi_6[k];

        t_10[k] = pb_y[k] * gi_7[k];

        t_11[k] = f_7 * gh0_4[k]
                  - f_8 * gh1_4[k]
                  + pb_z[k] * gi_7[k];

        t_12[k] = f_9 * gh0_5[k]
                  - f_10 * gh1_5[k]
                  + pb_y[k] * gi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, gh0_6, gh0_7, gh1_6, gh1_7, gi_9, \
                         gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gh0_6[k]
                  - f_6 * gh1_6[k]
                  + pb_y[k] * gi_9[k];

        t_14[k] = f_3 * gh0_7[k]
                  - f_4 * gh1_7[k]
                  + pb_y[k] * gi_10[k];

        t_15[k] = pb_y[k] * gi_11[k];

        t_16[k] = f_9 * gh0_7[k]
                  - f_10 * gh1_7[k]
                  + pb_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, gh0_8, gh0_10, gh0_11, gh1_8, gh1_9, gh1_10, \
                         gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * gh0_8[k]
                  - f_2 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_18[k] = f_9 * gh0_10[k]
                  - f_10 * gh1_9[k]
                  + pb_y[k] * gi_13[k];

        t_19[k] = f_7 * gh0_11[k]
                  - f_8 * gh1_10[k]
                  + pb_y[k] * gi_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, pb_z, fk_0, gh0_12, gh0_13, \
                         gh1_11, gh1_12, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * gh0_12[k]
                  - f_6 * gh1_11[k]
                  + pb_y[k] * gi_15[k];

        t_21[k] = f_3 * gh0_13[k]
                  - f_4 * gh1_12[k]
                  + pb_y[k] * gi_16[k];

        t_22[k] = pb_y[k] * gi_17[k];

        t_23[k] = f_1 * gh0_13[k]
                  - f_2 * gh1_12[k]
                  + pb_z[k] * gi_17[k];

        t_24[k] = pa_y[k] * fk_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_17, fk_0, fk_1, \
                         gh0_19, gh1_17, gi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * fk_0[k];

        t_26[k] = f_11 * dk0_0[k]
                  - f_12 * dk1_0[k]
                  + pa_y[k] * fk_1[k];

        t_27[k] = f_13 * fi_17[k]
                  + f_9 * gh0_19[k]
                  - f_10 * gh1_17[k]
                  + pb_x[k] * gi_21[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, fi_18, fi_19, fi_20, gh0_21, gh0_23, gh0_24, \
                         gh1_19, gh1_21, gh1_22, gi_22, gi_23, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_13 * fi_18[k]
                  + f_7 * gh0_21[k]
                  - f_8 * gh1_19[k]
                  + pb_x[k] * gi_22[k];

        t_29[k] = f_13 * fi_19[k]
                  + f_5 * gh0_23[k]
                  - f_6 * gh1_21[k]
                  + pb_x[k] * gi_23[k];

        t_30[k] = f_13 * fi_20[k]
                  + f_3 * gh0_24[k]
                  - f_4 * gh1_22[k]
                  + pb_x[k] * gi_24[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pa_z, pb_x, dk0_0, dk0_5, dk1_0, dk1_5, \
                         fi_21, fk_2, fk_3, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * fi_21[k]
                  + pb_x[k] * gi_25[k];

        t_32[k] = f_11 * dk0_5[k]
                  - f_12 * dk1_5[k]
                  + pa_x[k] * fk_3[k];

        t_33[k] = f_11 * dk0_0[k]
                  - f_12 * dk1_0[k]
                  + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, fi_22, fi_23, fi_24, gh0_32, gh0_34, gh0_35, \
                         gh1_30, gh1_32, gh1_33, gi_27, gi_28, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * fi_22[k]
                  + f_9 * gh0_32[k]
                  - f_10 * gh1_30[k]
                  + pb_x[k] * gi_27[k];

        t_35[k] = f_13 * fi_23[k]
                  + f_7 * gh0_34[k]
                  - f_8 * gh1_32[k]
                  + pb_x[k] * gi_28[k];

        t_36[k] = f_13 * fi_24[k]
                  + f_5 * gh0_35[k]
                  - f_6 * gh1_33[k]
                  + pb_x[k] * gi_29[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, dk0_14, dk1_14, fi_25, fi_26, \
                         fk_4, fk_5, gh0_40, gh1_38, gi_30, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * fi_25[k]
                  + f_3 * gh0_40[k]
                  - f_4 * gh1_38[k]
                  + pb_x[k] * gi_30[k];

        t_38[k] = f_13 * fi_26[k]
                  + pb_x[k] * gi_31[k];

        t_39[k] = f_11 * dk0_14[k]
                  - f_12 * dk1_14[k]
                  + pa_x[k] * fk_4[k];

        t_40[k] = pa_x[k] * fk_5[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_x, fk_8, gh0_49, gh0_51, gh0_52, \
                         gh1_47, gh1_48, gh1_49, gi_34, gi_35, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_x[k] * fk_8[k];

        t_42[k] = f_1 * gh0_49[k]
                  - f_2 * gh1_47[k]
                  + pb_x[k] * gi_34[k];

        t_43[k] = f_9 * gh0_51[k]
                  - f_10 * gh1_48[k]
                  + pb_x[k] * gi_35[k];

        t_44[k] = f_9 * gh0_52[k]
                  - f_10 * gh1_49[k]
                  + pb_x[k] * gi_36[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, gh0_53, gh0_54, gh0_55, gh1_50, gh1_51, \
                         gh1_52, gi_37, gi_38, gi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_7 * gh0_53[k]
                  - f_8 * gh1_50[k]
                  + pb_x[k] * gi_37[k];

        t_46[k] = f_7 * gh0_54[k]
                  - f_8 * gh1_51[k]
                  + pb_x[k] * gi_38[k];

        t_47[k] = f_5 * gh0_55[k]
                  - f_6 * gh1_52[k]
                  + pb_x[k] * gi_39[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, gh0_56, gh0_57, gh0_58, gh1_53, gh1_54, \
                         gh1_55, gi_40, gi_41, gi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * gh0_56[k]
                  - f_6 * gh1_53[k]
                  + pb_x[k] * gi_40[k];

        t_49[k] = f_5 * gh0_57[k]
                  - f_6 * gh1_54[k]
                  + pb_x[k] * gi_41[k];

        t_50[k] = f_3 * gh0_58[k]
                  - f_4 * gh1_55[k]
                  + pb_x[k] * gi_42[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, gh0_60, gh0_61, gh0_62, gh1_57, gh1_58, \
                         gh1_59, gi_43, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * gh0_60[k]
                  - f_4 * gh1_57[k]
                  + pb_x[k] * gi_43[k];

        t_52[k] = f_3 * gh0_61[k]
                  - f_4 * gh1_58[k]
                  + pb_x[k] * gi_44[k];

        t_53[k] = f_3 * gh0_62[k]
                  - f_4 * gh1_59[k]
                  + pb_x[k] * gi_45[k];

        t_54[k] = pb_x[k] * gi_46[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pb_x, pb_y, fi_36, gh0_58, gh1_55, \
                         gi_46, gi_48, gi_49, gi_50, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * gi_48[k];

        t_56[k] = pb_x[k] * gi_49[k];

        t_57[k] = pb_x[k] * gi_50[k];

        t_58[k] = pb_x[k] * gi_51[k];

        t_59[k] = f_0 * fi_36[k]
                  + f_1 * gh0_58[k]
                  - f_2 * gh1_55[k]
                  + pb_y[k] * gi_46[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_z, gh0_58, gh0_59, gh0_60, gh1_55, gh1_56, \
                         gh1_57, gi_46, gi_47, gi_48, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pb_z[k] * gi_46[k];

        t_61[k] = f_3 * gh0_58[k]
                  - f_4 * gh1_55[k]
                  + pb_z[k] * gi_47[k];

        t_62[k] = f_5 * gh0_59[k]
                  - f_6 * gh1_56[k]
                  + pb_z[k] * gi_48[k];

        t_63[k] = f_7 * gh0_60[k]
                  - f_8 * gh1_57[k]
                  + pb_z[k] * gi_49[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_z, dk0_5, dk1_5, fk_5, fk_6, gh0_61, \
                         gh0_62, gh1_58, gh1_59, gi_50, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_9 * gh0_61[k]
                  - f_10 * gh1_58[k]
                  + pb_z[k] * gi_50[k];

        t_65[k] = f_1 * gh0_62[k]
                  - f_2 * gh1_59[k]
                  + pb_z[k] * gi_51[k];

        t_66[k] = pa_z[k] * fk_5[k];

        t_67[k] = f_11 * dk0_5[k]
                  - f_12 * dk1_5[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_y, fi_43, fi_44, fi_45, gh0_74, gh0_75, gh0_76, \
                         gh1_70, gh1_71, gh1_72, gi_54, gi_55, gi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_13 * fi_43[k]
                  + f_9 * gh0_74[k]
                  - f_10 * gh1_70[k]
                  + pb_y[k] * gi_54[k];

        t_69[k] = f_13 * fi_44[k]
                  + f_7 * gh0_75[k]
                  - f_8 * gh1_71[k]
                  + pb_y[k] * gi_55[k];

        t_70[k] = f_13 * fi_45[k]
                  + f_5 * gh0_76[k]
                  - f_6 * gh1_72[k]
                  + pb_y[k] * gi_56[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, dk0_14, dk1_14, fi_46, fi_47, \
                         fk_7, fk_8, gh0_77, gh1_73, gi_57, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_13 * fi_46[k]
                  + f_3 * gh0_77[k]
                  - f_4 * gh1_73[k]
                  + pb_y[k] * gi_57[k];

        t_72[k] = f_13 * fi_47[k]
                  + pb_y[k] * gi_58[k];

        t_73[k] = f_11 * dk0_14[k]
                  - f_12 * dk1_14[k]
                  + pa_y[k] * fk_7[k];

        t_74[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, gh0_83, gh0_85, gh0_86, gh1_78, gh1_79, \
                         gh1_80, gi_60, gi_61, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * gh0_83[k]
                  - f_2 * gh1_78[k]
                  + pb_x[k] * gi_60[k];

        t_76[k] = f_9 * gh0_85[k]
                  - f_10 * gh1_79[k]
                  + pb_x[k] * gi_61[k];

        t_77[k] = f_9 * gh0_86[k]
                  - f_10 * gh1_80[k]
                  + pb_x[k] * gi_62[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_x, gh0_87, gh0_88, gh0_89, gh1_81, gh1_82, \
                         gh1_83, gi_63, gi_64, gi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_7 * gh0_87[k]
                  - f_8 * gh1_81[k]
                  + pb_x[k] * gi_63[k];

        t_79[k] = f_7 * gh0_88[k]
                  - f_8 * gh1_82[k]
                  + pb_x[k] * gi_64[k];

        t_80[k] = f_5 * gh0_89[k]
                  - f_6 * gh1_83[k]
                  + pb_x[k] * gi_65[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, gh0_90, gh0_91, gh0_92, gh1_84, gh1_85, \
                         gh1_86, gi_66, gi_67, gi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * gh0_90[k]
                  - f_6 * gh1_84[k]
                  + pb_x[k] * gi_66[k];

        t_82[k] = f_5 * gh0_91[k]
                  - f_6 * gh1_85[k]
                  + pb_x[k] * gi_67[k];

        t_83[k] = f_3 * gh0_92[k]
                  - f_4 * gh1_86[k]
                  + pb_x[k] * gi_68[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, gh0_93, gh0_94, gh0_96, gh1_87, gh1_88, \
                         gh1_90, gi_69, gi_70, gi_71, gi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_3 * gh0_93[k]
                  - f_4 * gh1_87[k]
                  + pb_x[k] * gi_69[k];

        t_85[k] = f_3 * gh0_94[k]
                  - f_4 * gh1_88[k]
                  + pb_x[k] * gi_70[k];

        t_86[k] = f_3 * gh0_96[k]
                  - f_4 * gh1_90[k]
                  + pb_x[k] * gi_71[k];

        t_87[k] = pb_x[k] * gi_72[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pb_x, pb_y, gh0_92, gh1_86, gi_72, \
                         gi_73, gi_74, gi_75, gi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_x[k] * gi_73[k];

        t_89[k] = pb_x[k] * gi_74[k];

        t_90[k] = pb_x[k] * gi_75[k];

        t_91[k] = pb_x[k] * gi_77[k];

        t_92[k] = f_1 * gh0_92[k]
                  - f_2 * gh1_86[k]
                  + pb_y[k] * gi_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_y, gh0_93, gh0_94, gh0_95, gh1_87, gh1_88, \
                         gh1_89, gi_73, gi_74, gi_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_9 * gh0_93[k]
                  - f_10 * gh1_87[k]
                  + pb_y[k] * gi_73[k];

        t_94[k] = f_7 * gh0_94[k]
                  - f_8 * gh1_88[k]
                  + pb_y[k] * gi_74[k];

        t_95[k] = f_5 * gh0_95[k]
                  - f_6 * gh1_89[k]
                  + pb_y[k] * gi_75[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pb_z, fi_62, gh0_96, gh1_90, gi_76, \
                         gi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_3 * gh0_96[k]
                  - f_4 * gh1_90[k]
                  + pb_y[k] * gi_76[k];

        t_97[k] = pb_y[k] * gi_77[k];

        t_98[k] = f_0 * fi_62[k]
                  + f_1 * gh0_96[k]
                  - f_2 * gh1_90[k]
                  + pb_z[k] * gi_77[k];
    }
}

auto
compute_prim_gk_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_29 = buffer.data(fi + 29);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_32 = buffer.data(gh1 + 32);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_38 = buffer.data(gi + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_7, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_9, gh1_10, gh1_11, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_12, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_16, fk_5, fk_6, \
                         gh0_15, gh1_21, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_16[k]
                  + f_1 * gh0_15[k]
                  - f_2 * gh1_21[k]
                  + pb_y[k] * gi_25[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_18, fi_19, fi_20, gh0_18, gh0_19, gh0_20, \
                         gh1_24, gh1_25, gh1_26, gi_28, gi_29, gi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_18[k]
                  + f_6 * gh0_18[k]
                  - f_7 * gh1_24[k]
                  + pb_y[k] * gi_28[k];

        t_21[k] = f_5 * fi_19[k]
                  + f_8 * gh0_19[k]
                  - f_9 * gh1_25[k]
                  + pb_y[k] * gi_29[k];

        t_22[k] = f_5 * fi_20[k]
                  + f_10 * gh0_20[k]
                  - f_11 * gh1_26[k]
                  + pb_y[k] * gi_30[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_21, fk_7, fk_8, \
                         gh0_21, gh1_27, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_21[k]
                  + f_12 * gh0_21[k]
                  - f_13 * gh1_27[k]
                  + pb_y[k] * gi_31[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_29, gh0_23, gh1_32, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_29[k]
                  + f_1 * gh0_23[k]
                  - f_2 * gh1_32[k]
                  + pb_z[k] * gi_38[k];
    }
}

auto
compute_prim_gk_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_29 = buffer.data(fi + 29);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_32 = buffer.data(gh0 + 32);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_56 = buffer.data(gh1 + 56);
    const auto *gh1_57 = buffer.data(gh1 + 57);
    const auto *gh1_74 = buffer.data(gh1 + 74);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_71 = buffer.data(gi + 71);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_17, \
                         gh1_18, gh1_19, gi_15, gi_16, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_17[k]
                 + pb_x[k] * gi_15[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_18[k]
                 + pb_x[k] * gi_16[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_19[k]
                 + pb_x[k] * gi_17[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_20, gi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_20[k]
                 + pb_x[k] * gi_18[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_23, gh1_24, gh1_25, gi_21, gi_22, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_23[k]
                  + pb_x[k] * gi_21[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_24[k]
                  + pb_x[k] * gi_22[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_25[k]
                  + pb_x[k] * gi_23[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_27, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_27[k]
                  + pb_x[k] * gi_24[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_16, fk_5, fk_6, \
                         gh0_21, gh1_45, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_16[k]
                  + f_1 * gh0_21[k]
                  - f_2 * gh1_45[k]
                  + pb_y[k] * gi_43[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_18, fi_19, fi_20, gh0_24, gh0_25, gh0_26, \
                         gh1_54, gh1_55, gh1_56, gi_50, gi_51, gi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_18[k]
                  + f_6 * gh0_24[k]
                  - f_7 * gh1_54[k]
                  + pb_y[k] * gi_50[k];

        t_21[k] = f_5 * fi_19[k]
                  + f_8 * gh0_25[k]
                  - f_9 * gh1_55[k]
                  + pb_y[k] * gi_51[k];

        t_22[k] = f_5 * fi_20[k]
                  + f_10 * gh0_26[k]
                  - f_11 * gh1_56[k]
                  + pb_y[k] * gi_52[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_21, fk_7, fk_8, \
                         gh0_27, gh1_57, gi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_21[k]
                  + f_12 * gh0_27[k]
                  - f_13 * gh1_57[k]
                  + pb_y[k] * gi_53[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_29, gh0_32, gh1_74, gi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_29[k]
                  + f_1 * gh0_32[k]
                  - f_2 * gh1_74[k]
                  + pb_z[k] * gi_71[k];
    }
}

auto
compute_prim_gk_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_29 = buffer.data(fi + 29);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_45 = buffer.data(gh0 + 45);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_56 = buffer.data(gh0 + 56);
    const auto *gh0_57 = buffer.data(gh0 + 57);
    const auto *gh0_74 = buffer.data(gh0 + 74);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);
    const auto *gh1_48 = buffer.data(gh1 + 48);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_65 = buffer.data(gh1 + 65);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_26 = buffer.data(gi + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_17, gh0_18, gh0_19, \
                         gh1_15, gh1_16, gh1_17, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_17[k]
                 - f_7 * gh1_15[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_18[k]
                 - f_9 * gh1_16[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_19[k]
                 - f_11 * gh1_17[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_20, gh1_18, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_20[k]
                 - f_13 * gh1_18[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_23, gh0_24, gh0_25, \
                         gh1_20, gh1_21, gh1_22, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_23[k]
                  - f_7 * gh1_20[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_24[k]
                  - f_9 * gh1_21[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_25[k]
                  - f_11 * gh1_22[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_27, gh1_23, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_27[k]
                  - f_13 * gh1_23[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_16, fk_5, fk_6, \
                         gh0_45, gh1_39, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_16[k]
                  + f_1 * gh0_45[k]
                  - f_2 * gh1_39[k]
                  + pb_y[k] * gi_17[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_18, fi_19, fi_20, gh0_54, gh0_55, gh0_56, \
                         gh1_46, gh1_47, gh1_48, gi_20, gi_21, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_18[k]
                  + f_6 * gh0_54[k]
                  - f_7 * gh1_46[k]
                  + pb_y[k] * gi_20[k];

        t_21[k] = f_5 * fi_19[k]
                  + f_8 * gh0_55[k]
                  - f_9 * gh1_47[k]
                  + pb_y[k] * gi_21[k];

        t_22[k] = f_5 * fi_20[k]
                  + f_10 * gh0_56[k]
                  - f_11 * gh1_48[k]
                  + pb_y[k] * gi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_21, fk_7, fk_8, \
                         gh0_57, gh1_49, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_21[k]
                  + f_12 * gh0_57[k]
                  - f_13 * gh1_49[k]
                  + pb_y[k] * gi_23[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_29, gh0_74, gh1_65, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_29[k]
                  + f_1 * gh0_74[k]
                  - f_2 * gh1_65[k]
                  + pb_z[k] * gi_26[k];
    }
}

auto
compute_prim_gk_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_32 = buffer.data(gh0 + 32);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_32 = buffer.data(gh1 + 32);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_1, dk1_1, fi_6, fi_7, fk_8, gh0_7, \
                         gh1_7, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_8[k];

        t_9[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_9, \
                         gh0_10, gh1_9, gh1_10, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_11, gh0_12, gh1_11, \
                         gh1_12, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, dk0_2, dk1_2, fi_13, fi_14, fi_15, \
                         fk_14, fk_15, fk_16, fk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = f_14 * fi_13[k]
                  + pa_x[k] * fk_15[k];

        t_18[k] = f_0 * fi_14[k]
                  + pa_x[k] * fk_16[k];

        t_19[k] = f_15 * fi_15[k]
                  + pa_x[k] * fk_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, fi_16, fi_17, fi_24, fi_25, \
                         fk_18, fk_19, fk_27, fk_28, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_16[k]
                  + pa_x[k] * fk_18[k];

        t_21[k] = f_16 * fi_17[k]
                  + pb_x[k] * gi_19[k];

        t_22[k] = pa_x[k] * fk_19[k];

        t_23[k] = f_14 * fi_24[k]
                  + pa_x[k] * fk_27[k];

        t_24[k] = f_0 * fi_25[k]
                  + pa_x[k] * fk_28[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, fi_26, fi_27, fi_32, fk_29, \
                         fk_30, fk_35, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_15 * fi_26[k]
                  + pa_x[k] * fk_29[k];

        t_26[k] = f_5 * fi_27[k]
                  + pa_x[k] * fk_30[k];

        t_27[k] = f_16 * fi_32[k]
                  + pb_x[k] * gi_24[k];

        t_28[k] = pa_x[k] * fk_35[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_z, pb_y, dk0_1, dk1_1, fi_17, fk_19, fk_20, \
                         gh0_21, gh1_21, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fi_17[k]
                  + f_1 * gh0_21[k]
                  - f_2 * gh1_21[k]
                  + pb_y[k] * gi_25[k];

        t_30[k] = pa_z[k] * fk_19[k];

        t_31[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_y, fi_19, fi_20, fi_21, gh0_24, gh0_25, gh0_26, \
                         gh1_24, gh1_25, gh1_26, gi_28, gi_29, gi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_5 * fi_19[k]
                  + f_6 * gh0_24[k]
                  - f_7 * gh1_24[k]
                  + pb_y[k] * gi_28[k];

        t_33[k] = f_5 * fi_20[k]
                  + f_8 * gh0_25[k]
                  - f_9 * gh1_25[k]
                  + pb_y[k] * gi_29[k];

        t_34[k] = f_5 * fi_21[k]
                  + f_10 * gh0_26[k]
                  - f_11 * gh1_26[k]
                  + pb_y[k] * gi_30[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_y, dk0_2, dk1_2, fi_22, fi_23, fk_26, \
                         gh0_27, gh1_27, gi_31, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * fi_22[k]
                  + f_12 * gh0_27[k]
                  - f_13 * gh1_27[k]
                  + pb_y[k] * gi_31[k];

        t_36[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_32[k];

        t_37[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_26[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, fi_28, fi_29, fi_30, fi_31, fk_31, \
                         fk_32, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_14 * fi_28[k]
                  + pa_y[k] * fk_31[k];

        t_39[k] = f_0 * fi_29[k]
                  + pa_y[k] * fk_32[k];

        t_40[k] = f_15 * fi_30[k]
                  + pa_y[k] * fk_33[k];

        t_41[k] = f_5 * fi_31[k]
                  + pa_y[k] * fk_34[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_y, pb_y, pb_z, fi_32, fk_35, gh0_32, gh1_32, \
                         gi_37, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_16 * fi_32[k]
                  + pb_y[k] * gi_37[k];

        t_43[k] = pa_y[k] * fk_35[k];

        t_44[k] = f_0 * fi_32[k]
                  + f_1 * gh0_32[k]
                  - f_2 * gh1_32[k]
                  + pb_z[k] * gi_38[k];
    }
}

auto
compute_prim_gk_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_32 = buffer.data(gh0 + 32);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_43 = buffer.data(gh1 + 43);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_54 = buffer.data(gh1 + 54);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_72 = buffer.data(gh1 + 72);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_90 = buffer.data(gi + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_16, \
                         gh1_17, gh1_18, gi_20, gi_21, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_16[k]
                 + pb_x[k] * gi_20[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_17[k]
                 + pb_x[k] * gi_21[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_18[k]
                 + pb_x[k] * gi_22[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_1, dk1_5, fi_6, fi_7, fk_8, gh0_7, \
                         gh1_19, gi_23, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_19[k]
                 + pb_x[k] * gi_23[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_24[k];

        t_9[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_9, \
                         gh0_10, gh1_22, gh1_23, gi_27, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_22[k]
                  + pb_x[k] * gi_27[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_23[k]
                  + pb_x[k] * gi_28[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_11, gh0_12, gh1_24, \
                         gh1_26, gi_29, gi_30, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_24[k]
                  + pb_x[k] * gi_29[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_26[k]
                  + pb_x[k] * gi_30[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_32[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, dk0_2, dk1_14, fi_13, fi_14, fi_15, \
                         fk_14, fk_15, fk_16, fk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = f_14 * fi_13[k]
                  + pa_x[k] * fk_15[k];

        t_18[k] = f_0 * fi_14[k]
                  + pa_x[k] * fk_16[k];

        t_19[k] = f_15 * fi_15[k]
                  + pa_x[k] * fk_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, fi_16, fi_17, fi_24, fi_25, \
                         fk_18, fk_19, fk_27, fk_28, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_16[k]
                  + pa_x[k] * fk_18[k];

        t_21[k] = f_16 * fi_17[k]
                  + pb_x[k] * gi_38[k];

        t_22[k] = pa_x[k] * fk_19[k];

        t_23[k] = f_14 * fi_24[k]
                  + pa_x[k] * fk_27[k];

        t_24[k] = f_0 * fi_25[k]
                  + pa_x[k] * fk_28[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, fi_26, fi_27, fi_32, fk_29, \
                         fk_30, fk_35, gi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_15 * fi_26[k]
                  + pa_x[k] * fk_29[k];

        t_26[k] = f_5 * fi_27[k]
                  + pa_x[k] * fk_30[k];

        t_27[k] = f_16 * fi_32[k]
                  + pb_x[k] * gi_44[k];

        t_28[k] = pa_x[k] * fk_35[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_z, pb_y, dk0_1, dk1_5, fi_17, fk_19, fk_20, \
                         gh0_21, gh1_43, gi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fi_17[k]
                  + f_1 * gh0_21[k]
                  - f_2 * gh1_43[k]
                  + pb_y[k] * gi_54[k];

        t_30[k] = pa_z[k] * fk_19[k];

        t_31[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_y, fi_19, fi_20, fi_21, gh0_24, gh0_25, gh0_26, \
                         gh1_52, gh1_53, gh1_54, gi_65, gi_66, gi_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_5 * fi_19[k]
                  + f_6 * gh0_24[k]
                  - f_7 * gh1_52[k]
                  + pb_y[k] * gi_65[k];

        t_33[k] = f_5 * fi_20[k]
                  + f_8 * gh0_25[k]
                  - f_9 * gh1_53[k]
                  + pb_y[k] * gi_66[k];

        t_34[k] = f_5 * fi_21[k]
                  + f_10 * gh0_26[k]
                  - f_11 * gh1_54[k]
                  + pb_y[k] * gi_67[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_y, dk0_2, dk1_14, fi_22, fi_23, fk_26, \
                         gh0_27, gh1_55, gi_68, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * fi_22[k]
                  + f_12 * gh0_27[k]
                  - f_13 * gh1_55[k]
                  + pb_y[k] * gi_68[k];

        t_36[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_69[k];

        t_37[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_26[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, fi_28, fi_29, fi_30, fi_31, fk_31, \
                         fk_32, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_14 * fi_28[k]
                  + pa_y[k] * fk_31[k];

        t_39[k] = f_0 * fi_29[k]
                  + pa_y[k] * fk_32[k];

        t_40[k] = f_15 * fi_30[k]
                  + pa_y[k] * fk_33[k];

        t_41[k] = f_5 * fi_31[k]
                  + pa_y[k] * fk_34[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_y, pb_y, pb_z, fi_32, fk_35, gh0_32, gh1_72, \
                         gi_75, gi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_16 * fi_32[k]
                  + pb_y[k] * gi_75[k];

        t_43[k] = pa_y[k] * fk_35[k];

        t_44[k] = f_0 * fi_32[k]
                  + f_1 * gh0_32[k]
                  - f_2 * gh1_72[k]
                  + pb_z[k] * gi_90[k];
    }
}

auto
compute_prim_gk_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);
    const auto f_13 = 1.0 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_1 = buffer.data(gh0 + 1);
    const auto *gh0_2 = buffer.data(gh0 + 2);
    const auto *gh0_3 = buffer.data(gh0 + 3);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_41 = buffer.data(gh0 + 41);
    const auto *gh0_42 = buffer.data(gh0 + 42);
    const auto *gh0_43 = buffer.data(gh0 + 43);
    const auto *gh0_44 = buffer.data(gh0 + 44);
    const auto *gh0_45 = buffer.data(gh0 + 45);
    const auto *gh0_46 = buffer.data(gh0 + 46);
    const auto *gh0_47 = buffer.data(gh0 + 47);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_54 = buffer.data(gh0 + 54);
    const auto *gh0_55 = buffer.data(gh0 + 55);
    const auto *gh0_60 = buffer.data(gh0 + 60);
    const auto *gh0_62 = buffer.data(gh0 + 62);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_64 = buffer.data(gh0 + 64);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_67 = buffer.data(gh0 + 67);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_71 = buffer.data(gh0 + 71);
    const auto *gh0_72 = buffer.data(gh0 + 72);

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
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_42 = buffer.data(gh1 + 42);
    const auto *gh1_43 = buffer.data(gh1 + 43);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_50 = buffer.data(gh1 + 50);
    const auto *gh1_51 = buffer.data(gh1 + 51);
    const auto *gh1_52 = buffer.data(gh1 + 52);
    const auto *gh1_53 = buffer.data(gh1 + 53);
    const auto *gh1_58 = buffer.data(gh1 + 58);
    const auto *gh1_59 = buffer.data(gh1 + 59);
    const auto *gh1_60 = buffer.data(gh1 + 60);
    const auto *gh1_61 = buffer.data(gh1 + 61);
    const auto *gh1_62 = buffer.data(gh1 + 62);
    const auto *gh1_63 = buffer.data(gh1 + 63);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_67 = buffer.data(gh1 + 67);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);

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
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_33 = buffer.data(gi + 33);
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
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fi_0, gh0_0, gh0_1, gh1_0, \
                         gh1_1, gi_0, gi_1, gi_2, gi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_2[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];

        t_3[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, gh0_2, gh0_3, gh0_4, gh1_2, gh1_3, gh1_4, \
                         gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_4[k];

        t_5[k] = f_7 * gh0_3[k]
                 - f_8 * gh1_3[k]
                 + pb_y[k] * gi_5[k];

        t_6[k] = f_7 * gh0_4[k]
                 - f_8 * gh1_4[k]
                 + pb_z[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, gh0_5, gh0_6, gh0_7, gh1_5, gh1_6, gh1_7, \
                         gi_7, gi_8, gi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_9 * gh0_5[k]
                 - f_10 * gh1_5[k]
                 + pb_y[k] * gi_7[k];

        t_8[k] = f_9 * gh0_6[k]
                 - f_10 * gh1_6[k]
                 + pb_z[k] * gi_8[k];

        t_9[k] = f_1 * gh0_7[k]
                 - f_2 * gh1_7[k]
                 + pb_y[k] * gi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, gh0_9, gh0_10, gh0_11, gh1_8, gh1_9, gh1_10, \
                         gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_8[k]
                  + pb_y[k] * gi_10[k];

        t_11[k] = f_7 * gh0_10[k]
                  - f_8 * gh1_9[k]
                  + pb_y[k] * gi_11[k];

        t_12[k] = f_5 * gh0_11[k]
                  - f_6 * gh1_10[k]
                  + pb_y[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, pb_z, fk_0, gh0_12, gh1_11, \
                         gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * gh0_12[k]
                  - f_4 * gh1_11[k]
                  + pb_y[k] * gi_13[k];

        t_14[k] = f_1 * gh0_12[k]
                  - f_2 * gh1_11[k]
                  + pb_z[k] * gi_14[k];

        t_15[k] = pa_y[k] * fk_0[k];

        t_16[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pb_x, dk0_0, dk1_0, fi_5, fi_6, fk_1, gh0_16, \
                         gh0_17, gh1_15, gh1_16, gi_18, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_11 * dk0_0[k]
                  - f_12 * dk1_0[k]
                  + pa_y[k] * fk_1[k];

        t_18[k] = f_13 * fi_5[k]
                  + f_9 * gh0_16[k]
                  - f_10 * gh1_15[k]
                  + pb_x[k] * gi_18[k];

        t_19[k] = f_13 * fi_6[k]
                  + f_7 * gh0_17[k]
                  - f_8 * gh1_16[k]
                  + pb_x[k] * gi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, fi_7, fi_8, fi_9, gh0_18, gh0_19, gh1_17, \
                         gh1_18, gi_20, gi_21, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_13 * fi_7[k]
                  + f_5 * gh0_18[k]
                  - f_6 * gh1_17[k]
                  + pb_x[k] * gi_20[k];

        t_21[k] = f_13 * fi_8[k]
                  + f_3 * gh0_19[k]
                  - f_4 * gh1_18[k]
                  + pb_x[k] * gi_21[k];

        t_22[k] = f_13 * fi_9[k]
                  + pb_x[k] * gi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_z, pb_x, dk0_0, dk0_5, dk1_0, dk1_5, \
                         fi_10, fk_2, fk_8, gh0_22, gh1_21, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * dk0_5[k]
                  - f_12 * dk1_5[k]
                  + pa_x[k] * fk_8[k];

        t_24[k] = f_11 * dk0_0[k]
                  - f_12 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_25[k] = f_13 * fi_10[k]
                  + f_9 * gh0_22[k]
                  - f_10 * gh1_21[k]
                  + pb_x[k] * gi_24[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_x, fi_11, fi_12, fi_13, gh0_23, gh0_24, gh0_26, \
                         gh1_22, gh1_23, gh1_25, gi_25, gi_26, gi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_13 * fi_11[k]
                  + f_7 * gh0_23[k]
                  - f_8 * gh1_22[k]
                  + pb_x[k] * gi_25[k];

        t_27[k] = f_13 * fi_12[k]
                  + f_5 * gh0_24[k]
                  - f_6 * gh1_23[k]
                  + pb_x[k] * gi_26[k];

        t_28[k] = f_13 * fi_13[k]
                  + f_3 * gh0_26[k]
                  - f_4 * gh1_25[k]
                  + pb_x[k] * gi_27[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pb_x, dk0_14, dk1_14, fi_14, fi_16, \
                         fi_17, fk_14, fk_15, fk_16, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_13 * fi_14[k]
                  + pb_x[k] * gi_28[k];

        t_30[k] = f_11 * dk0_14[k]
                  - f_12 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_31[k] = f_14 * fi_16[k]
                  + pa_x[k] * fk_15[k];

        t_32[k] = f_0 * fi_17[k]
                  + pa_x[k] * fk_16[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pb_x, fi_18, fi_19, fi_20, fi_29, \
                         fk_17, fk_18, fk_19, fk_27, gi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_15 * fi_18[k]
                  + pa_x[k] * fk_17[k];

        t_34[k] = f_13 * fi_19[k]
                  + pa_x[k] * fk_18[k];

        t_35[k] = f_16 * fi_20[k]
                  + pb_x[k] * gi_33[k];

        t_36[k] = pa_x[k] * fk_19[k];

        t_37[k] = f_14 * fi_29[k]
                  + pa_x[k] * fk_27[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_x, pb_x, fi_30, fi_31, fi_32, fi_38, \
                         fk_28, fk_29, fk_30, fk_35, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * fi_30[k]
                  + pa_x[k] * fk_28[k];

        t_39[k] = f_15 * fi_31[k]
                  + pa_x[k] * fk_29[k];

        t_40[k] = f_13 * fi_32[k]
                  + pa_x[k] * fk_30[k];

        t_41[k] = f_16 * fi_38[k]
                  + pb_x[k] * gi_38[k];

        t_42[k] = pa_x[k] * fk_35[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, gh0_35, gh0_37, gh0_38, gh1_34, gh1_35, \
                         gh1_36, gi_39, gi_40, gi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * gh0_35[k]
                  - f_2 * gh1_34[k]
                  + pb_x[k] * gi_39[k];

        t_44[k] = f_9 * gh0_37[k]
                  - f_10 * gh1_35[k]
                  + pb_x[k] * gi_40[k];

        t_45[k] = f_9 * gh0_38[k]
                  - f_10 * gh1_36[k]
                  + pb_x[k] * gi_41[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, gh0_39, gh0_40, gh0_41, gh1_37, gh1_38, \
                         gh1_39, gi_42, gi_43, gi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_7 * gh0_39[k]
                  - f_8 * gh1_37[k]
                  + pb_x[k] * gi_42[k];

        t_47[k] = f_7 * gh0_40[k]
                  - f_8 * gh1_38[k]
                  + pb_x[k] * gi_43[k];

        t_48[k] = f_5 * gh0_41[k]
                  - f_6 * gh1_39[k]
                  + pb_x[k] * gi_44[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, gh0_42, gh0_43, gh0_47, gh1_40, gh1_41, \
                         gh1_45, gi_45, gi_46, gi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * gh0_42[k]
                  - f_6 * gh1_40[k]
                  + pb_x[k] * gi_45[k];

        t_50[k] = f_3 * gh0_43[k]
                  - f_4 * gh1_41[k]
                  + pb_x[k] * gi_46[k];

        t_51[k] = f_3 * gh0_47[k]
                  - f_4 * gh1_45[k]
                  + pb_x[k] * gi_47[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_y, pb_z, fi_20, gh0_43, gh0_44, gh1_41, gh1_42, \
                         gi_48, gi_49, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * fi_20[k]
                  + f_1 * gh0_43[k]
                  - f_2 * gh1_41[k]
                  + pb_y[k] * gi_48[k];

        t_53[k] = f_3 * gh0_43[k]
                  - f_4 * gh1_41[k]
                  + pb_z[k] * gi_49[k];

        t_54[k] = f_5 * gh0_44[k]
                  - f_6 * gh1_42[k]
                  + pb_z[k] * gi_50[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_z, fk_19, gh0_45, gh0_46, gh0_47, \
                         gh1_43, gh1_44, gh1_45, gi_51, gi_52, gi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_7 * gh0_45[k]
                  - f_8 * gh1_43[k]
                  + pb_z[k] * gi_51[k];

        t_56[k] = f_9 * gh0_46[k]
                  - f_10 * gh1_44[k]
                  + pb_z[k] * gi_52[k];

        t_57[k] = f_1 * gh0_47[k]
                  - f_2 * gh1_45[k]
                  + pb_z[k] * gi_53[k];

        t_58[k] = pa_z[k] * fk_19[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_z, pb_y, dk0_5, dk1_5, fi_23, fi_24, fk_20, \
                         gh0_52, gh0_53, gh1_50, gh1_51, gi_56, gi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_11 * dk0_5[k]
                  - f_12 * dk1_5[k]
                  + pa_z[k] * fk_20[k];

        t_60[k] = f_13 * fi_23[k]
                  + f_9 * gh0_52[k]
                  - f_10 * gh1_50[k]
                  + pb_y[k] * gi_56[k];

        t_61[k] = f_13 * fi_24[k]
                  + f_7 * gh0_53[k]
                  - f_8 * gh1_51[k]
                  + pb_y[k] * gi_57[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pb_y, fi_25, fi_26, fi_27, gh0_54, gh0_55, gh1_52, \
                         gh1_53, gi_58, gi_59, gi_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_13 * fi_25[k]
                  + f_5 * gh0_54[k]
                  - f_6 * gh1_52[k]
                  + pb_y[k] * gi_58[k];

        t_63[k] = f_13 * fi_26[k]
                  + f_3 * gh0_55[k]
                  - f_4 * gh1_53[k]
                  + pb_y[k] * gi_59[k];

        t_64[k] = f_13 * fi_27[k]
                  + pb_y[k] * gi_60[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_y, dk0_14, dk1_14, fi_34, fi_35, fi_36, \
                         fk_26, fk_31, fk_32, fk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * dk0_14[k]
                  - f_12 * dk1_14[k]
                  + pa_y[k] * fk_26[k];

        t_66[k] = f_14 * fi_34[k]
                  + pa_y[k] * fk_31[k];

        t_67[k] = f_0 * fi_35[k]
                  + pa_y[k] * fk_32[k];

        t_68[k] = f_15 * fi_36[k]
                  + pa_y[k] * fk_33[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_x, pb_y, fi_37, fi_38, fk_34, fk_35, \
                         gh0_60, gh1_58, gi_65, gi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_13 * fi_37[k]
                  + pa_y[k] * fk_34[k];

        t_70[k] = f_16 * fi_38[k]
                  + pb_y[k] * gi_65[k];

        t_71[k] = pa_y[k] * fk_35[k];

        t_72[k] = f_1 * gh0_60[k]
                  - f_2 * gh1_58[k]
                  + pb_x[k] * gi_66[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, gh0_62, gh0_63, gh0_64, gh1_59, gh1_60, \
                         gh1_61, gi_67, gi_68, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_9 * gh0_62[k]
                  - f_10 * gh1_59[k]
                  + pb_x[k] * gi_67[k];

        t_74[k] = f_9 * gh0_63[k]
                  - f_10 * gh1_60[k]
                  + pb_x[k] * gi_68[k];

        t_75[k] = f_7 * gh0_64[k]
                  - f_8 * gh1_61[k]
                  + pb_x[k] * gi_69[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_x, gh0_65, gh0_66, gh0_67, gh1_62, gh1_63, \
                         gh1_64, gi_70, gi_71, gi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_7 * gh0_65[k]
                  - f_8 * gh1_62[k]
                  + pb_x[k] * gi_70[k];

        t_77[k] = f_5 * gh0_66[k]
                  - f_6 * gh1_63[k]
                  + pb_x[k] * gi_71[k];

        t_78[k] = f_5 * gh0_67[k]
                  - f_6 * gh1_64[k]
                  + pb_x[k] * gi_72[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pb_y, gh0_68, gh0_69, gh0_72, gh1_65, \
                         gh1_66, gh1_69, gi_73, gi_74, gi_75, gi_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * gh0_68[k]
                  - f_4 * gh1_65[k]
                  + pb_x[k] * gi_73[k];

        t_80[k] = f_3 * gh0_72[k]
                  - f_4 * gh1_69[k]
                  + pb_x[k] * gi_74[k];

        t_81[k] = f_1 * gh0_68[k]
                  - f_2 * gh1_65[k]
                  + pb_y[k] * gi_75[k];

        t_82[k] = f_9 * gh0_69[k]
                  - f_10 * gh1_66[k]
                  + pb_y[k] * gi_76[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, gh0_70, gh0_71, gh0_72, gh1_67, gh1_68, \
                         gh1_69, gi_77, gi_78, gi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_7 * gh0_70[k]
                  - f_8 * gh1_67[k]
                  + pb_y[k] * gi_77[k];

        t_84[k] = f_5 * gh0_71[k]
                  - f_6 * gh1_68[k]
                  + pb_y[k] * gi_78[k];

        t_85[k] = f_3 * gh0_72[k]
                  - f_4 * gh1_69[k]
                  + pb_y[k] * gi_79[k];
    }

#pragma omp simd aligned(t_86, pb_z, fi_38, gh0_72, gh1_69, gi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_0 * fi_38[k]
                  + f_1 * gh0_72[k]
                  - f_2 * gh1_69[k]
                  + pb_z[k] * gi_80[k];
    }
}

auto
compute_prim_gk_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_32 = buffer.data(fi + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_41 = buffer.data(gh0 + 41);
    const auto *gh0_50 = buffer.data(gh0 + 50);
    const auto *gh0_51 = buffer.data(gh0 + 51);
    const auto *gh0_52 = buffer.data(gh0 + 52);
    const auto *gh0_53 = buffer.data(gh0 + 53);
    const auto *gh0_69 = buffer.data(gh0 + 69);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);
    const auto *gh1_48 = buffer.data(gh1 + 48);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_65 = buffer.data(gh1 + 65);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_26 = buffer.data(gi + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_15, gh0_16, gh0_17, \
                         gh1_15, gh1_16, gh1_17, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_15[k]
                 - f_7 * gh1_15[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_16[k]
                 - f_9 * gh1_16[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_17[k]
                 - f_11 * gh1_17[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_5, dk1_5, fi_6, fi_7, fk_3, gh0_18, \
                         gh1_18, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_18[k]
                 - f_13 * gh1_18[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_8[k];

        t_9[k] = f_3 * dk0_5[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_21, \
                         gh0_22, gh1_20, gh1_21, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_21[k]
                  - f_7 * gh1_20[k]
                  + pb_x[k] * gi_10[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_22[k]
                  - f_9 * gh1_21[k]
                  + pb_x[k] * gi_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_23, gh0_25, gh1_22, \
                         gh1_23, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_23[k]
                  - f_11 * gh1_22[k]
                  + pb_x[k] * gi_12[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_25[k]
                  - f_13 * gh1_23[k]
                  + pb_x[k] * gi_13[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dk0_14, dk1_14, fi_17, fk_4, \
                         fk_5, fk_8, gh0_41, gh1_39, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_4[k];

        t_17[k] = pa_x[k] * fk_5[k];

        t_18[k] = pa_x[k] * fk_8[k];

        t_19[k] = f_0 * fi_17[k]
                  + f_1 * gh0_41[k]
                  - f_2 * gh1_39[k]
                  + pb_y[k] * gi_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dk0_5, dk1_5, fi_19, fk_5, fk_6, \
                         gh0_50, gh1_46, gi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fk_5[k];

        t_21[k] = f_3 * dk0_5[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_6[k];

        t_22[k] = f_5 * fi_19[k]
                  + f_6 * gh0_50[k]
                  - f_7 * gh1_46[k]
                  + pb_y[k] * gi_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fi_20, fi_21, fi_22, gh0_51, gh0_52, gh0_53, \
                         gh1_47, gh1_48, gh1_49, gi_21, gi_22, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_20[k]
                  + f_8 * gh0_51[k]
                  - f_9 * gh1_47[k]
                  + pb_y[k] * gi_21[k];

        t_24[k] = f_5 * fi_21[k]
                  + f_10 * gh0_52[k]
                  - f_11 * gh1_48[k]
                  + pb_y[k] * gi_22[k];

        t_25[k] = f_5 * fi_22[k]
                  + f_12 * gh0_53[k]
                  - f_13 * gh1_49[k]
                  + pb_y[k] * gi_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dk0_14, dk1_14, fi_23, fk_7, fk_8, \
                         gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_24[k];

        t_27[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_7[k];

        t_28[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fi_32, gh0_69, gh1_65, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fi_32[k]
                  + f_1 * gh0_69[k]
                  - f_2 * gh1_65[k]
                  + pb_z[k] * gi_26[k];
    }
}

auto
compute_prim_gk_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_20 = buffer.data(fi + 20);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_41 = buffer.data(gh1 + 41);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_47 = buffer.data(gi + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_7, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_9, gh1_10, gh1_11, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_12, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_13, fk_5, fk_6, \
                         gh0_15, gh1_24, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_13[k]
                  + f_1 * gh0_15[k]
                  - f_2 * gh1_24[k]
                  + pb_y[k] * gi_28[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_15, fi_16, fi_17, gh0_18, gh0_19, gh0_20, \
                         gh1_27, gh1_28, gh1_29, gi_31, gi_32, gi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_15[k]
                  + f_6 * gh0_18[k]
                  - f_7 * gh1_27[k]
                  + pb_y[k] * gi_31[k];

        t_21[k] = f_5 * fi_16[k]
                  + f_8 * gh0_19[k]
                  - f_9 * gh1_28[k]
                  + pb_y[k] * gi_32[k];

        t_22[k] = f_5 * fi_17[k]
                  + f_10 * gh0_20[k]
                  - f_11 * gh1_29[k]
                  + pb_y[k] * gi_33[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_18, fk_7, fk_8, \
                         gh0_21, gh1_30, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_18[k]
                  + f_12 * gh0_21[k]
                  - f_13 * gh1_30[k]
                  + pb_y[k] * gi_34[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_20, gh0_23, gh1_41, gi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_20[k]
                  + f_1 * gh0_23[k]
                  - f_2 * gh1_41[k]
                  + pb_z[k] * gi_47[k];
    }
}

auto
compute_prim_gk_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_29 = buffer.data(fi + 29);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_41 = buffer.data(gh0 + 41);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_41 = buffer.data(gh1 + 41);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_47 = buffer.data(gi + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_7, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_9, gh1_10, gh1_11, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_12, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_16, fk_5, fk_6, \
                         gh0_24, gh1_24, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_16[k]
                  + f_1 * gh0_24[k]
                  - f_2 * gh1_24[k]
                  + pb_y[k] * gi_28[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_18, fi_19, fi_20, gh0_27, gh0_28, gh0_29, \
                         gh1_27, gh1_28, gh1_29, gi_31, gi_32, gi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_18[k]
                  + f_6 * gh0_27[k]
                  - f_7 * gh1_27[k]
                  + pb_y[k] * gi_31[k];

        t_21[k] = f_5 * fi_19[k]
                  + f_8 * gh0_28[k]
                  - f_9 * gh1_28[k]
                  + pb_y[k] * gi_32[k];

        t_22[k] = f_5 * fi_20[k]
                  + f_10 * gh0_29[k]
                  - f_11 * gh1_29[k]
                  + pb_y[k] * gi_33[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_21, fk_7, fk_8, \
                         gh0_30, gh1_30, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_21[k]
                  + f_12 * gh0_30[k]
                  - f_13 * gh1_30[k]
                  + pb_y[k] * gi_34[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_29, gh0_41, gh1_41, gi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_29[k]
                  + f_1 * gh0_41[k]
                  - f_2 * gh1_41[k]
                  + pb_z[k] * gi_47[k];
    }
}

auto
compute_prim_gk_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_29 = buffer.data(fi + 29);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_41 = buffer.data(gh0 + 41);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_14 = buffer.data(gh1 + 14);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_47 = buffer.data(gh1 + 47);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_47 = buffer.data(gi + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_6, \
                         gh1_7, gh1_8, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_6[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_7[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_8[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_7, gh1_9, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_9[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_9, gh0_10, gh0_11, \
                         gh1_11, gh1_12, gh1_13, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_11[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_12[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_13[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_12, gh1_14, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_14[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_16, fk_5, fk_6, \
                         gh0_24, gh1_27, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_16[k]
                  + f_1 * gh0_24[k]
                  - f_2 * gh1_27[k]
                  + pb_y[k] * gi_28[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_18, fi_19, fi_20, gh0_27, gh0_28, gh0_29, \
                         gh1_31, gh1_32, gh1_33, gi_31, gi_32, gi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_18[k]
                  + f_6 * gh0_27[k]
                  - f_7 * gh1_31[k]
                  + pb_y[k] * gi_31[k];

        t_21[k] = f_5 * fi_19[k]
                  + f_8 * gh0_28[k]
                  - f_9 * gh1_32[k]
                  + pb_y[k] * gi_32[k];

        t_22[k] = f_5 * fi_20[k]
                  + f_10 * gh0_29[k]
                  - f_11 * gh1_33[k]
                  + pb_y[k] * gi_33[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_21, fk_7, fk_8, \
                         gh0_30, gh1_34, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_21[k]
                  + f_12 * gh0_30[k]
                  - f_13 * gh1_34[k]
                  + pb_y[k] * gi_34[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_29, gh0_41, gh1_47, gi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_29[k]
                  + f_1 * gh0_41[k]
                  - f_2 * gh1_47[k]
                  + pb_z[k] * gi_47[k];
    }
}

auto
compute_prim_gk_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_1 = buffer.data(dk1 + 1);
    const auto *dk1_2 = buffer.data(dk1 + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_29 = buffer.data(fi + 29);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_14 = buffer.data(gh0 + 14);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_47 = buffer.data(gh0 + 47);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_41 = buffer.data(gh1 + 41);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_26 = buffer.data(gi + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_6, gh0_7, gh0_8, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_6[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_7[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_8[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dk0_0, dk0_1, dk1_0, dk1_1, fi_6, \
                         fk_2, fk_3, gh0_9, gh1_7, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_9[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_1[k]
                 + pa_x[k] * fk_3[k];

        t_9[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, fi_8, fi_9, fi_10, gh0_11, gh0_12, gh0_13, \
                         gh1_9, gh1_10, gh1_11, gi_10, gi_11, gi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fi_8[k]
                  + f_6 * gh0_11[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_11[k] = f_5 * fi_9[k]
                  + f_8 * gh0_12[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];

        t_12[k] = f_5 * fi_10[k]
                  + f_10 * gh0_13[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dk0_2, dk1_2, fi_11, fk_4, fk_5, \
                         fk_8, gh0_14, gh1_12, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_11[k]
                  + f_12 * gh0_14[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_14[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_x[k] * fk_4[k];

        t_15[k] = pa_x[k] * fk_5[k];

        t_16[k] = pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, dk0_1, dk1_1, fi_16, fk_5, fk_6, \
                         gh0_27, gh1_24, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fi_16[k]
                  + f_1 * gh0_27[k]
                  - f_2 * gh1_24[k]
                  + pb_y[k] * gi_17[k];

        t_18[k] = pa_z[k] * fk_5[k];

        t_19[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_1[k]
                  + pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, fi_18, fi_19, fi_20, gh0_31, gh0_32, gh0_33, \
                         gh1_27, gh1_28, gh1_29, gi_20, gi_21, gi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_18[k]
                  + f_6 * gh0_31[k]
                  - f_7 * gh1_27[k]
                  + pb_y[k] * gi_20[k];

        t_21[k] = f_5 * fi_19[k]
                  + f_8 * gh0_32[k]
                  - f_9 * gh1_28[k]
                  + pb_y[k] * gi_21[k];

        t_22[k] = f_5 * fi_20[k]
                  + f_10 * gh0_33[k]
                  - f_11 * gh1_29[k]
                  + pb_y[k] * gi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_y, dk0_2, dk1_2, fi_21, fk_7, fk_8, \
                         gh0_34, gh1_30, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_21[k]
                  + f_12 * gh0_34[k]
                  - f_13 * gh1_30[k]
                  + pb_y[k] * gi_23[k];

        t_24[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_2[k]
                  + pa_y[k] * fk_7[k];

        t_25[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_26, pb_z, fi_29, gh0_47, gh1_41, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * fi_29[k]
                  + f_1 * gh0_47[k]
                  - f_2 * gh1_41[k]
                  + pb_z[k] * gi_26[k];
    }
}

auto
compute_prim_gk_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_1 = buffer.data(dk0 + 1);
    const auto *dk0_2 = buffer.data(dk0 + 2);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_41 = buffer.data(gh1 + 41);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_50 = buffer.data(gi + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_1, dk1_5, fi_6, fi_7, fk_8, gh0_7, \
                         gh1_7, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_8[k];

        t_9[k] = f_3 * dk0_1[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_9, \
                         gh0_10, gh1_9, gh1_10, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_11, gh0_12, gh1_11, \
                         gh1_12, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dk0_2, dk1_14, fi_13, fk_14, \
                         fk_15, fk_23, gh0_15, gh1_24, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = pa_x[k] * fk_15[k];

        t_18[k] = pa_x[k] * fk_23[k];

        t_19[k] = f_0 * fi_13[k]
                  + f_1 * gh0_15[k]
                  - f_2 * gh1_24[k]
                  + pb_y[k] * gi_29[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dk0_1, dk1_5, fi_15, fk_15, fk_16, \
                         gh0_18, gh1_27, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fk_15[k];

        t_21[k] = f_3 * dk0_1[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_16[k];

        t_22[k] = f_5 * fi_15[k]
                  + f_6 * gh0_18[k]
                  - f_7 * gh1_27[k]
                  + pb_y[k] * gi_32[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fi_16, fi_17, fi_18, gh0_19, gh0_20, gh0_21, \
                         gh1_28, gh1_29, gh1_30, gi_33, gi_34, gi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_16[k]
                  + f_8 * gh0_19[k]
                  - f_9 * gh1_28[k]
                  + pb_y[k] * gi_33[k];

        t_24[k] = f_5 * fi_17[k]
                  + f_10 * gh0_20[k]
                  - f_11 * gh1_29[k]
                  + pb_y[k] * gi_34[k];

        t_25[k] = f_5 * fi_18[k]
                  + f_12 * gh0_21[k]
                  - f_13 * gh1_30[k]
                  + pb_y[k] * gi_35[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dk0_2, dk1_14, fi_19, fk_22, fk_23, \
                         gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fi_19[k]
                  + pb_y[k] * gi_36[k];

        t_27[k] = f_3 * dk0_2[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_22[k];

        t_28[k] = pa_y[k] * fk_23[k];
    }

#pragma omp simd aligned(t_29, pb_z, fi_20, gh0_23, gh1_41, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fi_20[k]
                  + f_1 * gh0_23[k]
                  - f_2 * gh1_41[k]
                  + pb_z[k] * gi_50[k];
    }
}

auto
compute_prim_gk_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_41 = buffer.data(gh0 + 41);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_21 = buffer.data(gh1 + 21);
    const auto *gh1_22 = buffer.data(gh1 + 22);
    const auto *gh1_23 = buffer.data(gh1 + 23);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_35 = buffer.data(gh1 + 35);
    const auto *gh1_36 = buffer.data(gh1 + 36);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_38 = buffer.data(gh1 + 38);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_5, dk1_5, fi_6, fi_7, fk_8, gh0_7, \
                         gh1_7, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_8[k];

        t_9[k] = f_3 * dk0_5[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_9, \
                         gh0_10, gh1_9, gh1_10, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_11, gh0_12, gh1_11, \
                         gh1_12, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, dk0_14, dk1_14, fi_13, fi_14, fi_15, \
                         fk_14, fk_15, fk_16, fk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = f_14 * fi_13[k]
                  + pa_x[k] * fk_15[k];

        t_18[k] = f_0 * fi_14[k]
                  + pa_x[k] * fk_16[k];

        t_19[k] = f_15 * fi_15[k]
                  + pa_x[k] * fk_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, fi_16, fi_17, fi_24, fi_25, \
                         fk_18, fk_19, fk_27, fk_28, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_16[k]
                  + pa_x[k] * fk_18[k];

        t_21[k] = f_16 * fi_17[k]
                  + pb_x[k] * gi_19[k];

        t_22[k] = pa_x[k] * fk_19[k];

        t_23[k] = f_14 * fi_24[k]
                  + pa_x[k] * fk_27[k];

        t_24[k] = f_0 * fi_25[k]
                  + pa_x[k] * fk_28[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, fi_26, fi_27, fi_32, fk_29, \
                         fk_30, fk_35, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_15 * fi_26[k]
                  + pa_x[k] * fk_29[k];

        t_26[k] = f_5 * fi_27[k]
                  + pa_x[k] * fk_30[k];

        t_27[k] = f_16 * fi_32[k]
                  + pb_x[k] * gi_24[k];

        t_28[k] = pa_x[k] * fk_35[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, gh0_21, gh0_22, gh0_23, gh1_21, gh1_22, \
                         gh1_23, gi_25, gi_26, gi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_6 * gh0_21[k]
                  - f_7 * gh1_21[k]
                  + pb_x[k] * gi_25[k];

        t_30[k] = f_8 * gh0_22[k]
                  - f_9 * gh1_22[k]
                  + pb_x[k] * gi_26[k];

        t_31[k] = f_10 * gh0_23[k]
                  - f_11 * gh1_23[k]
                  + pb_x[k] * gi_27[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, dk0_5, dk1_5, fi_17, fk_19, \
                         fk_20, gh0_24, gh1_24, gi_28, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * gh0_24[k]
                  - f_13 * gh1_24[k]
                  + pb_x[k] * gi_28[k];

        t_33[k] = f_0 * fi_17[k]
                  + f_1 * gh0_24[k]
                  - f_2 * gh1_24[k]
                  + pb_y[k] * gi_29[k];

        t_34[k] = pa_z[k] * fk_19[k];

        t_35[k] = f_3 * dk0_5[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_20[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, fi_19, fi_20, fi_21, gh0_27, gh0_28, gh0_29, \
                         gh1_27, gh1_28, gh1_29, gi_32, gi_33, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * fi_19[k]
                  + f_6 * gh0_27[k]
                  - f_7 * gh1_27[k]
                  + pb_y[k] * gi_32[k];

        t_37[k] = f_5 * fi_20[k]
                  + f_8 * gh0_28[k]
                  - f_9 * gh1_28[k]
                  + pb_y[k] * gi_33[k];

        t_38[k] = f_5 * fi_21[k]
                  + f_10 * gh0_29[k]
                  - f_11 * gh1_29[k]
                  + pb_y[k] * gi_34[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, dk0_14, dk1_14, fi_22, fi_23, fk_26, \
                         gh0_30, gh1_30, gi_35, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * fi_22[k]
                  + f_12 * gh0_30[k]
                  - f_13 * gh1_30[k]
                  + pb_y[k] * gi_35[k];

        t_40[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_36[k];

        t_41[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, fi_28, fi_29, fi_30, fi_31, fk_31, \
                         fk_32, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_14 * fi_28[k]
                  + pa_y[k] * fk_31[k];

        t_43[k] = f_0 * fi_29[k]
                  + pa_y[k] * fk_32[k];

        t_44[k] = f_15 * fi_30[k]
                  + pa_y[k] * fk_33[k];

        t_45[k] = f_5 * fi_31[k]
                  + pa_y[k] * fk_34[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pb_x, pb_y, fi_32, fk_35, gh0_35, \
                         gh0_36, gh1_35, gh1_36, gi_41, gi_42, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_16 * fi_32[k]
                  + pb_y[k] * gi_41[k];

        t_47[k] = pa_y[k] * fk_35[k];

        t_48[k] = f_6 * gh0_35[k]
                  - f_7 * gh1_35[k]
                  + pb_x[k] * gi_42[k];

        t_49[k] = f_8 * gh0_36[k]
                  - f_9 * gh1_36[k]
                  + pb_x[k] * gi_43[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, pb_y, gh0_37, gh0_38, gh0_41, gh1_37, gh1_38, \
                         gh1_41, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_10 * gh0_37[k]
                  - f_11 * gh1_37[k]
                  + pb_x[k] * gi_44[k];

        t_51[k] = f_12 * gh0_41[k]
                  - f_13 * gh1_41[k]
                  + pb_x[k] * gi_45[k];

        t_52[k] = f_6 * gh0_38[k]
                  - f_7 * gh1_38[k]
                  + pb_y[k] * gi_46[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, gh0_39, gh0_40, gh0_41, gh1_39, gh1_40, \
                         gh1_41, gi_47, gi_48, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * gh0_39[k]
                  - f_9 * gh1_39[k]
                  + pb_y[k] * gi_47[k];

        t_54[k] = f_10 * gh0_40[k]
                  - f_11 * gh1_40[k]
                  + pb_y[k] * gi_48[k];

        t_55[k] = f_12 * gh0_41[k]
                  - f_13 * gh1_41[k]
                  + pb_y[k] * gi_49[k];
    }

#pragma omp simd aligned(t_56, pb_z, fi_32, gh0_41, gh1_41, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * fi_32[k]
                  + f_1 * gh0_41[k]
                  - f_2 * gh1_41[k]
                  + pb_z[k] * gi_50[k];
    }
}

auto
compute_prim_gk_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_28 = buffer.data(gh0 + 28);
    const auto *gh0_29 = buffer.data(gh0 + 29);
    const auto *gh0_30 = buffer.data(gh0 + 30);
    const auto *gh0_35 = buffer.data(gh0 + 35);
    const auto *gh0_36 = buffer.data(gh0 + 36);
    const auto *gh0_37 = buffer.data(gh0 + 37);
    const auto *gh0_38 = buffer.data(gh0 + 38);
    const auto *gh0_39 = buffer.data(gh0 + 39);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_41 = buffer.data(gh0 + 41);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_14 = buffer.data(gh1 + 14);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_42 = buffer.data(gh1 + 42);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_4, gh0_5, gh0_6, gh1_6, \
                         gh1_7, gh1_8, gi_6, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_4[k]
                 - f_7 * gh1_6[k]
                 + pb_x[k] * gi_6[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_5[k]
                 - f_9 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_6[k]
                 - f_11 * gh1_8[k]
                 + pb_x[k] * gi_8[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_5, dk1_5, fi_6, fi_7, fk_8, gh0_7, \
                         gh1_9, gi_9, gi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_7[k]
                 - f_13 * gh1_9[k]
                 + pb_x[k] * gi_9[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_10[k];

        t_9[k] = f_3 * dk0_5[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_9, \
                         gh0_10, gh1_11, gh1_12, gi_12, gi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_9[k]
                  - f_7 * gh1_11[k]
                  + pb_x[k] * gi_12[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_10[k]
                  - f_9 * gh1_12[k]
                  + pb_x[k] * gi_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_11, gh0_12, gh1_13, \
                         gh1_14, gi_14, gi_15, gi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_11[k]
                  - f_11 * gh1_13[k]
                  + pb_x[k] * gi_14[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_12[k]
                  - f_13 * gh1_14[k]
                  + pb_x[k] * gi_15[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_16[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, dk0_14, dk1_14, fi_13, fi_14, fi_15, \
                         fk_14, fk_15, fk_16, fk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = f_14 * fi_13[k]
                  + pa_x[k] * fk_15[k];

        t_18[k] = f_0 * fi_14[k]
                  + pa_x[k] * fk_16[k];

        t_19[k] = f_15 * fi_15[k]
                  + pa_x[k] * fk_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, fi_16, fi_17, fi_24, fi_25, \
                         fk_18, fk_19, fk_27, fk_28, gi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_16[k]
                  + pa_x[k] * fk_18[k];

        t_21[k] = f_16 * fi_17[k]
                  + pb_x[k] * gi_21[k];

        t_22[k] = pa_x[k] * fk_19[k];

        t_23[k] = f_14 * fi_24[k]
                  + pa_x[k] * fk_27[k];

        t_24[k] = f_0 * fi_25[k]
                  + pa_x[k] * fk_28[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, fi_26, fi_27, fi_32, fk_29, \
                         fk_30, fk_35, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_15 * fi_26[k]
                  + pa_x[k] * fk_29[k];

        t_26[k] = f_5 * fi_27[k]
                  + pa_x[k] * fk_30[k];

        t_27[k] = f_16 * fi_32[k]
                  + pb_x[k] * gi_26[k];

        t_28[k] = pa_x[k] * fk_35[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, gh0_21, gh0_22, gh0_23, gh1_24, gh1_25, \
                         gh1_26, gi_28, gi_29, gi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_6 * gh0_21[k]
                  - f_7 * gh1_24[k]
                  + pb_x[k] * gi_28[k];

        t_30[k] = f_8 * gh0_22[k]
                  - f_9 * gh1_25[k]
                  + pb_x[k] * gi_29[k];

        t_31[k] = f_10 * gh0_23[k]
                  - f_11 * gh1_26[k]
                  + pb_x[k] * gi_30[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, dk0_5, dk1_5, fi_17, fk_19, \
                         fk_20, gh0_24, gh1_27, gi_31, gi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * gh0_24[k]
                  - f_13 * gh1_27[k]
                  + pb_x[k] * gi_31[k];

        t_33[k] = f_0 * fi_17[k]
                  + f_1 * gh0_24[k]
                  - f_2 * gh1_27[k]
                  + pb_y[k] * gi_32[k];

        t_34[k] = pa_z[k] * fk_19[k];

        t_35[k] = f_3 * dk0_5[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_20[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, fi_19, fi_20, fi_21, gh0_27, gh0_28, gh0_29, \
                         gh1_31, gh1_32, gh1_33, gi_36, gi_37, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * fi_19[k]
                  + f_6 * gh0_27[k]
                  - f_7 * gh1_31[k]
                  + pb_y[k] * gi_36[k];

        t_37[k] = f_5 * fi_20[k]
                  + f_8 * gh0_28[k]
                  - f_9 * gh1_32[k]
                  + pb_y[k] * gi_37[k];

        t_38[k] = f_5 * fi_21[k]
                  + f_10 * gh0_29[k]
                  - f_11 * gh1_33[k]
                  + pb_y[k] * gi_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, dk0_14, dk1_14, fi_22, fi_23, fk_26, \
                         gh0_30, gh1_34, gi_39, gi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * fi_22[k]
                  + f_12 * gh0_30[k]
                  - f_13 * gh1_34[k]
                  + pb_y[k] * gi_39[k];

        t_40[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_40[k];

        t_41[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, fi_28, fi_29, fi_30, fi_31, fk_31, \
                         fk_32, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_14 * fi_28[k]
                  + pa_y[k] * fk_31[k];

        t_43[k] = f_0 * fi_29[k]
                  + pa_y[k] * fk_32[k];

        t_44[k] = f_15 * fi_30[k]
                  + pa_y[k] * fk_33[k];

        t_45[k] = f_5 * fi_31[k]
                  + pa_y[k] * fk_34[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pb_x, pb_y, fi_32, fk_35, gh0_35, \
                         gh0_36, gh1_40, gh1_41, gi_45, gi_47, gi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_16 * fi_32[k]
                  + pb_y[k] * gi_45[k];

        t_47[k] = pa_y[k] * fk_35[k];

        t_48[k] = f_6 * gh0_35[k]
                  - f_7 * gh1_40[k]
                  + pb_x[k] * gi_47[k];

        t_49[k] = f_8 * gh0_36[k]
                  - f_9 * gh1_41[k]
                  + pb_x[k] * gi_48[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, pb_y, gh0_37, gh0_38, gh0_41, gh1_42, gh1_44, \
                         gh1_47, gi_49, gi_50, gi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_10 * gh0_37[k]
                  - f_11 * gh1_42[k]
                  + pb_x[k] * gi_49[k];

        t_51[k] = f_12 * gh0_41[k]
                  - f_13 * gh1_47[k]
                  + pb_x[k] * gi_50[k];

        t_52[k] = f_6 * gh0_38[k]
                  - f_7 * gh1_44[k]
                  + pb_y[k] * gi_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, gh0_39, gh0_40, gh0_41, gh1_45, gh1_46, \
                         gh1_47, gi_53, gi_54, gi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * gh0_39[k]
                  - f_9 * gh1_45[k]
                  + pb_y[k] * gi_53[k];

        t_54[k] = f_10 * gh0_40[k]
                  - f_11 * gh1_46[k]
                  + pb_y[k] * gi_54[k];

        t_55[k] = f_12 * gh0_41[k]
                  - f_13 * gh1_47[k]
                  + pb_y[k] * gi_55[k];
    }

#pragma omp simd aligned(t_56, pb_z, fi_32, gh0_41, gh1_47, gi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * fi_32[k]
                  + f_1 * gh0_41[k]
                  - f_2 * gh1_47[k]
                  + pb_z[k] * gi_56[k];
    }
}

auto
compute_prim_gk_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_14 = buffer.data(gh0 + 14);
    const auto *gh0_24 = buffer.data(gh0 + 24);
    const auto *gh0_25 = buffer.data(gh0 + 25);
    const auto *gh0_26 = buffer.data(gh0 + 26);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_40 = buffer.data(gh0 + 40);
    const auto *gh0_41 = buffer.data(gh0 + 41);
    const auto *gh0_42 = buffer.data(gh0 + 42);
    const auto *gh0_44 = buffer.data(gh0 + 44);
    const auto *gh0_45 = buffer.data(gh0 + 45);
    const auto *gh0_46 = buffer.data(gh0 + 46);
    const auto *gh0_47 = buffer.data(gh0 + 47);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_14 = buffer.data(gh1 + 14);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_26 = buffer.data(gh1 + 26);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_31 = buffer.data(gh1 + 31);
    const auto *gh1_32 = buffer.data(gh1 + 32);
    const auto *gh1_33 = buffer.data(gh1 + 33);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_40 = buffer.data(gh1 + 40);
    const auto *gh1_41 = buffer.data(gh1 + 41);
    const auto *gh1_42 = buffer.data(gh1 + 42);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_45 = buffer.data(gh1 + 45);
    const auto *gh1_46 = buffer.data(gh1 + 46);
    const auto *gh1_47 = buffer.data(gh1 + 47);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_6, gh0_7, gh0_8, gh1_6, \
                         gh1_7, gh1_8, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_6[k]
                 - f_7 * gh1_6[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_7[k]
                 - f_9 * gh1_7[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_8[k]
                 - f_11 * gh1_8[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_5, dk1_5, fi_6, fi_7, fk_8, gh0_9, \
                         gh1_9, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_9[k]
                 - f_13 * gh1_9[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_8[k];

        t_9[k] = f_3 * dk0_5[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_11, \
                         gh0_12, gh1_11, gh1_12, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_11[k]
                  - f_7 * gh1_11[k]
                  + pb_x[k] * gi_10[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_12[k]
                  - f_9 * gh1_12[k]
                  + pb_x[k] * gi_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_13, gh0_14, gh1_13, \
                         gh1_14, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_13[k]
                  - f_11 * gh1_13[k]
                  + pb_x[k] * gi_12[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_14[k]
                  - f_13 * gh1_14[k]
                  + pb_x[k] * gi_13[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, dk0_14, dk1_14, fi_13, fi_14, fi_15, \
                         fk_14, fk_15, fk_16, fk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_14[k];

        t_17[k] = f_14 * fi_13[k]
                  + pa_x[k] * fk_15[k];

        t_18[k] = f_0 * fi_14[k]
                  + pa_x[k] * fk_16[k];

        t_19[k] = f_15 * fi_15[k]
                  + pa_x[k] * fk_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, fi_16, fi_17, fi_24, fi_25, \
                         fk_18, fk_19, fk_27, fk_28, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fi_16[k]
                  + pa_x[k] * fk_18[k];

        t_21[k] = f_16 * fi_17[k]
                  + pb_x[k] * gi_19[k];

        t_22[k] = pa_x[k] * fk_19[k];

        t_23[k] = f_14 * fi_24[k]
                  + pa_x[k] * fk_27[k];

        t_24[k] = f_0 * fi_25[k]
                  + pa_x[k] * fk_28[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, fi_26, fi_27, fi_32, fk_29, \
                         fk_30, fk_35, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_15 * fi_26[k]
                  + pa_x[k] * fk_29[k];

        t_26[k] = f_5 * fi_27[k]
                  + pa_x[k] * fk_30[k];

        t_27[k] = f_16 * fi_32[k]
                  + pb_x[k] * gi_24[k];

        t_28[k] = pa_x[k] * fk_35[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, gh0_24, gh0_25, gh0_26, gh1_24, gh1_25, \
                         gh1_26, gi_25, gi_26, gi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_6 * gh0_24[k]
                  - f_7 * gh1_24[k]
                  + pb_x[k] * gi_25[k];

        t_30[k] = f_8 * gh0_25[k]
                  - f_9 * gh1_25[k]
                  + pb_x[k] * gi_26[k];

        t_31[k] = f_10 * gh0_26[k]
                  - f_11 * gh1_26[k]
                  + pb_x[k] * gi_27[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, dk0_5, dk1_5, fi_17, fk_19, \
                         fk_20, gh0_27, gh1_27, gi_28, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * gh0_27[k]
                  - f_13 * gh1_27[k]
                  + pb_x[k] * gi_28[k];

        t_33[k] = f_0 * fi_17[k]
                  + f_1 * gh0_27[k]
                  - f_2 * gh1_27[k]
                  + pb_y[k] * gi_29[k];

        t_34[k] = pa_z[k] * fk_19[k];

        t_35[k] = f_3 * dk0_5[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_20[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, fi_19, fi_20, fi_21, gh0_31, gh0_32, gh0_33, \
                         gh1_31, gh1_32, gh1_33, gi_32, gi_33, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * fi_19[k]
                  + f_6 * gh0_31[k]
                  - f_7 * gh1_31[k]
                  + pb_y[k] * gi_32[k];

        t_37[k] = f_5 * fi_20[k]
                  + f_8 * gh0_32[k]
                  - f_9 * gh1_32[k]
                  + pb_y[k] * gi_33[k];

        t_38[k] = f_5 * fi_21[k]
                  + f_10 * gh0_33[k]
                  - f_11 * gh1_33[k]
                  + pb_y[k] * gi_34[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, dk0_14, dk1_14, fi_22, fi_23, fk_26, \
                         gh0_34, gh1_34, gi_35, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * fi_22[k]
                  + f_12 * gh0_34[k]
                  - f_13 * gh1_34[k]
                  + pb_y[k] * gi_35[k];

        t_40[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_36[k];

        t_41[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, fi_28, fi_29, fi_30, fi_31, fk_31, \
                         fk_32, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_14 * fi_28[k]
                  + pa_y[k] * fk_31[k];

        t_43[k] = f_0 * fi_29[k]
                  + pa_y[k] * fk_32[k];

        t_44[k] = f_15 * fi_30[k]
                  + pa_y[k] * fk_33[k];

        t_45[k] = f_5 * fi_31[k]
                  + pa_y[k] * fk_34[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pb_x, pb_y, fi_32, fk_35, gh0_40, \
                         gh0_41, gh1_40, gh1_41, gi_41, gi_42, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_16 * fi_32[k]
                  + pb_y[k] * gi_41[k];

        t_47[k] = pa_y[k] * fk_35[k];

        t_48[k] = f_6 * gh0_40[k]
                  - f_7 * gh1_40[k]
                  + pb_x[k] * gi_42[k];

        t_49[k] = f_8 * gh0_41[k]
                  - f_9 * gh1_41[k]
                  + pb_x[k] * gi_43[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, pb_y, gh0_42, gh0_44, gh0_47, gh1_42, gh1_44, \
                         gh1_47, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_10 * gh0_42[k]
                  - f_11 * gh1_42[k]
                  + pb_x[k] * gi_44[k];

        t_51[k] = f_12 * gh0_47[k]
                  - f_13 * gh1_47[k]
                  + pb_x[k] * gi_45[k];

        t_52[k] = f_6 * gh0_44[k]
                  - f_7 * gh1_44[k]
                  + pb_y[k] * gi_46[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, gh0_45, gh0_46, gh0_47, gh1_45, gh1_46, \
                         gh1_47, gi_47, gi_48, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * gh0_45[k]
                  - f_9 * gh1_45[k]
                  + pb_y[k] * gi_47[k];

        t_54[k] = f_10 * gh0_46[k]
                  - f_11 * gh1_46[k]
                  + pb_y[k] * gi_48[k];

        t_55[k] = f_12 * gh0_47[k]
                  - f_13 * gh1_47[k]
                  + pb_y[k] * gi_49[k];
    }

#pragma omp simd aligned(t_56, pb_z, fi_32, gh0_47, gh1_47, gi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * fi_32[k]
                  + f_1 * gh0_47[k]
                  - f_2 * gh1_47[k]
                  + pb_z[k] * gi_50[k];
    }
}

auto
compute_prim_gk_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk0, const size_t dk1,
                                      const size_t fi, const size_t fk, const size_t gh0,
                                      const size_t gh1, const size_t gi, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / beta;
    const auto f_7 = 2.0 * alpha / (beta * p);
    const auto f_8 = 1.5 / beta;
    const auto f_9 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / beta;
    const auto f_11 = alpha / (beta * p);
    const auto f_12 = 0.5 / beta;
    const auto f_13 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_5 = buffer.data(dk0 + 5);
    const auto *dk0_14 = buffer.data(dk0 + 14);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_5 = buffer.data(dk1 + 5);
    const auto *dk1_14 = buffer.data(dk1 + 14);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_32 = buffer.data(fi + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_14 = buffer.data(gh0 + 14);
    const auto *gh0_27 = buffer.data(gh0 + 27);
    const auto *gh0_31 = buffer.data(gh0 + 31);
    const auto *gh0_32 = buffer.data(gh0 + 32);
    const auto *gh0_33 = buffer.data(gh0 + 33);
    const auto *gh0_34 = buffer.data(gh0 + 34);
    const auto *gh0_47 = buffer.data(gh0 + 47);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_28 = buffer.data(gh1 + 28);
    const auto *gh1_29 = buffer.data(gh1 + 29);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_41 = buffer.data(gh1 + 41);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_26 = buffer.data(gi + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dk0_0, dk1_0, fi_0, fk_0, fk_1, \
                         gh0_0, gh1_0, gi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pa_y[k] * fk_0[k];

        t_2[k] = pa_z[k] * fk_0[k];

        t_3[k] = f_3 * dk0_0[k]
                 - f_4 * dk1_0[k]
                 + pa_y[k] * fk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fi_3, fi_4, fi_5, gh0_6, gh0_7, gh0_8, gh1_4, \
                         gh1_5, gh1_6, gi_4, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi_3[k]
                 + f_6 * gh0_6[k]
                 - f_7 * gh1_4[k]
                 + pb_x[k] * gi_4[k];

        t_5[k] = f_5 * fi_4[k]
                 + f_8 * gh0_7[k]
                 - f_9 * gh1_5[k]
                 + pb_x[k] * gi_5[k];

        t_6[k] = f_5 * fi_5[k]
                 + f_10 * gh0_8[k]
                 - f_11 * gh1_6[k]
                 + pb_x[k] * gi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dk0_5, dk1_5, fi_6, fi_7, fk_3, gh0_9, \
                         gh1_7, gi_7, gi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fi_6[k]
                 + f_12 * gh0_9[k]
                 - f_13 * gh1_7[k]
                 + pb_x[k] * gi_7[k];

        t_8[k] = f_5 * fi_7[k]
                 + pb_x[k] * gi_8[k];

        t_9[k] = f_3 * dk0_5[k]
                 - f_4 * dk1_5[k]
                 + pa_x[k] * fk_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dk0_0, dk1_0, fi_8, fi_9, fk_2, gh0_11, \
                         gh0_12, gh1_9, gh1_10, gi_10, gi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dk0_0[k]
                  - f_4 * dk1_0[k]
                  + pa_z[k] * fk_2[k];

        t_11[k] = f_5 * fi_8[k]
                  + f_6 * gh0_11[k]
                  - f_7 * gh1_9[k]
                  + pb_x[k] * gi_10[k];

        t_12[k] = f_5 * fi_9[k]
                  + f_8 * gh0_12[k]
                  - f_9 * gh1_10[k]
                  + pb_x[k] * gi_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fi_10, fi_11, fi_12, gh0_13, gh0_14, gh1_11, \
                         gh1_12, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fi_10[k]
                  + f_10 * gh0_13[k]
                  - f_11 * gh1_11[k]
                  + pb_x[k] * gi_12[k];

        t_14[k] = f_5 * fi_11[k]
                  + f_12 * gh0_14[k]
                  - f_13 * gh1_12[k]
                  + pb_x[k] * gi_13[k];

        t_15[k] = f_5 * fi_12[k]
                  + pb_x[k] * gi_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dk0_14, dk1_14, fi_17, fk_4, \
                         fk_5, fk_8, gh0_27, gh1_24, gi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_x[k] * fk_4[k];

        t_17[k] = pa_x[k] * fk_5[k];

        t_18[k] = pa_x[k] * fk_8[k];

        t_19[k] = f_0 * fi_17[k]
                  + f_1 * gh0_27[k]
                  - f_2 * gh1_24[k]
                  + pb_y[k] * gi_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dk0_5, dk1_5, fi_19, fk_5, fk_6, \
                         gh0_31, gh1_27, gi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fk_5[k];

        t_21[k] = f_3 * dk0_5[k]
                  - f_4 * dk1_5[k]
                  + pa_z[k] * fk_6[k];

        t_22[k] = f_5 * fi_19[k]
                  + f_6 * gh0_31[k]
                  - f_7 * gh1_27[k]
                  + pb_y[k] * gi_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fi_20, fi_21, fi_22, gh0_32, gh0_33, gh0_34, \
                         gh1_28, gh1_29, gh1_30, gi_21, gi_22, gi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fi_20[k]
                  + f_8 * gh0_32[k]
                  - f_9 * gh1_28[k]
                  + pb_y[k] * gi_21[k];

        t_24[k] = f_5 * fi_21[k]
                  + f_10 * gh0_33[k]
                  - f_11 * gh1_29[k]
                  + pb_y[k] * gi_22[k];

        t_25[k] = f_5 * fi_22[k]
                  + f_12 * gh0_34[k]
                  - f_13 * gh1_30[k]
                  + pb_y[k] * gi_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dk0_14, dk1_14, fi_23, fk_7, fk_8, \
                         gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fi_23[k]
                  + pb_y[k] * gi_24[k];

        t_27[k] = f_3 * dk0_14[k]
                  - f_4 * dk1_14[k]
                  + pa_y[k] * fk_7[k];

        t_28[k] = pa_y[k] * fk_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fi_32, gh0_47, gh1_41, gi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fi_32[k]
                  + f_1 * gh0_47[k]
                  - f_2 * gh1_41[k]
                  + pb_z[k] * gi_26[k];
    }
}

}  // namespace simdt2ceri
