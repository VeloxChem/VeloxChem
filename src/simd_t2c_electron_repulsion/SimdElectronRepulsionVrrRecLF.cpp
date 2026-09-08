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


#include "SimdElectronRepulsionVrrRecLF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_lf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t if0, const size_t if1,
                                     const size_t kd, const size_t kf, const size_t lp0,
                                     const size_t lp1, const size_t ld, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 3.0 / p;
    const auto f_10 = 2.5 / alpha;
    const auto f_11 = 2.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_10 = buffer.data(if0 + 10);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_30 = buffer.data(if0 + 30);
    const auto *if0_36 = buffer.data(if0 + 36);
    const auto *if0_50 = buffer.data(if0 + 50);
    const auto *if0_59 = buffer.data(if0 + 59);
    const auto *if0_60 = buffer.data(if0 + 60);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_80 = buffer.data(if0 + 80);
    const auto *if0_90 = buffer.data(if0 + 90);
    const auto *if0_99 = buffer.data(if0 + 99);
    const auto *if0_100 = buffer.data(if0 + 100);
    const auto *if0_106 = buffer.data(if0 + 106);
    const auto *if0_120 = buffer.data(if0 + 120);
    const auto *if0_126 = buffer.data(if0 + 126);
    const auto *if0_129 = buffer.data(if0 + 129);
    const auto *if0_130 = buffer.data(if0 + 130);
    const auto *if0_140 = buffer.data(if0 + 140);
    const auto *if0_149 = buffer.data(if0 + 149);
    const auto *if0_156 = buffer.data(if0 + 156);
    const auto *if0_176 = buffer.data(if0 + 176);
    const auto *if0_179 = buffer.data(if0 + 179);
    const auto *if0_186 = buffer.data(if0 + 186);
    const auto *if0_189 = buffer.data(if0 + 189);
    const auto *if0_209 = buffer.data(if0 + 209);
    const auto *if0_216 = buffer.data(if0 + 216);
    const auto *if0_226 = buffer.data(if0 + 226);
    const auto *if0_236 = buffer.data(if0 + 236);
    const auto *if0_239 = buffer.data(if0 + 239);
    const auto *if0_246 = buffer.data(if0 + 246);
    const auto *if0_249 = buffer.data(if0 + 249);
    const auto *if0_256 = buffer.data(if0 + 256);
    const auto *if0_259 = buffer.data(if0 + 259);
    const auto *if0_269 = buffer.data(if0 + 269);
    const auto *if0_279 = buffer.data(if0 + 279);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_10 = buffer.data(if1 + 10);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_30 = buffer.data(if1 + 30);
    const auto *if1_36 = buffer.data(if1 + 36);
    const auto *if1_50 = buffer.data(if1 + 50);
    const auto *if1_59 = buffer.data(if1 + 59);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_80 = buffer.data(if1 + 80);
    const auto *if1_90 = buffer.data(if1 + 90);
    const auto *if1_99 = buffer.data(if1 + 99);
    const auto *if1_100 = buffer.data(if1 + 100);
    const auto *if1_106 = buffer.data(if1 + 106);
    const auto *if1_120 = buffer.data(if1 + 120);
    const auto *if1_126 = buffer.data(if1 + 126);
    const auto *if1_129 = buffer.data(if1 + 129);
    const auto *if1_130 = buffer.data(if1 + 130);
    const auto *if1_140 = buffer.data(if1 + 140);
    const auto *if1_149 = buffer.data(if1 + 149);
    const auto *if1_156 = buffer.data(if1 + 156);
    const auto *if1_176 = buffer.data(if1 + 176);
    const auto *if1_179 = buffer.data(if1 + 179);
    const auto *if1_186 = buffer.data(if1 + 186);
    const auto *if1_189 = buffer.data(if1 + 189);
    const auto *if1_209 = buffer.data(if1 + 209);
    const auto *if1_216 = buffer.data(if1 + 216);
    const auto *if1_226 = buffer.data(if1 + 226);
    const auto *if1_236 = buffer.data(if1 + 236);
    const auto *if1_239 = buffer.data(if1 + 239);
    const auto *if1_246 = buffer.data(if1 + 246);
    const auto *if1_249 = buffer.data(if1 + 249);
    const auto *if1_256 = buffer.data(if1 + 256);
    const auto *if1_259 = buffer.data(if1 + 259);
    const auto *if1_269 = buffer.data(if1 + 269);
    const auto *if1_279 = buffer.data(if1 + 279);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_215 = buffer.data(kd + 215);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_359 = buffer.data(kf + 359);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_28 = buffer.data(lp0 + 28);
    const auto *lp0_32 = buffer.data(lp0 + 32);
    const auto *lp0_43 = buffer.data(lp0 + 43);
    const auto *lp0_47 = buffer.data(lp0 + 47);
    const auto *lp0_61 = buffer.data(lp0 + 61);
    const auto *lp0_65 = buffer.data(lp0 + 65);
    const auto *lp0_82 = buffer.data(lp0 + 82);
    const auto *lp0_108 = buffer.data(lp0 + 108);
    const auto *lp0_109 = buffer.data(lp0 + 109);
    const auto *lp0_110 = buffer.data(lp0 + 110);
    const auto *lp0_114 = buffer.data(lp0 + 114);
    const auto *lp0_117 = buffer.data(lp0 + 117);
    const auto *lp0_120 = buffer.data(lp0 + 120);
    const auto *lp0_123 = buffer.data(lp0 + 123);
    const auto *lp0_126 = buffer.data(lp0 + 126);
    const auto *lp0_132 = buffer.data(lp0 + 132);
    const auto *lp0_133 = buffer.data(lp0 + 133);
    const auto *lp0_134 = buffer.data(lp0 + 134);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_28 = buffer.data(lp1 + 28);
    const auto *lp1_32 = buffer.data(lp1 + 32);
    const auto *lp1_43 = buffer.data(lp1 + 43);
    const auto *lp1_47 = buffer.data(lp1 + 47);
    const auto *lp1_61 = buffer.data(lp1 + 61);
    const auto *lp1_65 = buffer.data(lp1 + 65);
    const auto *lp1_82 = buffer.data(lp1 + 82);
    const auto *lp1_108 = buffer.data(lp1 + 108);
    const auto *lp1_109 = buffer.data(lp1 + 109);
    const auto *lp1_110 = buffer.data(lp1 + 110);
    const auto *lp1_114 = buffer.data(lp1 + 114);
    const auto *lp1_117 = buffer.data(lp1 + 117);
    const auto *lp1_120 = buffer.data(lp1 + 120);
    const auto *lp1_123 = buffer.data(lp1 + 123);
    const auto *lp1_126 = buffer.data(lp1 + 126);
    const auto *lp1_132 = buffer.data(lp1 + 132);
    const auto *lp1_133 = buffer.data(lp1 + 133);
    const auto *lp1_134 = buffer.data(lp1 + 134);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_122 = buffer.data(ld + 122);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_127 = buffer.data(ld + 127);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_164 = buffer.data(ld + 164);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_169 = buffer.data(ld + 169);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_173 = buffer.data(ld + 173);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_212 = buffer.data(ld + 212);
    const auto *ld_213 = buffer.data(ld + 213);
    const auto *ld_215 = buffer.data(ld + 215);
    const auto *ld_216 = buffer.data(ld + 216);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_220 = buffer.data(ld + 220);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_222 = buffer.data(ld + 222);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_258 = buffer.data(ld + 258);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_264 = buffer.data(ld + 264);
    const auto *ld_267 = buffer.data(ld + 267);
    const auto *ld_268 = buffer.data(ld + 268);
    const auto *ld_269 = buffer.data(ld + 269);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, kd_3, lp0_0, lp1_0, \
                         ld_0, ld_2, ld_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_0 * kd_3[k]
                 + pb_x[k] * ld_3[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, kd_5, lp0_1, lp0_2, lp1_1, \
                         lp1_2, ld_3, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * kd_5[k]
                 + pb_x[k] * ld_5[k];

        t_6[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_3[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = pb_y[k] * ld_5[k];

        t_9[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, kd_0, kd_9, \
                         kf_0, ld_6, ld_7, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * kf_0[k];

        t_11[k] = f_3 * kd_0[k]
                  + pb_y[k] * ld_6[k];

        t_12[k] = pb_z[k] * ld_6[k];

        t_13[k] = f_4 * kd_9[k]
                  + pb_x[k] * ld_9[k];

        t_14[k] = pb_z[k] * ld_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, kd_3, kd_5, kf_5, \
                         kf_6, kf_9, ld_9, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * kf_5[k];

        t_16[k] = f_5 * kd_3[k]
                  + pa_y[k] * kf_6[k];

        t_17[k] = pb_z[k] * ld_9[k];

        t_18[k] = f_3 * kd_5[k]
                  + pb_y[k] * ld_11[k];

        t_19[k] = pa_y[k] * kf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, kd_0, kf_0, kf_3, \
                         ld_12, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * kf_0[k];

        t_21[k] = pb_y[k] * ld_12[k];

        t_22[k] = f_3 * kd_0[k]
                  + pb_z[k] * ld_12[k];

        t_23[k] = pa_z[k] * kf_3[k];

        t_24[k] = pb_y[k] * ld_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, kd_3, kd_5, \
                         kd_17, kf_6, kf_9, ld_15, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * kd_17[k]
                  + pb_x[k] * ld_17[k];

        t_26[k] = pa_z[k] * kf_6[k];

        t_27[k] = f_3 * kd_3[k]
                  + pb_z[k] * ld_15[k];

        t_28[k] = pb_y[k] * ld_17[k];

        t_29[k] = f_5 * kd_5[k]
                  + pa_z[k] * kf_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, if0_0, if1_0, kd_6, \
                         kd_21, kf_10, ld_18, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * if0_0[k]
                  - f_7 * if1_0[k]
                  + pa_y[k] * kf_10[k];

        t_31[k] = f_8 * kd_6[k]
                  + pb_y[k] * ld_18[k];

        t_32[k] = pb_z[k] * ld_18[k];

        t_33[k] = f_9 * kd_21[k]
                  + pb_x[k] * ld_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, if0_36, if1_36, kd_23, \
                         kf_36, ld_19, ld_21, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * ld_19[k];

        t_35[k] = f_9 * kd_23[k]
                  + pb_x[k] * ld_23[k];

        t_36[k] = f_10 * if0_36[k]
                  - f_11 * if1_36[k]
                  + pa_x[k] * kf_36[k];

        t_37[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, kd_11, kf_11, \
                         kf_20, kf_22, lp0_11, lp1_11, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * kd_11[k]
                  + pb_y[k] * ld_23[k];

        t_39[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_23[k];

        t_40[k] = pa_y[k] * kf_20[k];

        t_41[k] = pa_z[k] * kf_11[k];

        t_42[k] = pa_y[k] * kf_22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, kd_9, kd_28, \
                         kf_13, kf_16, kf_25, ld_27, ld_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * kf_13[k];

        t_44[k] = f_9 * kd_28[k]
                  + pb_x[k] * ld_28[k];

        t_45[k] = pa_y[k] * kf_25[k];

        t_46[k] = pa_z[k] * kf_16[k];

        t_47[k] = f_3 * kd_9[k]
                  + pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, if0_0, if1_0, kd_17, kf_20, \
                         kf_29, ld_29, ld_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * kd_17[k]
                  + pb_y[k] * ld_29[k];

        t_49[k] = pa_y[k] * kf_29[k];

        t_50[k] = f_6 * if0_0[k]
                  - f_7 * if1_0[k]
                  + pa_z[k] * kf_20[k];

        t_51[k] = pb_y[k] * ld_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, kd_12, kd_33, kd_35, ld_30, \
                         ld_32, ld_33, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * kd_12[k]
                  + pb_z[k] * ld_30[k];

        t_53[k] = f_9 * kd_33[k]
                  + pb_x[k] * ld_33[k];

        t_54[k] = pb_y[k] * ld_32[k];

        t_55[k] = f_9 * kd_35[k]
                  + pb_x[k] * ld_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, if0_59, if1_59, kd_15, \
                         kf_59, lp0_16, lp1_16, ld_33, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_y[k] * ld_33[k];

        t_57[k] = f_8 * kd_15[k]
                  + pb_z[k] * ld_33[k];

        t_58[k] = pb_y[k] * ld_35[k];

        t_59[k] = f_10 * if0_59[k]
                  - f_11 * if1_59[k]
                  + pa_x[k] * kf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, if0_10, if1_10, \
                         kd_18, kd_39, kf_30, ld_36, ld_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * if0_10[k]
                  - f_13 * if1_10[k]
                  + pa_y[k] * kf_30[k];

        t_61[k] = f_5 * kd_18[k]
                  + pb_y[k] * ld_36[k];

        t_62[k] = pb_z[k] * ld_36[k];

        t_63[k] = f_14 * kd_39[k]
                  + pb_x[k] * ld_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, if0_66, if1_66, kd_41, \
                         kf_66, ld_37, ld_39, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * ld_37[k];

        t_65[k] = f_14 * kd_41[k]
                  + pb_x[k] * ld_41[k];

        t_66[k] = f_15 * if0_66[k]
                  - f_16 * if1_66[k]
                  + pa_x[k] * kf_66[k];

        t_67[k] = pb_z[k] * ld_39[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, kd_18, kd_23, kf_30, \
                         kf_31, lp0_20, lp1_20, ld_41, ld_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * kd_23[k]
                  + pb_y[k] * ld_41[k];

        t_69[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_z[k] * ld_41[k];

        t_70[k] = pa_z[k] * kf_30[k];

        t_71[k] = pa_z[k] * kf_31[k];

        t_72[k] = f_3 * kd_18[k]
                  + pb_z[k] * ld_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, kd_21, kd_46, kd_47, \
                         kf_33, kf_36, ld_45, ld_46, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * kf_33[k];

        t_74[k] = f_14 * kd_46[k]
                  + pb_x[k] * ld_46[k];

        t_75[k] = f_14 * kd_47[k]
                  + pb_x[k] * ld_47[k];

        t_76[k] = pa_z[k] * kf_36[k];

        t_77[k] = f_3 * kd_21[k]
                  + pb_z[k] * ld_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, kd_23, kd_29, kd_30, \
                         kf_39, kf_50, kf_52, ld_47, ld_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * kd_29[k]
                  + pb_y[k] * ld_47[k];

        t_79[k] = f_5 * kd_23[k]
                  + pa_z[k] * kf_39[k];

        t_80[k] = pa_y[k] * kf_50[k];

        t_81[k] = f_3 * kd_30[k]
                  + pb_y[k] * ld_48[k];

        t_82[k] = pa_y[k] * kf_52[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, kd_27, kd_33, kd_51, \
                         kd_52, kf_55, kf_56, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_14 * kd_51[k]
                  + pb_x[k] * ld_51[k];

        t_84[k] = f_14 * kd_52[k]
                  + pb_x[k] * ld_52[k];

        t_85[k] = pa_y[k] * kf_55[k];

        t_86[k] = f_5 * kd_33[k]
                  + pa_y[k] * kf_56[k];

        t_87[k] = f_8 * kd_27[k]
                  + pb_z[k] * ld_51[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, if0_20, if1_20, kd_35, \
                         kf_50, kf_59, ld_53, ld_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * kd_35[k]
                  + pb_y[k] * ld_53[k];

        t_89[k] = pa_y[k] * kf_59[k];

        t_90[k] = f_12 * if0_20[k]
                  - f_13 * if1_20[k]
                  + pa_z[k] * kf_50[k];

        t_91[k] = pb_y[k] * ld_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, kd_30, kd_57, kd_59, ld_54, \
                         ld_56, ld_57, ld_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * kd_30[k]
                  + pb_z[k] * ld_54[k];

        t_93[k] = f_14 * kd_57[k]
                  + pb_x[k] * ld_57[k];

        t_94[k] = pb_y[k] * ld_56[k];

        t_95[k] = f_14 * kd_59[k]
                  + pb_x[k] * ld_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, if0_99, if1_99, kd_33, \
                         kf_99, lp0_28, lp1_28, ld_57, ld_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_28[k]
                  - f_2 * lp1_28[k]
                  + pb_y[k] * ld_57[k];

        t_97[k] = f_5 * kd_33[k]
                  + pb_z[k] * ld_57[k];

        t_98[k] = pb_y[k] * ld_59[k];

        t_99[k] = f_15 * if0_99[k]
                  - f_16 * if1_99[k]
                  + pa_x[k] * kf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, pb_z, if0_30, if1_30, \
                         kd_36, kd_63, kf_60, ld_60, ld_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * if0_30[k]
                   - f_18 * if1_30[k]
                   + pa_y[k] * kf_60[k];

        t_101[k] = f_19 * kd_36[k]
                   + pb_y[k] * ld_60[k];

        t_102[k] = pb_z[k] * ld_60[k];

        t_103[k] = f_19 * kd_63[k]
                   + pb_x[k] * ld_63[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, pb_z, if0_106, if1_106, \
                         kd_65, kf_106, ld_61, ld_63, ld_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * ld_61[k];

        t_105[k] = f_19 * kd_65[k]
                   + pb_x[k] * ld_65[k];

        t_106[k] = f_17 * if0_106[k]
                   - f_18 * if1_106[k]
                   + pa_x[k] * kf_106[k];

        t_107[k] = pb_z[k] * ld_63[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, kd_36, kd_41, \
                         kf_60, kf_61, lp0_32, lp1_32, ld_65, ld_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_19 * kd_41[k]
                   + pb_y[k] * ld_65[k];

        t_109[k] = f_1 * lp0_32[k]
                   - f_2 * lp1_32[k]
                   + pb_z[k] * ld_65[k];

        t_110[k] = pa_z[k] * kf_60[k];

        t_111[k] = pa_z[k] * kf_61[k];

        t_112[k] = f_3 * kd_36[k]
                   + pb_z[k] * ld_66[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, kd_39, kd_70, \
                         kd_71, kf_63, kf_66, ld_69, ld_70, ld_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * kf_63[k];

        t_114[k] = f_19 * kd_70[k]
                   + pb_x[k] * ld_70[k];

        t_115[k] = f_19 * kd_71[k]
                   + pb_x[k] * ld_71[k];

        t_116[k] = pa_z[k] * kf_66[k];

        t_117[k] = f_3 * kd_39[k]
                   + pb_z[k] * ld_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pa_z, pb_y, if0_50, if1_50, kd_41, \
                         kd_47, kd_48, kf_69, kf_80, ld_71, ld_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * kd_47[k]
                   + pb_y[k] * ld_71[k];

        t_119[k] = f_5 * kd_41[k]
                   + pa_z[k] * kf_69[k];

        t_120[k] = f_6 * if0_50[k]
                   - f_7 * if1_50[k]
                   + pa_y[k] * kf_80[k];

        t_121[k] = f_8 * kd_48[k]
                   + pb_y[k] * ld_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pb_x, pb_z, kd_42, kd_75, kd_76, kd_77, \
                         ld_72, ld_75, ld_76, ld_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * kd_42[k]
                   + pb_z[k] * ld_72[k];

        t_123[k] = f_19 * kd_75[k]
                   + pb_x[k] * ld_75[k];

        t_124[k] = f_19 * kd_76[k]
                   + pb_x[k] * ld_76[k];

        t_125[k] = f_19 * kd_77[k]
                   + pb_x[k] * ld_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, if0_126, if1_126, kd_45, \
                         kd_53, kf_126, ld_75, ld_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_17 * if0_126[k]
                   - f_18 * if1_126[k]
                   + pa_x[k] * kf_126[k];

        t_127[k] = f_8 * kd_45[k]
                   + pb_z[k] * ld_75[k];

        t_128[k] = f_8 * kd_53[k]
                   + pb_y[k] * ld_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, if0_129, if1_129, \
                         kd_54, kf_90, kf_92, kf_129, ld_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_17 * if0_129[k]
                   - f_18 * if1_129[k]
                   + pa_x[k] * kf_129[k];

        t_130[k] = pa_y[k] * kf_90[k];

        t_131[k] = f_3 * kd_54[k]
                   + pb_y[k] * ld_78[k];

        t_132[k] = pa_y[k] * kf_92[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pa_y, pb_x, pb_z, kd_51, kd_57, \
                         kd_81, kd_82, kf_95, kf_96, ld_81, ld_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_19 * kd_81[k]
                   + pb_x[k] * ld_81[k];

        t_134[k] = f_19 * kd_82[k]
                   + pb_x[k] * ld_82[k];

        t_135[k] = pa_y[k] * kf_95[k];

        t_136[k] = f_5 * kd_57[k]
                   + pa_y[k] * kf_96[k];

        t_137[k] = f_5 * kd_51[k]
                   + pb_z[k] * ld_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pa_z, pb_y, if0_50, if1_50, kd_59, \
                         kf_90, kf_99, ld_83, ld_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * kd_59[k]
                   + pb_y[k] * ld_83[k];

        t_139[k] = pa_y[k] * kf_99[k];

        t_140[k] = f_17 * if0_50[k]
                   - f_18 * if1_50[k]
                   + pa_z[k] * kf_90[k];

        t_141[k] = pb_y[k] * ld_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, kd_54, kd_87, kd_89, \
                         ld_84, ld_86, ld_87, ld_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_19 * kd_54[k]
                   + pb_z[k] * ld_84[k];

        t_143[k] = f_19 * kd_87[k]
                   + pb_x[k] * ld_87[k];

        t_144[k] = pb_y[k] * ld_86[k];

        t_145[k] = f_19 * kd_89[k]
                   + pb_x[k] * ld_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, pb_z, if0_149, if1_149, \
                         kd_57, kf_149, lp0_43, lp1_43, ld_87, ld_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * lp0_43[k]
                   - f_2 * lp1_43[k]
                   + pb_y[k] * ld_87[k];

        t_147[k] = f_19 * kd_57[k]
                   + pb_z[k] * ld_87[k];

        t_148[k] = pb_y[k] * ld_89[k];

        t_149[k] = f_17 * if0_149[k]
                   - f_18 * if1_149[k]
                   + pa_x[k] * kf_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pb_x, pb_y, pb_z, if0_60, if1_60, \
                         kd_60, kd_93, kf_100, ld_90, ld_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_15 * if0_60[k]
                   - f_16 * if1_60[k]
                   + pa_y[k] * kf_100[k];

        t_151[k] = f_14 * kd_60[k]
                   + pb_y[k] * ld_90[k];

        t_152[k] = pb_z[k] * ld_90[k];

        t_153[k] = f_5 * kd_93[k]
                   + pb_x[k] * ld_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, if0_156, if1_156, \
                         kd_95, kf_156, ld_91, ld_93, ld_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * ld_91[k];

        t_155[k] = f_5 * kd_95[k]
                   + pb_x[k] * ld_95[k];

        t_156[k] = f_12 * if0_156[k]
                   - f_13 * if1_156[k]
                   + pa_x[k] * kf_156[k];

        t_157[k] = pb_z[k] * ld_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_y, pb_z, kd_60, kd_65, \
                         kf_100, kf_101, lp0_47, lp1_47, ld_95, ld_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_14 * kd_65[k]
                   + pb_y[k] * ld_95[k];

        t_159[k] = f_1 * lp0_47[k]
                   - f_2 * lp1_47[k]
                   + pb_z[k] * ld_95[k];

        t_160[k] = pa_z[k] * kf_100[k];

        t_161[k] = pa_z[k] * kf_101[k];

        t_162[k] = f_3 * kd_60[k]
                   + pb_z[k] * ld_96[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_z, pb_x, pb_z, kd_63, kd_100, \
                         kd_101, kf_103, kf_106, ld_99, ld_100, \
                         ld_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * kf_103[k];

        t_164[k] = f_5 * kd_100[k]
                   + pb_x[k] * ld_100[k];

        t_165[k] = f_5 * kd_101[k]
                   + pb_x[k] * ld_101[k];

        t_166[k] = pa_z[k] * kf_106[k];

        t_167[k] = f_3 * kd_63[k]
                   + pb_z[k] * ld_99[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, if0_80, if1_80, kd_65, \
                         kd_71, kd_72, kf_109, kf_120, ld_101, ld_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * kd_71[k]
                   + pb_y[k] * ld_101[k];

        t_169[k] = f_5 * kd_65[k]
                   + pa_z[k] * kf_109[k];

        t_170[k] = f_12 * if0_80[k]
                   - f_13 * if1_80[k]
                   + pa_y[k] * kf_120[k];

        t_171[k] = f_5 * kd_72[k]
                   + pb_y[k] * ld_102[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, pb_z, kd_66, kd_105, kd_106, \
                         kd_107, ld_102, ld_105, ld_106, ld_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_8 * kd_66[k]
                   + pb_z[k] * ld_102[k];

        t_173[k] = f_5 * kd_105[k]
                   + pb_x[k] * ld_105[k];

        t_174[k] = f_5 * kd_106[k]
                   + pb_x[k] * ld_106[k];

        t_175[k] = f_5 * kd_107[k]
                   + pb_x[k] * ld_107[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, pb_y, pb_z, if0_176, if1_176, kd_69, \
                         kd_77, kf_176, ld_105, ld_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_12 * if0_176[k]
                   - f_13 * if1_176[k]
                   + pa_x[k] * kf_176[k];

        t_177[k] = f_8 * kd_69[k]
                   + pb_z[k] * ld_105[k];

        t_178[k] = f_5 * kd_77[k]
                   + pb_y[k] * ld_107[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pb_y, if0_90, if0_179, if1_90, \
                         if1_179, kd_78, kf_130, kf_179, ld_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_12 * if0_179[k]
                   - f_13 * if1_179[k]
                   + pa_x[k] * kf_179[k];

        t_180[k] = f_6 * if0_90[k]
                   - f_7 * if1_90[k]
                   + pa_y[k] * kf_130[k];

        t_181[k] = f_8 * kd_78[k]
                   + pb_y[k] * ld_108[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_x, pb_z, kd_72, kd_111, kd_112, \
                         kd_113, ld_108, ld_111, ld_112, ld_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_5 * kd_72[k]
                   + pb_z[k] * ld_108[k];

        t_183[k] = f_5 * kd_111[k]
                   + pb_x[k] * ld_111[k];

        t_184[k] = f_5 * kd_112[k]
                   + pb_x[k] * ld_112[k];

        t_185[k] = f_5 * kd_113[k]
                   + pb_x[k] * ld_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_x, pb_y, pb_z, if0_186, if1_186, kd_75, \
                         kd_83, kf_186, ld_111, ld_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_12 * if0_186[k]
                   - f_13 * if1_186[k]
                   + pa_x[k] * kf_186[k];

        t_187[k] = f_5 * kd_75[k]
                   + pb_z[k] * ld_111[k];

        t_188[k] = f_8 * kd_83[k]
                   + pb_y[k] * ld_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pa_y, pb_y, if0_189, if1_189, \
                         kd_84, kf_140, kf_142, kf_189, ld_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * if0_189[k]
                   - f_13 * if1_189[k]
                   + pa_x[k] * kf_189[k];

        t_190[k] = pa_y[k] * kf_140[k];

        t_191[k] = f_3 * kd_84[k]
                   + pb_y[k] * ld_114[k];

        t_192[k] = pa_y[k] * kf_142[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_z, kd_81, kd_87, \
                         kd_117, kd_118, kf_145, kf_146, ld_117, \
                         ld_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_5 * kd_117[k]
                   + pb_x[k] * ld_117[k];

        t_194[k] = f_5 * kd_118[k]
                   + pb_x[k] * ld_118[k];

        t_195[k] = pa_y[k] * kf_145[k];

        t_196[k] = f_5 * kd_87[k]
                   + pa_y[k] * kf_146[k];

        t_197[k] = f_19 * kd_81[k]
                   + pb_z[k] * ld_117[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_y, pa_z, pb_y, if0_90, if1_90, kd_89, \
                         kf_140, kf_149, ld_119, ld_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * kd_89[k]
                   + pb_y[k] * ld_119[k];

        t_199[k] = pa_y[k] * kf_149[k];

        t_200[k] = f_15 * if0_90[k]
                   - f_16 * if1_90[k]
                   + pa_z[k] * kf_140[k];

        t_201[k] = pb_y[k] * ld_120[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, kd_84, kd_123, kd_125, \
                         ld_120, ld_122, ld_123, ld_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_14 * kd_84[k]
                   + pb_z[k] * ld_120[k];

        t_203[k] = f_5 * kd_123[k]
                   + pb_x[k] * ld_123[k];

        t_204[k] = pb_y[k] * ld_122[k];

        t_205[k] = f_5 * kd_125[k]
                   + pb_x[k] * ld_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pb_y, pb_z, if0_209, if1_209, \
                         kd_87, kf_209, lp0_61, lp1_61, ld_123, \
                         ld_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * lp0_61[k]
                   - f_2 * lp1_61[k]
                   + pb_y[k] * ld_123[k];

        t_207[k] = f_14 * kd_87[k]
                   + pb_z[k] * ld_123[k];

        t_208[k] = pb_y[k] * ld_125[k];

        t_209[k] = f_12 * if0_209[k]
                   - f_13 * if1_209[k]
                   + pa_x[k] * kf_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_y, pb_x, pb_y, pb_z, if0_100, if1_100, \
                         kd_90, kd_129, kf_150, ld_126, ld_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * if0_100[k]
                   - f_11 * if1_100[k]
                   + pa_y[k] * kf_150[k];

        t_211[k] = f_9 * kd_90[k]
                   + pb_y[k] * ld_126[k];

        t_212[k] = pb_z[k] * ld_126[k];

        t_213[k] = f_8 * kd_129[k]
                   + pb_x[k] * ld_129[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_x, pb_x, pb_z, if0_216, if1_216, \
                         kd_131, kf_216, ld_127, ld_129, ld_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_z[k] * ld_127[k];

        t_215[k] = f_8 * kd_131[k]
                   + pb_x[k] * ld_131[k];

        t_216[k] = f_6 * if0_216[k]
                   - f_7 * if1_216[k]
                   + pa_x[k] * kf_216[k];

        t_217[k] = pb_z[k] * ld_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, pa_z, pb_y, pb_z, kd_90, kd_95, \
                         kf_150, kf_151, lp0_65, lp1_65, ld_131, \
                         ld_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_9 * kd_95[k]
                   + pb_y[k] * ld_131[k];

        t_219[k] = f_1 * lp0_65[k]
                   - f_2 * lp1_65[k]
                   + pb_z[k] * ld_131[k];

        t_220[k] = pa_z[k] * kf_150[k];

        t_221[k] = pa_z[k] * kf_151[k];

        t_222[k] = f_3 * kd_90[k]
                   + pb_z[k] * ld_132[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, kd_93, kd_136, \
                         kd_137, kf_153, kf_156, ld_135, ld_136, \
                         ld_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pa_z[k] * kf_153[k];

        t_224[k] = f_8 * kd_136[k]
                   + pb_x[k] * ld_136[k];

        t_225[k] = f_8 * kd_137[k]
                   + pb_x[k] * ld_137[k];

        t_226[k] = pa_z[k] * kf_156[k];

        t_227[k] = f_3 * kd_93[k]
                   + pb_z[k] * ld_135[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, pa_z, pb_y, if0_120, if1_120, \
                         kd_95, kd_101, kd_102, kf_159, kf_170, ld_137, \
                         ld_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_14 * kd_101[k]
                   + pb_y[k] * ld_137[k];

        t_229[k] = f_5 * kd_95[k]
                   + pa_z[k] * kf_159[k];

        t_230[k] = f_17 * if0_120[k]
                   - f_18 * if1_120[k]
                   + pa_y[k] * kf_170[k];

        t_231[k] = f_19 * kd_102[k]
                   + pb_y[k] * ld_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pb_x, pb_z, kd_96, kd_141, kd_142, \
                         kd_143, ld_138, ld_141, ld_142, ld_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * kd_96[k]
                   + pb_z[k] * ld_138[k];

        t_233[k] = f_8 * kd_141[k]
                   + pb_x[k] * ld_141[k];

        t_234[k] = f_8 * kd_142[k]
                   + pb_x[k] * ld_142[k];

        t_235[k] = f_8 * kd_143[k]
                   + pb_x[k] * ld_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_x, pb_y, pb_z, if0_236, if1_236, kd_99, \
                         kd_107, kf_236, ld_141, ld_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_6 * if0_236[k]
                   - f_7 * if1_236[k]
                   + pa_x[k] * kf_236[k];

        t_237[k] = f_8 * kd_99[k]
                   + pb_z[k] * ld_141[k];

        t_238[k] = f_19 * kd_107[k]
                   + pb_y[k] * ld_143[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pa_x, pa_y, pb_y, if0_130, if0_239, if1_130, \
                         if1_239, kd_108, kf_180, kf_239, ld_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_6 * if0_239[k]
                   - f_7 * if1_239[k]
                   + pa_x[k] * kf_239[k];

        t_240[k] = f_12 * if0_130[k]
                   - f_13 * if1_130[k]
                   + pa_y[k] * kf_180[k];

        t_241[k] = f_5 * kd_108[k]
                   + pb_y[k] * ld_144[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pb_x, pb_z, kd_102, kd_147, kd_148, \
                         kd_149, ld_144, ld_147, ld_148, ld_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_5 * kd_102[k]
                   + pb_z[k] * ld_144[k];

        t_243[k] = f_8 * kd_147[k]
                   + pb_x[k] * ld_147[k];

        t_244[k] = f_8 * kd_148[k]
                   + pb_x[k] * ld_148[k];

        t_245[k] = f_8 * kd_149[k]
                   + pb_x[k] * ld_149[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pb_y, pb_z, if0_246, if1_246, kd_105, \
                         kd_113, kf_246, ld_147, ld_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_6 * if0_246[k]
                   - f_7 * if1_246[k]
                   + pa_x[k] * kf_246[k];

        t_247[k] = f_5 * kd_105[k]
                   + pb_z[k] * ld_147[k];

        t_248[k] = f_5 * kd_113[k]
                   + pb_y[k] * ld_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pa_y, pb_y, if0_140, if0_249, if1_140, \
                         if1_249, kd_114, kf_190, kf_249, ld_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * if0_249[k]
                   - f_7 * if1_249[k]
                   + pa_x[k] * kf_249[k];

        t_250[k] = f_6 * if0_140[k]
                   - f_7 * if1_140[k]
                   + pa_y[k] * kf_190[k];

        t_251[k] = f_8 * kd_114[k]
                   + pb_y[k] * ld_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, pb_z, kd_108, kd_153, kd_154, \
                         kd_155, ld_150, ld_153, ld_154, ld_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_19 * kd_108[k]
                   + pb_z[k] * ld_150[k];

        t_253[k] = f_8 * kd_153[k]
                   + pb_x[k] * ld_153[k];

        t_254[k] = f_8 * kd_154[k]
                   + pb_x[k] * ld_154[k];

        t_255[k] = f_8 * kd_155[k]
                   + pb_x[k] * ld_155[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pa_x, pb_y, pb_z, if0_256, if1_256, kd_111, \
                         kd_119, kf_256, ld_153, ld_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_6 * if0_256[k]
                   - f_7 * if1_256[k]
                   + pa_x[k] * kf_256[k];

        t_257[k] = f_19 * kd_111[k]
                   + pb_z[k] * ld_153[k];

        t_258[k] = f_8 * kd_119[k]
                   + pb_y[k] * ld_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pa_x, pa_y, pb_y, if0_259, if1_259, \
                         kd_120, kf_200, kf_202, kf_259, ld_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_6 * if0_259[k]
                   - f_7 * if1_259[k]
                   + pa_x[k] * kf_259[k];

        t_260[k] = pa_y[k] * kf_200[k];

        t_261[k] = f_3 * kd_120[k]
                   + pb_y[k] * ld_156[k];

        t_262[k] = pa_y[k] * kf_202[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, pa_y, pb_x, pb_z, kd_117, kd_123, \
                         kd_159, kd_160, kf_205, kf_206, ld_159, \
                         ld_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_8 * kd_159[k]
                   + pb_x[k] * ld_159[k];

        t_264[k] = f_8 * kd_160[k]
                   + pb_x[k] * ld_160[k];

        t_265[k] = pa_y[k] * kf_205[k];

        t_266[k] = f_5 * kd_123[k]
                   + pa_y[k] * kf_206[k];

        t_267[k] = f_14 * kd_117[k]
                   + pb_z[k] * ld_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pa_z, pb_y, if0_140, if1_140, \
                         kd_125, kf_200, kf_209, ld_161, ld_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_3 * kd_125[k]
                   + pb_y[k] * ld_161[k];

        t_269[k] = pa_y[k] * kf_209[k];

        t_270[k] = f_10 * if0_140[k]
                   - f_11 * if1_140[k]
                   + pa_z[k] * kf_200[k];

        t_271[k] = pb_y[k] * ld_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, kd_120, kd_165, kd_167, \
                         ld_162, ld_164, ld_165, ld_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_9 * kd_120[k]
                   + pb_z[k] * ld_162[k];

        t_273[k] = f_8 * kd_165[k]
                   + pb_x[k] * ld_165[k];

        t_274[k] = pb_y[k] * ld_164[k];

        t_275[k] = f_8 * kd_167[k]
                   + pb_x[k] * ld_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pb_y, pb_z, if0_279, if1_279, \
                         kd_123, kf_279, lp0_82, lp1_82, ld_165, \
                         ld_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * lp0_82[k]
                   - f_2 * lp1_82[k]
                   + pb_y[k] * ld_165[k];

        t_277[k] = f_9 * kd_123[k]
                   + pb_z[k] * ld_165[k];

        t_278[k] = pb_y[k] * ld_167[k];

        t_279[k] = f_6 * if0_279[k]
                   - f_7 * if1_279[k]
                   + pa_x[k] * kf_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pa_x, pb_x, pb_y, pb_z, kd_126, \
                         kd_168, kd_171, kf_280, ld_168, ld_169, \
                         ld_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_5 * kd_168[k]
                   + pa_x[k] * kf_280[k];

        t_281[k] = f_4 * kd_126[k]
                   + pb_y[k] * ld_168[k];

        t_282[k] = pb_z[k] * ld_168[k];

        t_283[k] = f_3 * kd_171[k]
                   + pb_x[k] * ld_171[k];

        t_284[k] = pb_z[k] * ld_169[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, pa_x, pb_x, pb_z, kd_173, kf_286, \
                         kf_288, kf_289, ld_171, ld_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_3 * kd_173[k]
                   + pb_x[k] * ld_173[k];

        t_286[k] = pa_x[k] * kf_286[k];

        t_287[k] = pb_z[k] * ld_171[k];

        t_288[k] = pa_x[k] * kf_288[k];

        t_289[k] = pa_x[k] * kf_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pa_z, pb_x, pb_z, kd_126, kd_178, \
                         kf_210, kf_211, kf_213, ld_174, ld_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_z[k] * kf_210[k];

        t_291[k] = pa_z[k] * kf_211[k];

        t_292[k] = f_3 * kd_126[k]
                   + pb_z[k] * ld_174[k];

        t_293[k] = pa_z[k] * kf_213[k];

        t_294[k] = f_3 * kd_178[k]
                   + pb_x[k] * ld_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, t_300, pa_x, pb_x, kd_179, kd_180, \
                         kf_296, kf_297, kf_298, kf_299, kf_300, \
                         ld_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_3 * kd_179[k]
                   + pb_x[k] * ld_179[k];

        t_296[k] = pa_x[k] * kf_296[k];

        t_297[k] = pa_x[k] * kf_297[k];

        t_298[k] = pa_x[k] * kf_298[k];

        t_299[k] = pa_x[k] * kf_299[k];

        t_300[k] = f_5 * kd_180[k]
                   + pa_x[k] * kf_300[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_x, pb_y, pb_z, kd_132, kd_138, kd_183, \
                         kd_184, ld_180, ld_183, ld_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * kd_138[k]
                   + pb_y[k] * ld_180[k];

        t_302[k] = f_8 * kd_132[k]
                   + pb_z[k] * ld_180[k];

        t_303[k] = f_3 * kd_183[k]
                   + pb_x[k] * ld_183[k];

        t_304[k] = f_3 * kd_184[k]
                   + pb_x[k] * ld_184[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, t_310, pa_x, pb_x, kd_185, kd_186, \
                         kf_306, kf_307, kf_308, kf_309, kf_310, \
                         ld_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_3 * kd_185[k]
                   + pb_x[k] * ld_185[k];

        t_306[k] = pa_x[k] * kf_306[k];

        t_307[k] = pa_x[k] * kf_307[k];

        t_308[k] = pa_x[k] * kf_308[k];

        t_309[k] = pa_x[k] * kf_309[k];

        t_310[k] = f_5 * kd_186[k]
                   + pa_x[k] * kf_310[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_x, pb_y, pb_z, kd_138, kd_144, kd_189, \
                         kd_190, ld_186, ld_189, ld_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_19 * kd_144[k]
                   + pb_y[k] * ld_186[k];

        t_312[k] = f_5 * kd_138[k]
                   + pb_z[k] * ld_186[k];

        t_313[k] = f_3 * kd_189[k]
                   + pb_x[k] * ld_189[k];

        t_314[k] = f_3 * kd_190[k]
                   + pb_x[k] * ld_190[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, pa_x, pb_x, kd_191, kd_192, \
                         kf_316, kf_317, kf_318, kf_319, kf_320, \
                         ld_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_3 * kd_191[k]
                   + pb_x[k] * ld_191[k];

        t_316[k] = pa_x[k] * kf_316[k];

        t_317[k] = pa_x[k] * kf_317[k];

        t_318[k] = pa_x[k] * kf_318[k];

        t_319[k] = pa_x[k] * kf_319[k];

        t_320[k] = f_5 * kd_192[k]
                   + pa_x[k] * kf_320[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_y, pb_z, kd_144, kd_150, kd_195, \
                         kd_196, ld_192, ld_195, ld_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_5 * kd_150[k]
                   + pb_y[k] * ld_192[k];

        t_322[k] = f_19 * kd_144[k]
                   + pb_z[k] * ld_192[k];

        t_323[k] = f_3 * kd_195[k]
                   + pb_x[k] * ld_195[k];

        t_324[k] = f_3 * kd_196[k]
                   + pb_x[k] * ld_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, t_330, pa_x, pb_x, kd_197, kd_198, \
                         kf_326, kf_327, kf_328, kf_329, kf_330, \
                         ld_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_3 * kd_197[k]
                   + pb_x[k] * ld_197[k];

        t_326[k] = pa_x[k] * kf_326[k];

        t_327[k] = pa_x[k] * kf_327[k];

        t_328[k] = pa_x[k] * kf_328[k];

        t_329[k] = pa_x[k] * kf_329[k];

        t_330[k] = f_5 * kd_198[k]
                   + pa_x[k] * kf_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_x, pb_y, pb_z, kd_150, kd_156, kd_201, \
                         kd_202, ld_198, ld_201, ld_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_8 * kd_156[k]
                   + pb_y[k] * ld_198[k];

        t_332[k] = f_14 * kd_150[k]
                   + pb_z[k] * ld_198[k];

        t_333[k] = f_3 * kd_201[k]
                   + pb_x[k] * ld_201[k];

        t_334[k] = f_3 * kd_202[k]
                   + pb_x[k] * ld_202[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, t_340, pa_x, pa_y, pb_x, kd_203, \
                         kf_270, kf_336, kf_337, kf_338, kf_339, \
                         ld_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_3 * kd_203[k]
                   + pb_x[k] * ld_203[k];

        t_336[k] = pa_x[k] * kf_336[k];

        t_337[k] = pa_x[k] * kf_337[k];

        t_338[k] = pa_x[k] * kf_338[k];

        t_339[k] = pa_x[k] * kf_339[k];

        t_340[k] = pa_y[k] * kf_270[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, pa_y, pb_x, pb_y, kd_162, kd_207, \
                         kd_208, kf_272, kf_275, ld_204, ld_207, \
                         ld_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_3 * kd_162[k]
                   + pb_y[k] * ld_204[k];

        t_342[k] = pa_y[k] * kf_272[k];

        t_343[k] = f_3 * kd_207[k]
                   + pb_x[k] * ld_207[k];

        t_344[k] = f_3 * kd_208[k]
                   + pb_x[k] * ld_208[k];

        t_345[k] = pa_y[k] * kf_275[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, t_351, pa_x, pb_y, kd_210, kf_346, \
                         kf_347, kf_348, kf_349, kf_350, ld_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pa_x[k] * kf_346[k];

        t_347[k] = pa_x[k] * kf_347[k];

        t_348[k] = pa_x[k] * kf_348[k];

        t_349[k] = pa_x[k] * kf_349[k];

        t_350[k] = f_5 * kd_210[k]
                   + pa_x[k] * kf_350[k];

        t_351[k] = pb_y[k] * ld_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_x, pb_y, pb_z, kd_162, kd_213, kd_215, \
                         ld_210, ld_212, ld_213, ld_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_4 * kd_162[k]
                   + pb_z[k] * ld_210[k];

        t_353[k] = f_3 * kd_213[k]
                   + pb_x[k] * ld_213[k];

        t_354[k] = pb_y[k] * ld_212[k];

        t_355[k] = f_3 * kd_215[k]
                   + pb_x[k] * ld_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, t_360, pa_x, pb_x, pb_y, kf_356, kf_357, \
                         kf_359, lp0_108, lp1_108, ld_215, ld_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * kf_356[k];

        t_357[k] = pa_x[k] * kf_357[k];

        t_358[k] = pb_y[k] * ld_215[k];

        t_359[k] = pa_x[k] * kf_359[k];

        t_360[k] = f_1 * lp0_108[k]
                   - f_2 * lp1_108[k]
                   + pb_x[k] * ld_216[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, t_365, pb_x, pb_y, pb_z, kd_168, ld_216, \
                         ld_219, ld_220, ld_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_0 * kd_168[k]
                   + pb_y[k] * ld_216[k];

        t_362[k] = pb_z[k] * ld_216[k];

        t_363[k] = pb_x[k] * ld_219[k];

        t_364[k] = pb_x[k] * ld_220[k];

        t_365[k] = pb_x[k] * ld_221[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pb_z, kd_171, kd_173, lp0_109, \
                         lp0_110, lp1_109, lp1_110, ld_219, ld_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_0 * kd_171[k]
                   + f_1 * lp0_109[k]
                   - f_2 * lp1_109[k]
                   + pb_y[k] * ld_219[k];

        t_367[k] = pb_z[k] * ld_219[k];

        t_368[k] = f_0 * kd_173[k]
                   + pb_y[k] * ld_221[k];

        t_369[k] = f_1 * lp0_110[k]
                   - f_2 * lp1_110[k]
                   + pb_z[k] * ld_221[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, t_375, pa_z, pb_x, pb_z, kd_168, \
                         kf_280, kf_281, ld_222, ld_225, ld_226, \
                         ld_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_z[k] * kf_280[k];

        t_371[k] = pa_z[k] * kf_281[k];

        t_372[k] = f_3 * kd_168[k]
                   + pb_z[k] * ld_222[k];

        t_373[k] = pb_x[k] * ld_225[k];

        t_374[k] = pb_x[k] * ld_226[k];

        t_375[k] = pb_x[k] * ld_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pb_y, pb_z, kd_171, kd_173, kd_179, \
                         kf_286, kf_289, ld_225, ld_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_z[k] * kf_286[k];

        t_377[k] = f_3 * kd_171[k]
                   + pb_z[k] * ld_225[k];

        t_378[k] = f_4 * kd_179[k]
                   + pb_y[k] * ld_227[k];

        t_379[k] = f_5 * kd_173[k]
                   + pa_z[k] * kf_289[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, pb_x, pb_y, pb_z, kd_174, kd_180, \
                         lp0_114, lp1_114, ld_228, ld_231, ld_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_1 * lp0_114[k]
                   - f_2 * lp1_114[k]
                   + pb_x[k] * ld_228[k];

        t_381[k] = f_9 * kd_180[k]
                   + pb_y[k] * ld_228[k];

        t_382[k] = f_8 * kd_174[k]
                   + pb_z[k] * ld_228[k];

        t_383[k] = pb_x[k] * ld_231[k];

        t_384[k] = pb_x[k] * ld_232[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_z, pb_x, pb_y, pb_z, if0_216, if1_216, \
                         kd_177, kd_185, kf_296, ld_231, ld_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_x[k] * ld_233[k];

        t_386[k] = f_6 * if0_216[k]
                   - f_7 * if1_216[k]
                   + pa_z[k] * kf_296[k];

        t_387[k] = f_8 * kd_177[k]
                   + pb_z[k] * ld_231[k];

        t_388[k] = f_9 * kd_185[k]
                   + pb_y[k] * ld_233[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_y, pb_x, pb_y, pb_z, if0_239, if1_239, \
                         kd_180, kd_186, kf_309, lp0_117, lp1_117, \
                         ld_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_10 * if0_239[k]
                   - f_11 * if1_239[k]
                   + pa_y[k] * kf_309[k];

        t_390[k] = f_1 * lp0_117[k]
                   - f_2 * lp1_117[k]
                   + pb_x[k] * ld_234[k];

        t_391[k] = f_14 * kd_186[k]
                   + pb_y[k] * ld_234[k];

        t_392[k] = f_5 * kd_180[k]
                   + pb_z[k] * ld_234[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, pa_z, pb_x, pb_z, if0_226, \
                         if1_226, kd_183, kf_306, ld_237, ld_238, \
                         ld_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pb_x[k] * ld_237[k];

        t_394[k] = pb_x[k] * ld_238[k];

        t_395[k] = pb_x[k] * ld_239[k];

        t_396[k] = f_12 * if0_226[k]
                   - f_13 * if1_226[k]
                   + pa_z[k] * kf_306[k];

        t_397[k] = f_5 * kd_183[k]
                   + pb_z[k] * ld_237[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, pa_y, pb_x, pb_y, if0_249, if1_249, \
                         kd_191, kd_192, kf_319, lp0_120, lp1_120, ld_239, \
                         ld_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_14 * kd_191[k]
                   + pb_y[k] * ld_239[k];

        t_399[k] = f_15 * if0_249[k]
                   - f_16 * if1_249[k]
                   + pa_y[k] * kf_319[k];

        t_400[k] = f_1 * lp0_120[k]
                   - f_2 * lp1_120[k]
                   + pb_x[k] * ld_240[k];

        t_401[k] = f_19 * kd_192[k]
                   + pb_y[k] * ld_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_z, pb_x, pb_z, if0_236, \
                         if1_236, kd_186, kf_316, ld_240, ld_243, ld_244, \
                         ld_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_19 * kd_186[k]
                   + pb_z[k] * ld_240[k];

        t_403[k] = pb_x[k] * ld_243[k];

        t_404[k] = pb_x[k] * ld_244[k];

        t_405[k] = pb_x[k] * ld_245[k];

        t_406[k] = f_17 * if0_236[k]
                   - f_18 * if1_236[k]
                   + pa_z[k] * kf_316[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_y, pb_y, pb_z, if0_259, if1_259, kd_189, \
                         kd_197, kf_329, ld_243, ld_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_19 * kd_189[k]
                   + pb_z[k] * ld_243[k];

        t_408[k] = f_19 * kd_197[k]
                   + pb_y[k] * ld_245[k];

        t_409[k] = f_17 * if0_259[k]
                   - f_18 * if1_259[k]
                   + pa_y[k] * kf_329[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pb_x, pb_y, pb_z, kd_192, kd_198, \
                         lp0_123, lp1_123, ld_246, ld_249, ld_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_1 * lp0_123[k]
                   - f_2 * lp1_123[k]
                   + pb_x[k] * ld_246[k];

        t_411[k] = f_5 * kd_198[k]
                   + pb_y[k] * ld_246[k];

        t_412[k] = f_14 * kd_192[k]
                   + pb_z[k] * ld_246[k];

        t_413[k] = pb_x[k] * ld_249[k];

        t_414[k] = pb_x[k] * ld_250[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_z, pb_x, pb_y, pb_z, if0_246, if1_246, \
                         kd_195, kd_203, kf_326, ld_249, ld_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = pb_x[k] * ld_251[k];

        t_416[k] = f_15 * if0_246[k]
                   - f_16 * if1_246[k]
                   + pa_z[k] * kf_326[k];

        t_417[k] = f_14 * kd_195[k]
                   + pb_z[k] * ld_249[k];

        t_418[k] = f_5 * kd_203[k]
                   + pb_y[k] * ld_251[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pa_y, pb_x, pb_y, pb_z, if0_269, if1_269, \
                         kd_198, kd_204, kf_339, lp0_126, lp1_126, \
                         ld_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_12 * if0_269[k]
                   - f_13 * if1_269[k]
                   + pa_y[k] * kf_339[k];

        t_420[k] = f_1 * lp0_126[k]
                   - f_2 * lp1_126[k]
                   + pb_x[k] * ld_252[k];

        t_421[k] = f_8 * kd_204[k]
                   + pb_y[k] * ld_252[k];

        t_422[k] = f_9 * kd_198[k]
                   + pb_z[k] * ld_252[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, pa_z, pb_x, pb_z, if0_256, \
                         if1_256, kd_201, kf_336, ld_255, ld_256, \
                         ld_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = pb_x[k] * ld_255[k];

        t_424[k] = pb_x[k] * ld_256[k];

        t_425[k] = pb_x[k] * ld_257[k];

        t_426[k] = f_10 * if0_256[k]
                   - f_11 * if1_256[k]
                   + pa_z[k] * kf_336[k];

        t_427[k] = f_9 * kd_201[k]
                   + pb_z[k] * ld_255[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pa_y, pb_y, if0_279, if1_279, \
                         kd_209, kd_210, kf_349, kf_350, kf_352, ld_257, \
                         ld_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_8 * kd_209[k]
                   + pb_y[k] * ld_257[k];

        t_429[k] = f_6 * if0_279[k]
                   - f_7 * if1_279[k]
                   + pa_y[k] * kf_349[k];

        t_430[k] = pa_y[k] * kf_350[k];

        t_431[k] = f_3 * kd_210[k]
                   + pb_y[k] * ld_258[k];

        t_432[k] = pa_y[k] * kf_352[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, pa_y, pb_x, pb_z, kd_207, kd_213, \
                         kf_356, ld_261, ld_262, ld_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = pb_x[k] * ld_261[k];

        t_434[k] = pb_x[k] * ld_262[k];

        t_435[k] = pb_x[k] * ld_263[k];

        t_436[k] = f_5 * kd_213[k]
                   + pa_y[k] * kf_356[k];

        t_437[k] = f_4 * kd_207[k]
                   + pb_z[k] * ld_261[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pa_y, pb_x, pb_y, pb_z, kd_210, \
                         kd_215, kf_359, lp0_132, lp1_132, ld_263, \
                         ld_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_3 * kd_215[k]
                   + pb_y[k] * ld_263[k];

        t_439[k] = pa_y[k] * kf_359[k];

        t_440[k] = f_1 * lp0_132[k]
                   - f_2 * lp1_132[k]
                   + pb_x[k] * ld_264[k];

        t_441[k] = pb_y[k] * ld_264[k];

        t_442[k] = f_0 * kd_210[k]
                   + pb_z[k] * ld_264[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, pb_x, pb_y, pb_z, kd_213, \
                         lp0_133, lp1_133, ld_267, ld_268, ld_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pb_x[k] * ld_267[k];

        t_444[k] = pb_x[k] * ld_268[k];

        t_445[k] = pb_x[k] * ld_269[k];

        t_446[k] = f_1 * lp0_133[k]
                   - f_2 * lp1_133[k]
                   + pb_y[k] * ld_267[k];

        t_447[k] = f_0 * kd_213[k]
                   + pb_z[k] * ld_267[k];

        t_448[k] = pb_y[k] * ld_269[k];
    }

#pragma omp simd aligned(t_449, pb_z, kd_215, lp0_134, lp1_134, \
                         ld_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * kd_215[k]
                   + f_1 * lp0_134[k]
                   - f_2 * lp1_134[k]
                   + pb_z[k] * ld_269[k];
    }
}

}  // namespace simdt2ceri
