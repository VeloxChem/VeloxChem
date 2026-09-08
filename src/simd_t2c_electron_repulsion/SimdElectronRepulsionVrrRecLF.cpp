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

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_4 = buffer.data(kf + 4);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_55 = buffer.data(ld + 55);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_58 = buffer.data(ld + 58);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_62 = buffer.data(ld + 62);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_64 = buffer.data(ld + 64);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_67 = buffer.data(ld + 67);
    const auto *ld_68 = buffer.data(ld + 68);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_73 = buffer.data(ld + 73);
    const auto *ld_74 = buffer.data(ld + 74);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_79 = buffer.data(ld + 79);
    const auto *ld_80 = buffer.data(ld + 80);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_85 = buffer.data(ld + 85);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_88 = buffer.data(ld + 88);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_92 = buffer.data(ld + 92);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_94 = buffer.data(ld + 94);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_97 = buffer.data(ld + 97);
    const auto *ld_98 = buffer.data(ld + 98);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_103 = buffer.data(ld + 103);
    const auto *ld_104 = buffer.data(ld + 104);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_109 = buffer.data(ld + 109);
    const auto *ld_110 = buffer.data(ld + 110);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_115 = buffer.data(ld + 115);
    const auto *ld_116 = buffer.data(ld + 116);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_121 = buffer.data(ld + 121);
    const auto *ld_122 = buffer.data(ld + 122);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_124 = buffer.data(ld + 124);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_127 = buffer.data(ld + 127);
    const auto *ld_128 = buffer.data(ld + 128);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_130 = buffer.data(ld + 130);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_133 = buffer.data(ld + 133);
    const auto *ld_134 = buffer.data(ld + 134);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_139 = buffer.data(ld + 139);
    const auto *ld_140 = buffer.data(ld + 140);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_145 = buffer.data(ld + 145);
    const auto *ld_146 = buffer.data(ld + 146);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_151 = buffer.data(ld + 151);
    const auto *ld_152 = buffer.data(ld + 152);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_157 = buffer.data(ld + 157);
    const auto *ld_158 = buffer.data(ld + 158);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_163 = buffer.data(ld + 163);
    const auto *ld_164 = buffer.data(ld + 164);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_166 = buffer.data(ld + 166);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_169 = buffer.data(ld + 169);
    const auto *ld_170 = buffer.data(ld + 170);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_172 = buffer.data(ld + 172);
    const auto *ld_173 = buffer.data(ld + 173);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_175 = buffer.data(ld + 175);
    const auto *ld_176 = buffer.data(ld + 176);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, kd_1, lp0_0, lp1_0, \
                         ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_0 * kd_1[k]
                 + pb_x[k] * ld_2[k];

        t_4[k] = pb_y[k] * ld_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, kd_2, lp0_1, lp0_2, lp1_1, \
                         lp1_2, ld_2, ld_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * kd_2[k]
                 + pb_x[k] * ld_3[k];

        t_6[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_2[k];

        t_7[k] = pb_z[k] * ld_2[k];

        t_8[k] = pb_y[k] * ld_3[k];

        t_9[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, kd_0, kd_4, \
                         kf_0, ld_4, ld_5, ld_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * kf_0[k];

        t_11[k] = f_3 * kd_0[k]
                  + pb_y[k] * ld_4[k];

        t_12[k] = pb_z[k] * ld_4[k];

        t_13[k] = f_4 * kd_4[k]
                  + pb_x[k] * ld_6[k];

        t_14[k] = pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, kd_1, kd_2, kf_2, \
                         kf_3, kf_4, ld_6, ld_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * kf_2[k];

        t_16[k] = f_5 * kd_1[k]
                  + pa_y[k] * kf_3[k];

        t_17[k] = pb_z[k] * ld_6[k];

        t_18[k] = f_3 * kd_2[k]
                  + pb_y[k] * ld_7[k];

        t_19[k] = pa_y[k] * kf_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, kd_0, kf_0, kf_1, \
                         ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * kf_0[k];

        t_21[k] = pb_y[k] * ld_8[k];

        t_22[k] = f_3 * kd_0[k]
                  + pb_z[k] * ld_8[k];

        t_23[k] = pa_z[k] * kf_1[k];

        t_24[k] = pb_y[k] * ld_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, kd_1, kd_2, \
                         kd_8, kf_3, kf_4, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * kd_8[k]
                  + pb_x[k] * ld_11[k];

        t_26[k] = pa_z[k] * kf_3[k];

        t_27[k] = f_3 * kd_1[k]
                  + pb_z[k] * ld_10[k];

        t_28[k] = pb_y[k] * ld_11[k];

        t_29[k] = f_5 * kd_2[k]
                  + pa_z[k] * kf_4[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, if0_0, if1_0, kd_3, \
                         kd_10, kf_5, ld_12, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * if0_0[k]
                  - f_7 * if1_0[k]
                  + pa_y[k] * kf_5[k];

        t_31[k] = f_8 * kd_3[k]
                  + pb_y[k] * ld_12[k];

        t_32[k] = pb_z[k] * ld_12[k];

        t_33[k] = f_9 * kd_10[k]
                  + pb_x[k] * ld_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, if0_4, if1_4, kd_11, kf_16, \
                         ld_13, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * ld_13[k];

        t_35[k] = f_9 * kd_11[k]
                  + pb_x[k] * ld_15[k];

        t_36[k] = f_10 * if0_4[k]
                  - f_11 * if1_4[k]
                  + pa_x[k] * kf_16[k];

        t_37[k] = pb_z[k] * ld_14[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, kd_5, kf_6, \
                         kf_9, kf_10, lp0_3, lp1_3, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * kd_5[k]
                  + pb_y[k] * ld_15[k];

        t_39[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_15[k];

        t_40[k] = pa_y[k] * kf_9[k];

        t_41[k] = pa_z[k] * kf_6[k];

        t_42[k] = pa_y[k] * kf_10[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, kd_4, kd_13, \
                         kf_7, kf_8, kf_11, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * kf_7[k];

        t_44[k] = f_9 * kd_13[k]
                  + pb_x[k] * ld_17[k];

        t_45[k] = pa_y[k] * kf_11[k];

        t_46[k] = pa_z[k] * kf_8[k];

        t_47[k] = f_3 * kd_4[k]
                  + pb_z[k] * ld_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, if0_0, if1_0, kd_8, kf_9, \
                         kf_12, ld_18, ld_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * kd_8[k]
                  + pb_y[k] * ld_18[k];

        t_49[k] = pa_y[k] * kf_12[k];

        t_50[k] = f_6 * if0_0[k]
                  - f_7 * if1_0[k]
                  + pa_z[k] * kf_9[k];

        t_51[k] = pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, kd_6, kd_16, kd_17, ld_19, \
                         ld_20, ld_21, ld_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * kd_6[k]
                  + pb_z[k] * ld_19[k];

        t_53[k] = f_9 * kd_16[k]
                  + pb_x[k] * ld_21[k];

        t_54[k] = pb_y[k] * ld_20[k];

        t_55[k] = f_9 * kd_17[k]
                  + pb_x[k] * ld_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, if0_6, if1_6, kd_7, kf_22, \
                         lp0_4, lp1_4, ld_21, ld_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_21[k];

        t_57[k] = f_8 * kd_7[k]
                  + pb_z[k] * ld_21[k];

        t_58[k] = pb_y[k] * ld_22[k];

        t_59[k] = f_10 * if0_6[k]
                  - f_11 * if1_6[k]
                  + pa_x[k] * kf_22[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, if0_1, if1_1, kd_9, \
                         kd_19, kf_13, ld_23, ld_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * if0_1[k]
                  - f_13 * if1_1[k]
                  + pa_y[k] * kf_13[k];

        t_61[k] = f_5 * kd_9[k]
                  + pb_y[k] * ld_23[k];

        t_62[k] = pb_z[k] * ld_23[k];

        t_63[k] = f_14 * kd_19[k]
                  + pb_x[k] * ld_25[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, if0_8, if1_8, kd_20, kf_26, \
                         ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * ld_24[k];

        t_65[k] = f_14 * kd_20[k]
                  + pb_x[k] * ld_26[k];

        t_66[k] = f_15 * if0_8[k]
                  - f_16 * if1_8[k]
                  + pa_x[k] * kf_26[k];

        t_67[k] = pb_z[k] * ld_25[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, kd_9, kd_11, kf_13, \
                         kf_14, lp0_5, lp1_5, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * kd_11[k]
                  + pb_y[k] * ld_26[k];

        t_69[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_26[k];

        t_70[k] = pa_z[k] * kf_13[k];

        t_71[k] = pa_z[k] * kf_14[k];

        t_72[k] = f_3 * kd_9[k]
                  + pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, kd_10, kd_23, kd_24, \
                         kf_15, kf_16, ld_28, ld_29, ld_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * kf_15[k];

        t_74[k] = f_14 * kd_23[k]
                  + pb_x[k] * ld_29[k];

        t_75[k] = f_14 * kd_24[k]
                  + pb_x[k] * ld_30[k];

        t_76[k] = pa_z[k] * kf_16[k];

        t_77[k] = f_3 * kd_10[k]
                  + pb_z[k] * ld_28[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, kd_11, kd_14, kd_15, \
                         kf_17, kf_18, kf_19, ld_30, ld_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * kd_14[k]
                  + pb_y[k] * ld_30[k];

        t_79[k] = f_5 * kd_11[k]
                  + pa_z[k] * kf_17[k];

        t_80[k] = pa_y[k] * kf_18[k];

        t_81[k] = f_3 * kd_15[k]
                  + pb_y[k] * ld_31[k];

        t_82[k] = pa_y[k] * kf_19[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, kd_12, kd_16, kd_26, \
                         kd_27, kf_20, kf_21, ld_32, ld_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_14 * kd_26[k]
                  + pb_x[k] * ld_32[k];

        t_84[k] = f_14 * kd_27[k]
                  + pb_x[k] * ld_33[k];

        t_85[k] = pa_y[k] * kf_20[k];

        t_86[k] = f_5 * kd_16[k]
                  + pa_y[k] * kf_21[k];

        t_87[k] = f_8 * kd_12[k]
                  + pb_z[k] * ld_32[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, if0_2, if1_2, kd_17, kf_18, \
                         kf_22, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * kd_17[k]
                  + pb_y[k] * ld_34[k];

        t_89[k] = pa_y[k] * kf_22[k];

        t_90[k] = f_12 * if0_2[k]
                  - f_13 * if1_2[k]
                  + pa_z[k] * kf_18[k];

        t_91[k] = pb_y[k] * ld_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, kd_15, kd_30, kd_31, ld_35, \
                         ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * kd_15[k]
                  + pb_z[k] * ld_35[k];

        t_93[k] = f_14 * kd_30[k]
                  + pb_x[k] * ld_37[k];

        t_94[k] = pb_y[k] * ld_36[k];

        t_95[k] = f_14 * kd_31[k]
                  + pb_x[k] * ld_38[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, if0_11, if1_11, kd_16, \
                         kf_33, lp0_6, lp1_6, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_37[k];

        t_97[k] = f_5 * kd_16[k]
                  + pb_z[k] * ld_37[k];

        t_98[k] = pb_y[k] * ld_38[k];

        t_99[k] = f_15 * if0_11[k]
                  - f_16 * if1_11[k]
                  + pa_x[k] * kf_33[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, pb_z, if0_3, if1_3, \
                         kd_18, kd_33, kf_23, ld_39, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * if0_3[k]
                   - f_18 * if1_3[k]
                   + pa_y[k] * kf_23[k];

        t_101[k] = f_19 * kd_18[k]
                   + pb_y[k] * ld_39[k];

        t_102[k] = pb_z[k] * ld_39[k];

        t_103[k] = f_19 * kd_33[k]
                   + pb_x[k] * ld_41[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, pb_z, if0_13, if1_13, kd_34, \
                         kf_37, ld_40, ld_41, ld_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * ld_40[k];

        t_105[k] = f_19 * kd_34[k]
                   + pb_x[k] * ld_42[k];

        t_106[k] = f_17 * if0_13[k]
                   - f_18 * if1_13[k]
                   + pa_x[k] * kf_37[k];

        t_107[k] = pb_z[k] * ld_41[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, kd_18, kd_20, \
                         kf_23, kf_24, lp0_7, lp1_7, ld_42, ld_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_19 * kd_20[k]
                   + pb_y[k] * ld_42[k];

        t_109[k] = f_1 * lp0_7[k]
                   - f_2 * lp1_7[k]
                   + pb_z[k] * ld_42[k];

        t_110[k] = pa_z[k] * kf_23[k];

        t_111[k] = pa_z[k] * kf_24[k];

        t_112[k] = f_3 * kd_18[k]
                   + pb_z[k] * ld_43[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, kd_19, kd_37, \
                         kd_38, kf_25, kf_26, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * kf_25[k];

        t_114[k] = f_19 * kd_37[k]
                   + pb_x[k] * ld_45[k];

        t_115[k] = f_19 * kd_38[k]
                   + pb_x[k] * ld_46[k];

        t_116[k] = pa_z[k] * kf_26[k];

        t_117[k] = f_3 * kd_19[k]
                   + pb_z[k] * ld_44[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pa_z, pb_y, if0_5, if1_5, kd_20, \
                         kd_24, kd_25, kf_27, kf_28, ld_46, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * kd_24[k]
                   + pb_y[k] * ld_46[k];

        t_119[k] = f_5 * kd_20[k]
                   + pa_z[k] * kf_27[k];

        t_120[k] = f_6 * if0_5[k]
                   - f_7 * if1_5[k]
                   + pa_y[k] * kf_28[k];

        t_121[k] = f_8 * kd_25[k]
                   + pb_y[k] * ld_47[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pb_x, pb_z, kd_21, kd_40, kd_41, kd_42, \
                         ld_47, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * kd_21[k]
                   + pb_z[k] * ld_47[k];

        t_123[k] = f_19 * kd_40[k]
                   + pb_x[k] * ld_48[k];

        t_124[k] = f_19 * kd_41[k]
                   + pb_x[k] * ld_49[k];

        t_125[k] = f_19 * kd_42[k]
                   + pb_x[k] * ld_50[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, if0_15, if1_15, kd_22, kd_28, \
                         kf_40, ld_48, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_17 * if0_15[k]
                   - f_18 * if1_15[k]
                   + pa_x[k] * kf_40[k];

        t_127[k] = f_8 * kd_22[k]
                   + pb_z[k] * ld_48[k];

        t_128[k] = f_8 * kd_28[k]
                   + pb_y[k] * ld_50[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, if0_16, if1_16, kd_29, \
                         kf_29, kf_30, kf_41, ld_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_17 * if0_16[k]
                   - f_18 * if1_16[k]
                   + pa_x[k] * kf_41[k];

        t_130[k] = pa_y[k] * kf_29[k];

        t_131[k] = f_3 * kd_29[k]
                   + pb_y[k] * ld_51[k];

        t_132[k] = pa_y[k] * kf_30[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pa_y, pb_x, pb_z, kd_26, kd_30, \
                         kd_44, kd_45, kf_31, kf_32, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_19 * kd_44[k]
                   + pb_x[k] * ld_52[k];

        t_134[k] = f_19 * kd_45[k]
                   + pb_x[k] * ld_53[k];

        t_135[k] = pa_y[k] * kf_31[k];

        t_136[k] = f_5 * kd_30[k]
                   + pa_y[k] * kf_32[k];

        t_137[k] = f_5 * kd_26[k]
                   + pb_z[k] * ld_52[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pa_z, pb_y, if0_5, if1_5, kd_31, \
                         kf_29, kf_33, ld_54, ld_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * kd_31[k]
                   + pb_y[k] * ld_54[k];

        t_139[k] = pa_y[k] * kf_33[k];

        t_140[k] = f_17 * if0_5[k]
                   - f_18 * if1_5[k]
                   + pa_z[k] * kf_29[k];

        t_141[k] = pb_y[k] * ld_55[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, kd_29, kd_48, kd_49, \
                         ld_55, ld_56, ld_57, ld_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_19 * kd_29[k]
                   + pb_z[k] * ld_55[k];

        t_143[k] = f_19 * kd_48[k]
                   + pb_x[k] * ld_57[k];

        t_144[k] = pb_y[k] * ld_56[k];

        t_145[k] = f_19 * kd_49[k]
                   + pb_x[k] * ld_58[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, pb_z, if0_19, if1_19, kd_30, \
                         kf_47, lp0_8, lp1_8, ld_57, ld_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * lp0_8[k]
                   - f_2 * lp1_8[k]
                   + pb_y[k] * ld_57[k];

        t_147[k] = f_19 * kd_30[k]
                   + pb_z[k] * ld_57[k];

        t_148[k] = pb_y[k] * ld_58[k];

        t_149[k] = f_17 * if0_19[k]
                   - f_18 * if1_19[k]
                   + pa_x[k] * kf_47[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pb_x, pb_y, pb_z, if0_7, if1_7, \
                         kd_32, kd_51, kf_34, ld_59, ld_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_15 * if0_7[k]
                   - f_16 * if1_7[k]
                   + pa_y[k] * kf_34[k];

        t_151[k] = f_14 * kd_32[k]
                   + pb_y[k] * ld_59[k];

        t_152[k] = pb_z[k] * ld_59[k];

        t_153[k] = f_5 * kd_51[k]
                   + pb_x[k] * ld_61[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, if0_20, if1_20, kd_52, \
                         kf_51, ld_60, ld_61, ld_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * ld_60[k];

        t_155[k] = f_5 * kd_52[k]
                   + pb_x[k] * ld_62[k];

        t_156[k] = f_12 * if0_20[k]
                   - f_13 * if1_20[k]
                   + pa_x[k] * kf_51[k];

        t_157[k] = pb_z[k] * ld_61[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_y, pb_z, kd_32, kd_34, \
                         kf_34, kf_35, lp0_9, lp1_9, ld_62, ld_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_14 * kd_34[k]
                   + pb_y[k] * ld_62[k];

        t_159[k] = f_1 * lp0_9[k]
                   - f_2 * lp1_9[k]
                   + pb_z[k] * ld_62[k];

        t_160[k] = pa_z[k] * kf_34[k];

        t_161[k] = pa_z[k] * kf_35[k];

        t_162[k] = f_3 * kd_32[k]
                   + pb_z[k] * ld_63[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_z, pb_x, pb_z, kd_33, kd_55, \
                         kd_56, kf_36, kf_37, ld_64, ld_65, ld_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * kf_36[k];

        t_164[k] = f_5 * kd_55[k]
                   + pb_x[k] * ld_65[k];

        t_165[k] = f_5 * kd_56[k]
                   + pb_x[k] * ld_66[k];

        t_166[k] = pa_z[k] * kf_37[k];

        t_167[k] = f_3 * kd_33[k]
                   + pb_z[k] * ld_64[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, if0_9, if1_9, kd_34, \
                         kd_38, kd_39, kf_38, kf_39, ld_66, ld_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * kd_38[k]
                   + pb_y[k] * ld_66[k];

        t_169[k] = f_5 * kd_34[k]
                   + pa_z[k] * kf_38[k];

        t_170[k] = f_12 * if0_9[k]
                   - f_13 * if1_9[k]
                   + pa_y[k] * kf_39[k];

        t_171[k] = f_5 * kd_39[k]
                   + pb_y[k] * ld_67[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, pb_z, kd_35, kd_58, kd_59, kd_60, \
                         ld_67, ld_68, ld_69, ld_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_8 * kd_35[k]
                   + pb_z[k] * ld_67[k];

        t_173[k] = f_5 * kd_58[k]
                   + pb_x[k] * ld_68[k];

        t_174[k] = f_5 * kd_59[k]
                   + pb_x[k] * ld_69[k];

        t_175[k] = f_5 * kd_60[k]
                   + pb_x[k] * ld_70[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, pb_y, pb_z, if0_21, if1_21, kd_36, kd_42, \
                         kf_54, ld_68, ld_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_12 * if0_21[k]
                   - f_13 * if1_21[k]
                   + pa_x[k] * kf_54[k];

        t_177[k] = f_8 * kd_36[k]
                   + pb_z[k] * ld_68[k];

        t_178[k] = f_5 * kd_42[k]
                   + pb_y[k] * ld_70[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pb_y, if0_10, if0_22, if1_10, \
                         if1_22, kd_43, kf_42, kf_55, ld_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_12 * if0_22[k]
                   - f_13 * if1_22[k]
                   + pa_x[k] * kf_55[k];

        t_180[k] = f_6 * if0_10[k]
                   - f_7 * if1_10[k]
                   + pa_y[k] * kf_42[k];

        t_181[k] = f_8 * kd_43[k]
                   + pb_y[k] * ld_71[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_x, pb_z, kd_39, kd_62, kd_63, kd_64, \
                         ld_71, ld_72, ld_73, ld_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_5 * kd_39[k]
                   + pb_z[k] * ld_71[k];

        t_183[k] = f_5 * kd_62[k]
                   + pb_x[k] * ld_72[k];

        t_184[k] = f_5 * kd_63[k]
                   + pb_x[k] * ld_73[k];

        t_185[k] = f_5 * kd_64[k]
                   + pb_x[k] * ld_74[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_x, pb_y, pb_z, if0_23, if1_23, kd_40, kd_46, \
                         kf_57, ld_72, ld_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_12 * if0_23[k]
                   - f_13 * if1_23[k]
                   + pa_x[k] * kf_57[k];

        t_187[k] = f_5 * kd_40[k]
                   + pb_z[k] * ld_72[k];

        t_188[k] = f_8 * kd_46[k]
                   + pb_y[k] * ld_74[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pa_y, pb_y, if0_24, if1_24, kd_47, \
                         kf_43, kf_44, kf_58, ld_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * if0_24[k]
                   - f_13 * if1_24[k]
                   + pa_x[k] * kf_58[k];

        t_190[k] = pa_y[k] * kf_43[k];

        t_191[k] = f_3 * kd_47[k]
                   + pb_y[k] * ld_75[k];

        t_192[k] = pa_y[k] * kf_44[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_z, kd_44, kd_48, \
                         kd_66, kd_67, kf_45, kf_46, ld_76, ld_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_5 * kd_66[k]
                   + pb_x[k] * ld_76[k];

        t_194[k] = f_5 * kd_67[k]
                   + pb_x[k] * ld_77[k];

        t_195[k] = pa_y[k] * kf_45[k];

        t_196[k] = f_5 * kd_48[k]
                   + pa_y[k] * kf_46[k];

        t_197[k] = f_19 * kd_44[k]
                   + pb_z[k] * ld_76[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_y, pa_z, pb_y, if0_10, if1_10, kd_49, \
                         kf_43, kf_47, ld_78, ld_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * kd_49[k]
                   + pb_y[k] * ld_78[k];

        t_199[k] = pa_y[k] * kf_47[k];

        t_200[k] = f_15 * if0_10[k]
                   - f_16 * if1_10[k]
                   + pa_z[k] * kf_43[k];

        t_201[k] = pb_y[k] * ld_79[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, kd_47, kd_70, kd_71, \
                         ld_79, ld_80, ld_81, ld_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_14 * kd_47[k]
                   + pb_z[k] * ld_79[k];

        t_203[k] = f_5 * kd_70[k]
                   + pb_x[k] * ld_81[k];

        t_204[k] = pb_y[k] * ld_80[k];

        t_205[k] = f_5 * kd_71[k]
                   + pb_x[k] * ld_82[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pb_y, pb_z, if0_25, if1_25, kd_48, \
                         kf_64, lp0_10, lp1_10, ld_81, ld_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * lp0_10[k]
                   - f_2 * lp1_10[k]
                   + pb_y[k] * ld_81[k];

        t_207[k] = f_14 * kd_48[k]
                   + pb_z[k] * ld_81[k];

        t_208[k] = pb_y[k] * ld_82[k];

        t_209[k] = f_12 * if0_25[k]
                   - f_13 * if1_25[k]
                   + pa_x[k] * kf_64[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_y, pb_x, pb_y, pb_z, if0_12, if1_12, \
                         kd_50, kd_73, kf_48, ld_83, ld_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * if0_12[k]
                   - f_11 * if1_12[k]
                   + pa_y[k] * kf_48[k];

        t_211[k] = f_9 * kd_50[k]
                   + pb_y[k] * ld_83[k];

        t_212[k] = pb_z[k] * ld_83[k];

        t_213[k] = f_8 * kd_73[k]
                   + pb_x[k] * ld_85[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_x, pb_x, pb_z, if0_26, if1_26, kd_74, \
                         kf_68, ld_84, ld_85, ld_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_z[k] * ld_84[k];

        t_215[k] = f_8 * kd_74[k]
                   + pb_x[k] * ld_86[k];

        t_216[k] = f_6 * if0_26[k]
                   - f_7 * if1_26[k]
                   + pa_x[k] * kf_68[k];

        t_217[k] = pb_z[k] * ld_85[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, pa_z, pb_y, pb_z, kd_50, kd_52, \
                         kf_48, kf_49, lp0_11, lp1_11, ld_86, ld_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_9 * kd_52[k]
                   + pb_y[k] * ld_86[k];

        t_219[k] = f_1 * lp0_11[k]
                   - f_2 * lp1_11[k]
                   + pb_z[k] * ld_86[k];

        t_220[k] = pa_z[k] * kf_48[k];

        t_221[k] = pa_z[k] * kf_49[k];

        t_222[k] = f_3 * kd_50[k]
                   + pb_z[k] * ld_87[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, kd_51, kd_76, \
                         kd_77, kf_50, kf_51, ld_88, ld_89, ld_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pa_z[k] * kf_50[k];

        t_224[k] = f_8 * kd_76[k]
                   + pb_x[k] * ld_89[k];

        t_225[k] = f_8 * kd_77[k]
                   + pb_x[k] * ld_90[k];

        t_226[k] = pa_z[k] * kf_51[k];

        t_227[k] = f_3 * kd_51[k]
                   + pb_z[k] * ld_88[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, pa_z, pb_y, if0_14, if1_14, kd_52, \
                         kd_56, kd_57, kf_52, kf_53, ld_90, ld_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_14 * kd_56[k]
                   + pb_y[k] * ld_90[k];

        t_229[k] = f_5 * kd_52[k]
                   + pa_z[k] * kf_52[k];

        t_230[k] = f_17 * if0_14[k]
                   - f_18 * if1_14[k]
                   + pa_y[k] * kf_53[k];

        t_231[k] = f_19 * kd_57[k]
                   + pb_y[k] * ld_91[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pb_x, pb_z, kd_53, kd_79, kd_80, kd_81, \
                         ld_91, ld_92, ld_93, ld_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * kd_53[k]
                   + pb_z[k] * ld_91[k];

        t_233[k] = f_8 * kd_79[k]
                   + pb_x[k] * ld_92[k];

        t_234[k] = f_8 * kd_80[k]
                   + pb_x[k] * ld_93[k];

        t_235[k] = f_8 * kd_81[k]
                   + pb_x[k] * ld_94[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_x, pb_y, pb_z, if0_28, if1_28, kd_54, kd_60, \
                         kf_69, ld_92, ld_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_6 * if0_28[k]
                   - f_7 * if1_28[k]
                   + pa_x[k] * kf_69[k];

        t_237[k] = f_8 * kd_54[k]
                   + pb_z[k] * ld_92[k];

        t_238[k] = f_19 * kd_60[k]
                   + pb_y[k] * ld_94[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pa_x, pa_y, pb_y, if0_17, if0_29, if1_17, \
                         if1_29, kd_61, kf_56, kf_70, ld_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_6 * if0_29[k]
                   - f_7 * if1_29[k]
                   + pa_x[k] * kf_70[k];

        t_240[k] = f_12 * if0_17[k]
                   - f_13 * if1_17[k]
                   + pa_y[k] * kf_56[k];

        t_241[k] = f_5 * kd_61[k]
                   + pb_y[k] * ld_95[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pb_x, pb_z, kd_57, kd_83, kd_84, kd_85, \
                         ld_95, ld_96, ld_97, ld_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_5 * kd_57[k]
                   + pb_z[k] * ld_95[k];

        t_243[k] = f_8 * kd_83[k]
                   + pb_x[k] * ld_96[k];

        t_244[k] = f_8 * kd_84[k]
                   + pb_x[k] * ld_97[k];

        t_245[k] = f_8 * kd_85[k]
                   + pb_x[k] * ld_98[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pb_y, pb_z, if0_30, if1_30, kd_58, kd_64, \
                         kf_71, ld_96, ld_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_6 * if0_30[k]
                   - f_7 * if1_30[k]
                   + pa_x[k] * kf_71[k];

        t_247[k] = f_5 * kd_58[k]
                   + pb_z[k] * ld_96[k];

        t_248[k] = f_5 * kd_64[k]
                   + pb_y[k] * ld_98[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pa_y, pb_y, if0_18, if0_31, if1_18, \
                         if1_31, kd_65, kf_59, kf_72, ld_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_6 * if0_31[k]
                   - f_7 * if1_31[k]
                   + pa_x[k] * kf_72[k];

        t_250[k] = f_6 * if0_18[k]
                   - f_7 * if1_18[k]
                   + pa_y[k] * kf_59[k];

        t_251[k] = f_8 * kd_65[k]
                   + pb_y[k] * ld_99[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, pb_z, kd_61, kd_87, kd_88, kd_89, \
                         ld_99, ld_100, ld_101, ld_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_19 * kd_61[k]
                   + pb_z[k] * ld_99[k];

        t_253[k] = f_8 * kd_87[k]
                   + pb_x[k] * ld_100[k];

        t_254[k] = f_8 * kd_88[k]
                   + pb_x[k] * ld_101[k];

        t_255[k] = f_8 * kd_89[k]
                   + pb_x[k] * ld_102[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pa_x, pb_y, pb_z, if0_32, if1_32, kd_62, kd_68, \
                         kf_73, ld_100, ld_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_6 * if0_32[k]
                   - f_7 * if1_32[k]
                   + pa_x[k] * kf_73[k];

        t_257[k] = f_19 * kd_62[k]
                   + pb_z[k] * ld_100[k];

        t_258[k] = f_8 * kd_68[k]
                   + pb_y[k] * ld_102[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pa_x, pa_y, pb_y, if0_33, if1_33, kd_69, \
                         kf_60, kf_61, kf_74, ld_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_6 * if0_33[k]
                   - f_7 * if1_33[k]
                   + pa_x[k] * kf_74[k];

        t_260[k] = pa_y[k] * kf_60[k];

        t_261[k] = f_3 * kd_69[k]
                   + pb_y[k] * ld_103[k];

        t_262[k] = pa_y[k] * kf_61[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, pa_y, pb_x, pb_z, kd_66, kd_70, \
                         kd_91, kd_92, kf_62, kf_63, ld_104, ld_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_8 * kd_91[k]
                   + pb_x[k] * ld_104[k];

        t_264[k] = f_8 * kd_92[k]
                   + pb_x[k] * ld_105[k];

        t_265[k] = pa_y[k] * kf_62[k];

        t_266[k] = f_5 * kd_70[k]
                   + pa_y[k] * kf_63[k];

        t_267[k] = f_14 * kd_66[k]
                   + pb_z[k] * ld_104[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pa_z, pb_y, if0_18, if1_18, kd_71, \
                         kf_60, kf_64, ld_106, ld_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_3 * kd_71[k]
                   + pb_y[k] * ld_106[k];

        t_269[k] = pa_y[k] * kf_64[k];

        t_270[k] = f_10 * if0_18[k]
                   - f_11 * if1_18[k]
                   + pa_z[k] * kf_60[k];

        t_271[k] = pb_y[k] * ld_107[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, kd_69, kd_94, kd_95, \
                         ld_107, ld_108, ld_109, ld_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_9 * kd_69[k]
                   + pb_z[k] * ld_107[k];

        t_273[k] = f_8 * kd_94[k]
                   + pb_x[k] * ld_109[k];

        t_274[k] = pb_y[k] * ld_108[k];

        t_275[k] = f_8 * kd_95[k]
                   + pb_x[k] * ld_110[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pb_y, pb_z, if0_35, if1_35, kd_70, \
                         kf_78, lp0_12, lp1_12, ld_109, ld_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * lp0_12[k]
                   - f_2 * lp1_12[k]
                   + pb_y[k] * ld_109[k];

        t_277[k] = f_9 * kd_70[k]
                   + pb_z[k] * ld_109[k];

        t_278[k] = pb_y[k] * ld_110[k];

        t_279[k] = f_6 * if0_35[k]
                   - f_7 * if1_35[k]
                   + pa_x[k] * kf_78[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pa_x, pb_x, pb_y, pb_z, kd_72, \
                         kd_96, kd_97, kf_79, ld_111, ld_112, ld_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_5 * kd_96[k]
                   + pa_x[k] * kf_79[k];

        t_281[k] = f_4 * kd_72[k]
                   + pb_y[k] * ld_111[k];

        t_282[k] = pb_z[k] * ld_111[k];

        t_283[k] = f_3 * kd_97[k]
                   + pb_x[k] * ld_113[k];

        t_284[k] = pb_z[k] * ld_112[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, pa_x, pb_x, pb_z, kd_98, kf_81, \
                         kf_82, kf_83, ld_113, ld_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_3 * kd_98[k]
                   + pb_x[k] * ld_114[k];

        t_286[k] = pa_x[k] * kf_81[k];

        t_287[k] = pb_z[k] * ld_113[k];

        t_288[k] = pa_x[k] * kf_82[k];

        t_289[k] = pa_x[k] * kf_83[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pa_z, pb_x, pb_z, kd_72, kd_101, \
                         kf_65, kf_66, kf_67, ld_115, ld_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_z[k] * kf_65[k];

        t_291[k] = pa_z[k] * kf_66[k];

        t_292[k] = f_3 * kd_72[k]
                   + pb_z[k] * ld_115[k];

        t_293[k] = pa_z[k] * kf_67[k];

        t_294[k] = f_3 * kd_101[k]
                   + pb_x[k] * ld_116[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, t_300, pa_x, pb_x, kd_102, kd_103, \
                         kf_84, kf_85, kf_86, kf_87, kf_88, ld_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_3 * kd_102[k]
                   + pb_x[k] * ld_117[k];

        t_296[k] = pa_x[k] * kf_84[k];

        t_297[k] = pa_x[k] * kf_85[k];

        t_298[k] = pa_x[k] * kf_86[k];

        t_299[k] = pa_x[k] * kf_87[k];

        t_300[k] = f_5 * kd_103[k]
                   + pa_x[k] * kf_88[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_x, pb_y, pb_z, kd_75, kd_78, kd_104, \
                         kd_105, ld_118, ld_119, ld_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * kd_78[k]
                   + pb_y[k] * ld_118[k];

        t_302[k] = f_8 * kd_75[k]
                   + pb_z[k] * ld_118[k];

        t_303[k] = f_3 * kd_104[k]
                   + pb_x[k] * ld_119[k];

        t_304[k] = f_3 * kd_105[k]
                   + pb_x[k] * ld_120[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, t_310, pa_x, pb_x, kd_106, kd_107, \
                         kf_89, kf_90, kf_91, kf_92, kf_93, ld_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_3 * kd_106[k]
                   + pb_x[k] * ld_121[k];

        t_306[k] = pa_x[k] * kf_89[k];

        t_307[k] = pa_x[k] * kf_90[k];

        t_308[k] = pa_x[k] * kf_91[k];

        t_309[k] = pa_x[k] * kf_92[k];

        t_310[k] = f_5 * kd_107[k]
                   + pa_x[k] * kf_93[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_x, pb_y, pb_z, kd_78, kd_82, kd_108, \
                         kd_109, ld_122, ld_123, ld_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_19 * kd_82[k]
                   + pb_y[k] * ld_122[k];

        t_312[k] = f_5 * kd_78[k]
                   + pb_z[k] * ld_122[k];

        t_313[k] = f_3 * kd_108[k]
                   + pb_x[k] * ld_123[k];

        t_314[k] = f_3 * kd_109[k]
                   + pb_x[k] * ld_124[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, pa_x, pb_x, kd_110, kd_111, \
                         kf_94, kf_95, kf_96, kf_97, kf_98, ld_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_3 * kd_110[k]
                   + pb_x[k] * ld_125[k];

        t_316[k] = pa_x[k] * kf_94[k];

        t_317[k] = pa_x[k] * kf_95[k];

        t_318[k] = pa_x[k] * kf_96[k];

        t_319[k] = pa_x[k] * kf_97[k];

        t_320[k] = f_5 * kd_111[k]
                   + pa_x[k] * kf_98[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_y, pb_z, kd_82, kd_86, kd_112, \
                         kd_113, ld_126, ld_127, ld_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_5 * kd_86[k]
                   + pb_y[k] * ld_126[k];

        t_322[k] = f_19 * kd_82[k]
                   + pb_z[k] * ld_126[k];

        t_323[k] = f_3 * kd_112[k]
                   + pb_x[k] * ld_127[k];

        t_324[k] = f_3 * kd_113[k]
                   + pb_x[k] * ld_128[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, t_330, pa_x, pb_x, kd_114, kd_115, \
                         kf_99, kf_100, kf_101, kf_102, kf_103, \
                         ld_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_3 * kd_114[k]
                   + pb_x[k] * ld_129[k];

        t_326[k] = pa_x[k] * kf_99[k];

        t_327[k] = pa_x[k] * kf_100[k];

        t_328[k] = pa_x[k] * kf_101[k];

        t_329[k] = pa_x[k] * kf_102[k];

        t_330[k] = f_5 * kd_115[k]
                   + pa_x[k] * kf_103[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_x, pb_y, pb_z, kd_86, kd_90, kd_116, \
                         kd_117, ld_130, ld_131, ld_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_8 * kd_90[k]
                   + pb_y[k] * ld_130[k];

        t_332[k] = f_14 * kd_86[k]
                   + pb_z[k] * ld_130[k];

        t_333[k] = f_3 * kd_116[k]
                   + pb_x[k] * ld_131[k];

        t_334[k] = f_3 * kd_117[k]
                   + pb_x[k] * ld_132[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, t_340, pa_x, pa_y, pb_x, kd_118, \
                         kf_75, kf_104, kf_105, kf_106, kf_107, \
                         ld_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_3 * kd_118[k]
                   + pb_x[k] * ld_133[k];

        t_336[k] = pa_x[k] * kf_104[k];

        t_337[k] = pa_x[k] * kf_105[k];

        t_338[k] = pa_x[k] * kf_106[k];

        t_339[k] = pa_x[k] * kf_107[k];

        t_340[k] = pa_y[k] * kf_75[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, pa_y, pb_x, pb_y, kd_93, kd_120, \
                         kd_121, kf_76, kf_77, ld_134, ld_135, ld_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_3 * kd_93[k]
                   + pb_y[k] * ld_134[k];

        t_342[k] = pa_y[k] * kf_76[k];

        t_343[k] = f_3 * kd_120[k]
                   + pb_x[k] * ld_135[k];

        t_344[k] = f_3 * kd_121[k]
                   + pb_x[k] * ld_136[k];

        t_345[k] = pa_y[k] * kf_77[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, t_351, pa_x, pb_y, kd_123, kf_108, \
                         kf_109, kf_110, kf_111, kf_112, ld_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pa_x[k] * kf_108[k];

        t_347[k] = pa_x[k] * kf_109[k];

        t_348[k] = pa_x[k] * kf_110[k];

        t_349[k] = pa_x[k] * kf_111[k];

        t_350[k] = f_5 * kd_123[k]
                   + pa_x[k] * kf_112[k];

        t_351[k] = pb_y[k] * ld_137[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_x, pb_y, pb_z, kd_93, kd_124, kd_125, \
                         ld_137, ld_138, ld_139, ld_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_4 * kd_93[k]
                   + pb_z[k] * ld_137[k];

        t_353[k] = f_3 * kd_124[k]
                   + pb_x[k] * ld_139[k];

        t_354[k] = pb_y[k] * ld_138[k];

        t_355[k] = f_3 * kd_125[k]
                   + pb_x[k] * ld_140[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, t_360, pa_x, pb_x, pb_y, kf_114, kf_115, \
                         kf_116, lp0_13, lp1_13, ld_140, ld_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * kf_114[k];

        t_357[k] = pa_x[k] * kf_115[k];

        t_358[k] = pb_y[k] * ld_140[k];

        t_359[k] = pa_x[k] * kf_116[k];

        t_360[k] = f_1 * lp0_13[k]
                   - f_2 * lp1_13[k]
                   + pb_x[k] * ld_141[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, t_365, pb_x, pb_y, pb_z, kd_96, ld_141, \
                         ld_142, ld_143, ld_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_0 * kd_96[k]
                   + pb_y[k] * ld_141[k];

        t_362[k] = pb_z[k] * ld_141[k];

        t_363[k] = pb_x[k] * ld_142[k];

        t_364[k] = pb_x[k] * ld_143[k];

        t_365[k] = pb_x[k] * ld_144[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pb_z, kd_97, kd_98, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_142, ld_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_0 * kd_97[k]
                   + f_1 * lp0_14[k]
                   - f_2 * lp1_14[k]
                   + pb_y[k] * ld_142[k];

        t_367[k] = pb_z[k] * ld_142[k];

        t_368[k] = f_0 * kd_98[k]
                   + pb_y[k] * ld_144[k];

        t_369[k] = f_1 * lp0_15[k]
                   - f_2 * lp1_15[k]
                   + pb_z[k] * ld_144[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, t_375, pa_z, pb_x, pb_z, kd_96, \
                         kf_79, kf_80, ld_145, ld_146, ld_147, ld_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_z[k] * kf_79[k];

        t_371[k] = pa_z[k] * kf_80[k];

        t_372[k] = f_3 * kd_96[k]
                   + pb_z[k] * ld_145[k];

        t_373[k] = pb_x[k] * ld_146[k];

        t_374[k] = pb_x[k] * ld_147[k];

        t_375[k] = pb_x[k] * ld_148[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pb_y, pb_z, kd_97, kd_98, kd_102, \
                         kf_81, kf_83, ld_146, ld_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_z[k] * kf_81[k];

        t_377[k] = f_3 * kd_97[k]
                   + pb_z[k] * ld_146[k];

        t_378[k] = f_4 * kd_102[k]
                   + pb_y[k] * ld_148[k];

        t_379[k] = f_5 * kd_98[k]
                   + pa_z[k] * kf_83[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, pb_x, pb_y, pb_z, kd_99, kd_103, \
                         lp0_16, lp1_16, ld_149, ld_150, ld_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_1 * lp0_16[k]
                   - f_2 * lp1_16[k]
                   + pb_x[k] * ld_149[k];

        t_381[k] = f_9 * kd_103[k]
                   + pb_y[k] * ld_149[k];

        t_382[k] = f_8 * kd_99[k]
                   + pb_z[k] * ld_149[k];

        t_383[k] = pb_x[k] * ld_150[k];

        t_384[k] = pb_x[k] * ld_151[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_z, pb_x, pb_y, pb_z, if0_26, if1_26, \
                         kd_100, kd_106, kf_84, ld_150, ld_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_x[k] * ld_152[k];

        t_386[k] = f_6 * if0_26[k]
                   - f_7 * if1_26[k]
                   + pa_z[k] * kf_84[k];

        t_387[k] = f_8 * kd_100[k]
                   + pb_z[k] * ld_150[k];

        t_388[k] = f_9 * kd_106[k]
                   + pb_y[k] * ld_152[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_y, pb_x, pb_y, pb_z, if0_29, if1_29, \
                         kd_103, kd_107, kf_92, lp0_17, lp1_17, \
                         ld_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_10 * if0_29[k]
                   - f_11 * if1_29[k]
                   + pa_y[k] * kf_92[k];

        t_390[k] = f_1 * lp0_17[k]
                   - f_2 * lp1_17[k]
                   + pb_x[k] * ld_153[k];

        t_391[k] = f_14 * kd_107[k]
                   + pb_y[k] * ld_153[k];

        t_392[k] = f_5 * kd_103[k]
                   + pb_z[k] * ld_153[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, pa_z, pb_x, pb_z, if0_27, if1_27, \
                         kd_104, kf_89, ld_154, ld_155, ld_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pb_x[k] * ld_154[k];

        t_394[k] = pb_x[k] * ld_155[k];

        t_395[k] = pb_x[k] * ld_156[k];

        t_396[k] = f_12 * if0_27[k]
                   - f_13 * if1_27[k]
                   + pa_z[k] * kf_89[k];

        t_397[k] = f_5 * kd_104[k]
                   + pb_z[k] * ld_154[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, pa_y, pb_x, pb_y, if0_31, if1_31, kd_110, \
                         kd_111, kf_97, lp0_18, lp1_18, ld_156, \
                         ld_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_14 * kd_110[k]
                   + pb_y[k] * ld_156[k];

        t_399[k] = f_15 * if0_31[k]
                   - f_16 * if1_31[k]
                   + pa_y[k] * kf_97[k];

        t_400[k] = f_1 * lp0_18[k]
                   - f_2 * lp1_18[k]
                   + pb_x[k] * ld_157[k];

        t_401[k] = f_19 * kd_111[k]
                   + pb_y[k] * ld_157[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_z, pb_x, pb_z, if0_28, if1_28, \
                         kd_107, kf_94, ld_157, ld_158, ld_159, \
                         ld_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_19 * kd_107[k]
                   + pb_z[k] * ld_157[k];

        t_403[k] = pb_x[k] * ld_158[k];

        t_404[k] = pb_x[k] * ld_159[k];

        t_405[k] = pb_x[k] * ld_160[k];

        t_406[k] = f_17 * if0_28[k]
                   - f_18 * if1_28[k]
                   + pa_z[k] * kf_94[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pa_y, pb_y, pb_z, if0_33, if1_33, kd_108, \
                         kd_114, kf_102, ld_158, ld_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_19 * kd_108[k]
                   + pb_z[k] * ld_158[k];

        t_408[k] = f_19 * kd_114[k]
                   + pb_y[k] * ld_160[k];

        t_409[k] = f_17 * if0_33[k]
                   - f_18 * if1_33[k]
                   + pa_y[k] * kf_102[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pb_x, pb_y, pb_z, kd_111, kd_115, \
                         lp0_19, lp1_19, ld_161, ld_162, ld_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_1 * lp0_19[k]
                   - f_2 * lp1_19[k]
                   + pb_x[k] * ld_161[k];

        t_411[k] = f_5 * kd_115[k]
                   + pb_y[k] * ld_161[k];

        t_412[k] = f_14 * kd_111[k]
                   + pb_z[k] * ld_161[k];

        t_413[k] = pb_x[k] * ld_162[k];

        t_414[k] = pb_x[k] * ld_163[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_z, pb_x, pb_y, pb_z, if0_30, if1_30, \
                         kd_112, kd_118, kf_99, ld_162, ld_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = pb_x[k] * ld_164[k];

        t_416[k] = f_15 * if0_30[k]
                   - f_16 * if1_30[k]
                   + pa_z[k] * kf_99[k];

        t_417[k] = f_14 * kd_112[k]
                   + pb_z[k] * ld_162[k];

        t_418[k] = f_5 * kd_118[k]
                   + pb_y[k] * ld_164[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pa_y, pb_x, pb_y, pb_z, if0_34, if1_34, \
                         kd_115, kd_119, kf_107, lp0_20, lp1_20, \
                         ld_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_12 * if0_34[k]
                   - f_13 * if1_34[k]
                   + pa_y[k] * kf_107[k];

        t_420[k] = f_1 * lp0_20[k]
                   - f_2 * lp1_20[k]
                   + pb_x[k] * ld_165[k];

        t_421[k] = f_8 * kd_119[k]
                   + pb_y[k] * ld_165[k];

        t_422[k] = f_9 * kd_115[k]
                   + pb_z[k] * ld_165[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, pa_z, pb_x, pb_z, if0_32, if1_32, \
                         kd_116, kf_104, ld_166, ld_167, ld_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = pb_x[k] * ld_166[k];

        t_424[k] = pb_x[k] * ld_167[k];

        t_425[k] = pb_x[k] * ld_168[k];

        t_426[k] = f_10 * if0_32[k]
                   - f_11 * if1_32[k]
                   + pa_z[k] * kf_104[k];

        t_427[k] = f_9 * kd_116[k]
                   + pb_z[k] * ld_166[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pa_y, pb_y, if0_35, if1_35, \
                         kd_122, kd_123, kf_111, kf_112, kf_113, ld_168, \
                         ld_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_8 * kd_122[k]
                   + pb_y[k] * ld_168[k];

        t_429[k] = f_6 * if0_35[k]
                   - f_7 * if1_35[k]
                   + pa_y[k] * kf_111[k];

        t_430[k] = pa_y[k] * kf_112[k];

        t_431[k] = f_3 * kd_123[k]
                   + pb_y[k] * ld_169[k];

        t_432[k] = pa_y[k] * kf_113[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, pa_y, pb_x, pb_z, kd_120, kd_124, \
                         kf_114, ld_170, ld_171, ld_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = pb_x[k] * ld_170[k];

        t_434[k] = pb_x[k] * ld_171[k];

        t_435[k] = pb_x[k] * ld_172[k];

        t_436[k] = f_5 * kd_124[k]
                   + pa_y[k] * kf_114[k];

        t_437[k] = f_4 * kd_120[k]
                   + pb_z[k] * ld_170[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pa_y, pb_x, pb_y, pb_z, kd_123, \
                         kd_125, kf_116, lp0_21, lp1_21, ld_172, \
                         ld_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_3 * kd_125[k]
                   + pb_y[k] * ld_172[k];

        t_439[k] = pa_y[k] * kf_116[k];

        t_440[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_173[k];

        t_441[k] = pb_y[k] * ld_173[k];

        t_442[k] = f_0 * kd_123[k]
                   + pb_z[k] * ld_173[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, pb_x, pb_y, pb_z, kd_124, \
                         lp0_22, lp1_22, ld_174, ld_175, ld_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pb_x[k] * ld_174[k];

        t_444[k] = pb_x[k] * ld_175[k];

        t_445[k] = pb_x[k] * ld_176[k];

        t_446[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_174[k];

        t_447[k] = f_0 * kd_124[k]
                   + pb_z[k] * ld_174[k];

        t_448[k] = pb_y[k] * ld_176[k];
    }

#pragma omp simd aligned(t_449, pb_z, kd_125, lp0_23, lp1_23, ld_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * kd_125[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_176[k];
    }
}

auto
compute_prim_lf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_1 = buffer.data(if0 + 1);
    const auto *if0_2 = buffer.data(if0 + 2);
    const auto *if0_3 = buffer.data(if0 + 3);
    const auto *if0_5 = buffer.data(if0 + 5);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_9 = buffer.data(if0 + 9);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_12 = buffer.data(if0 + 12);
    const auto *if0_13 = buffer.data(if0 + 13);
    const auto *if0_15 = buffer.data(if0 + 15);
    const auto *if0_16 = buffer.data(if0 + 16);
    const auto *if0_18 = buffer.data(if0 + 18);
    const auto *if0_19 = buffer.data(if0 + 19);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_23 = buffer.data(if0 + 23);
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
    const auto *if0_36 = buffer.data(if0 + 36);
    const auto *if0_37 = buffer.data(if0 + 37);
    const auto *if0_39 = buffer.data(if0 + 39);
    const auto *if0_40 = buffer.data(if0 + 40);
    const auto *if0_42 = buffer.data(if0 + 42);
    const auto *if0_43 = buffer.data(if0 + 43);
    const auto *if0_44 = buffer.data(if0 + 44);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_13 = buffer.data(if1 + 13);
    const auto *if1_15 = buffer.data(if1 + 15);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_24 = buffer.data(if1 + 24);
    const auto *if1_25 = buffer.data(if1 + 25);
    const auto *if1_29 = buffer.data(if1 + 29);
    const auto *if1_30 = buffer.data(if1 + 30);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_34 = buffer.data(if1 + 34);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_36 = buffer.data(if1 + 36);
    const auto *if1_37 = buffer.data(if1 + 37);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_42 = buffer.data(if1 + 42);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_46 = buffer.data(if1 + 46);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_48 = buffer.data(if1 + 48);
    const auto *if1_49 = buffer.data(if1 + 49);
    const auto *if1_53 = buffer.data(if1 + 53);
    const auto *if1_57 = buffer.data(if1 + 57);
    const auto *if1_61 = buffer.data(if1 + 61);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_69 = buffer.data(if1 + 69);
    const auto *if1_71 = buffer.data(if1 + 71);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_76 = buffer.data(if1 + 76);
    const auto *if1_79 = buffer.data(if1 + 79);
    const auto *if1_83 = buffer.data(if1 + 83);
    const auto *if1_91 = buffer.data(if1 + 91);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_136 = buffer.data(kf + 136);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_55 = buffer.data(ld + 55);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_58 = buffer.data(ld + 58);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_62 = buffer.data(ld + 62);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_64 = buffer.data(ld + 64);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_67 = buffer.data(ld + 67);
    const auto *ld_68 = buffer.data(ld + 68);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_73 = buffer.data(ld + 73);
    const auto *ld_74 = buffer.data(ld + 74);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_79 = buffer.data(ld + 79);
    const auto *ld_80 = buffer.data(ld + 80);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_85 = buffer.data(ld + 85);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_88 = buffer.data(ld + 88);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_92 = buffer.data(ld + 92);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_94 = buffer.data(ld + 94);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_97 = buffer.data(ld + 97);
    const auto *ld_98 = buffer.data(ld + 98);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_103 = buffer.data(ld + 103);
    const auto *ld_104 = buffer.data(ld + 104);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_109 = buffer.data(ld + 109);
    const auto *ld_110 = buffer.data(ld + 110);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, kd_1, kd_2, lp0_0, \
                         lp1_0, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_0 * kd_1[k]
                 + pb_x[k] * ld_1[k];

        t_4[k] = f_0 * kd_2[k]
                 + pb_x[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, kf_0, lp0_1, lp0_2, lp1_1, \
                         lp1_2, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_6[k] = pb_y[k] * ld_2[k];

        t_7[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_8[k] = pa_y[k] * kf_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, kd_0, kd_1, kd_2, kd_4, \
                         kf_3, ld_3, ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * kd_0[k]
                 + pb_y[k] * ld_3[k];

        t_10[k] = f_4 * kd_4[k]
                  + pb_x[k] * ld_4[k];

        t_11[k] = f_5 * kd_1[k]
                  + pa_y[k] * kf_3[k];

        t_12[k] = f_3 * kd_2[k]
                  + pb_y[k] * ld_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_z, kd_0, kd_8, \
                         kf_0, kf_3, kf_5, ld_6, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * kf_5[k];

        t_14[k] = pa_z[k] * kf_0[k];

        t_15[k] = f_3 * kd_0[k]
                  + pb_z[k] * ld_6[k];

        t_16[k] = f_4 * kd_8[k]
                  + pb_x[k] * ld_8[k];

        t_17[k] = pa_z[k] * kf_3[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_z, if0_0, if1_0, kd_1, kd_2, kf_5, \
                         kf_6, ld_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * kd_1[k]
                  + pb_z[k] * ld_7[k];

        t_19[k] = f_5 * kd_2[k]
                  + pa_z[k] * kf_5[k];

        t_20[k] = f_6 * if0_0[k]
                  - f_7 * if1_0[k]
                  + pa_y[k] * kf_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pb_x, pb_y, pb_z, if0_5, if1_13, \
                         kd_3, kd_10, kf_14, ld_9, ld_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * kd_3[k]
                  + pb_y[k] * ld_9[k];

        t_22[k] = pb_z[k] * ld_9[k];

        t_23[k] = f_9 * kd_10[k]
                  + pb_x[k] * ld_10[k];

        t_24[k] = f_10 * if0_5[k]
                  - f_11 * if1_13[k]
                  + pa_x[k] * kf_14[k];

        t_25[k] = pb_z[k] * ld_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, pb_z, kd_5, kf_7, kf_9, \
                         lp0_3, lp1_3, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * kd_5[k]
                  + pb_y[k] * ld_11[k];

        t_27[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_11[k];

        t_28[k] = pa_y[k] * kf_9[k];

        t_29[k] = pa_z[k] * kf_7[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_y, pb_z, if0_0, if1_0, kd_4, \
                         kd_8, kf_8, kf_10, ld_12, ld_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * kd_4[k]
                  + pb_z[k] * ld_12[k];

        t_31[k] = f_3 * kd_8[k]
                  + pb_y[k] * ld_13[k];

        t_32[k] = pa_y[k] * kf_10[k];

        t_33[k] = f_6 * if0_0[k]
                  - f_7 * if1_0[k]
                  + pa_z[k] * kf_8[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, kd_6, kd_7, kd_16, \
                         lp0_4, lp1_4, ld_14, ld_15, ld_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * ld_14[k];

        t_35[k] = f_8 * kd_6[k]
                  + pb_z[k] * ld_14[k];

        t_36[k] = f_9 * kd_16[k]
                  + pb_x[k] * ld_16[k];

        t_37[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_15[k];

        t_38[k] = f_8 * kd_7[k]
                  + pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pa_y, pb_y, if0_1, if0_8, if1_6, \
                         if1_19, kd_9, kf_11, kf_23, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_y[k] * ld_16[k];

        t_40[k] = f_10 * if0_8[k]
                  - f_11 * if1_19[k]
                  + pa_x[k] * kf_23[k];

        t_41[k] = f_12 * if0_1[k]
                  - f_13 * if1_6[k]
                  + pa_y[k] * kf_11[k];

        t_42[k] = f_5 * kd_9[k]
                  + pb_y[k] * ld_17[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_z, if0_11, if1_22, kd_18, \
                         kf_27, ld_17, ld_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_z[k] * ld_17[k];

        t_44[k] = f_14 * kd_18[k]
                  + pb_x[k] * ld_18[k];

        t_45[k] = f_15 * if0_11[k]
                  - f_16 * if1_22[k]
                  + pa_x[k] * kf_27[k];

        t_46[k] = pb_z[k] * ld_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, kd_9, kd_11, kf_11, \
                         kf_14, lp0_5, lp1_5, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * kd_11[k]
                  + pb_y[k] * ld_19[k];

        t_48[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_19[k];

        t_49[k] = pa_z[k] * kf_11[k];

        t_50[k] = f_3 * kd_9[k]
                  + pb_z[k] * ld_20[k];

        t_51[k] = pa_z[k] * kf_14[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pa_z, pb_y, pb_z, kd_10, kd_11, kd_13, \
                         kf_16, kf_17, ld_21, ld_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * kd_10[k]
                  + pb_z[k] * ld_21[k];

        t_53[k] = f_8 * kd_13[k]
                  + pb_y[k] * ld_22[k];

        t_54[k] = f_5 * kd_11[k]
                  + pa_z[k] * kf_16[k];

        t_55[k] = pa_y[k] * kf_17[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_y, pb_z, kd_12, kd_15, kd_16, \
                         kf_19, kf_21, kf_23, ld_23, ld_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * kf_19[k];

        t_57[k] = f_5 * kd_15[k]
                  + pa_y[k] * kf_21[k];

        t_58[k] = f_8 * kd_12[k]
                  + pb_z[k] * ld_23[k];

        t_59[k] = f_3 * kd_16[k]
                  + pb_y[k] * ld_24[k];

        t_60[k] = pa_y[k] * kf_23[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_z, pb_x, pb_y, pb_z, if0_2, if1_8, kd_14, \
                         kd_27, kf_17, ld_25, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_12 * if0_2[k]
                  - f_13 * if1_8[k]
                  + pa_z[k] * kf_17[k];

        t_62[k] = pb_y[k] * ld_25[k];

        t_63[k] = f_5 * kd_14[k]
                  + pb_z[k] * ld_25[k];

        t_64[k] = f_14 * kd_27[k]
                  + pb_x[k] * ld_27[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_y, pb_z, if0_15, if1_29, kd_15, \
                         kf_37, lp0_6, lp1_6, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_26[k];

        t_66[k] = f_5 * kd_15[k]
                  + pb_z[k] * ld_26[k];

        t_67[k] = pb_y[k] * ld_27[k];

        t_68[k] = f_15 * if0_15[k]
                  - f_16 * if1_29[k]
                  + pa_x[k] * kf_37[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_x, pb_y, pb_z, if0_3, if1_11, kd_17, \
                         kd_29, kf_24, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_17 * if0_3[k]
                  - f_18 * if1_11[k]
                  + pa_y[k] * kf_24[k];

        t_70[k] = f_19 * kd_17[k]
                  + pb_y[k] * ld_28[k];

        t_71[k] = pb_z[k] * ld_28[k];

        t_72[k] = f_19 * kd_29[k]
                  + pb_x[k] * ld_29[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_y, pb_z, if0_18, if1_32, kd_19, \
                         kf_41, lp0_7, lp1_7, ld_29, ld_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_17 * if0_18[k]
                  - f_18 * if1_32[k]
                  + pa_x[k] * kf_41[k];

        t_74[k] = pb_z[k] * ld_29[k];

        t_75[k] = f_19 * kd_19[k]
                  + pb_y[k] * ld_30[k];

        t_76[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_30[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_z, pb_y, pb_z, kd_17, kd_18, kd_22, \
                         kf_24, kf_27, ld_31, ld_32, ld_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * kf_24[k];

        t_78[k] = f_3 * kd_17[k]
                  + pb_z[k] * ld_31[k];

        t_79[k] = pa_z[k] * kf_27[k];

        t_80[k] = f_3 * kd_18[k]
                  + pb_z[k] * ld_32[k];

        t_81[k] = f_5 * kd_22[k]
                  + pb_y[k] * ld_33[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pa_z, pb_z, if0_6, if1_15, kd_19, kd_20, \
                         kf_29, kf_30, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_5 * kd_19[k]
                  + pa_z[k] * kf_29[k];

        t_83[k] = f_6 * if0_6[k]
                  - f_7 * if1_15[k]
                  + pa_y[k] * kf_30[k];

        t_84[k] = f_8 * kd_20[k]
                  + pb_z[k] * ld_34[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pa_x, pb_y, pb_z, if0_20, if1_35, kd_21, kd_24, \
                         kf_45, ld_35, ld_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_17 * if0_20[k]
                  - f_18 * if1_35[k]
                  + pa_x[k] * kf_45[k];

        t_86[k] = f_8 * kd_21[k]
                  + pb_z[k] * ld_35[k];

        t_87[k] = f_8 * kd_24[k]
                  + pb_y[k] * ld_36[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, if0_21, if1_36, kd_26, kf_31, \
                         kf_33, kf_35, kf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_17 * if0_21[k]
                  - f_18 * if1_36[k]
                  + pa_x[k] * kf_46[k];

        t_89[k] = pa_y[k] * kf_31[k];

        t_90[k] = pa_y[k] * kf_33[k];

        t_91[k] = f_5 * kd_26[k]
                  + pa_y[k] * kf_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_y, pb_z, if0_6, if1_15, kd_23, \
                         kd_27, kf_31, kf_37, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * kd_23[k]
                  + pb_z[k] * ld_37[k];

        t_93[k] = f_3 * kd_27[k]
                  + pb_y[k] * ld_38[k];

        t_94[k] = pa_y[k] * kf_37[k];

        t_95[k] = f_17 * if0_6[k]
                  - f_18 * if1_15[k]
                  + pa_z[k] * kf_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, kd_25, kd_26, kd_41, \
                         lp0_8, lp1_8, ld_39, ld_40, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_y[k] * ld_39[k];

        t_97[k] = f_19 * kd_25[k]
                  + pb_z[k] * ld_39[k];

        t_98[k] = f_19 * kd_41[k]
                  + pb_x[k] * ld_41[k];

        t_99[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_40[k];

        t_100[k] = f_19 * kd_26[k]
                   + pb_z[k] * ld_40[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pa_y, pb_y, if0_9, if0_25, if1_20, \
                         if1_42, kd_28, kf_38, kf_54, ld_41, ld_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_y[k] * ld_41[k];

        t_102[k] = f_17 * if0_25[k]
                   - f_18 * if1_42[k]
                   + pa_x[k] * kf_54[k];

        t_103[k] = f_15 * if0_9[k]
                   - f_16 * if1_20[k]
                   + pa_y[k] * kf_38[k];

        t_104[k] = f_14 * kd_28[k]
                   + pb_y[k] * ld_42[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_x, pb_z, if0_26, if1_45, kd_43, \
                         kf_58, ld_42, ld_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_z[k] * ld_42[k];

        t_106[k] = f_5 * kd_43[k]
                   + pb_x[k] * ld_43[k];

        t_107[k] = f_12 * if0_26[k]
                   - f_13 * if1_45[k]
                   + pa_x[k] * kf_58[k];

        t_108[k] = pb_z[k] * ld_43[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pa_z, pb_y, pb_z, kd_28, kd_30, \
                         kf_38, kf_41, lp0_9, lp1_9, ld_44, ld_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * kd_30[k]
                   + pb_y[k] * ld_44[k];

        t_110[k] = f_1 * lp0_9[k]
                   - f_2 * lp1_9[k]
                   + pb_z[k] * ld_44[k];

        t_111[k] = pa_z[k] * kf_38[k];

        t_112[k] = f_3 * kd_28[k]
                   + pb_z[k] * ld_45[k];

        t_113[k] = pa_z[k] * kf_41[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_y, pb_z, kd_29, kd_30, kd_33, kf_43, \
                         ld_46, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * kd_29[k]
                   + pb_z[k] * ld_46[k];

        t_115[k] = f_19 * kd_33[k]
                   + pb_y[k] * ld_47[k];

        t_116[k] = f_5 * kd_30[k]
                   + pa_z[k] * kf_43[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_x, pa_y, pb_z, if0_12, if0_27, if1_24, \
                         if1_46, kd_31, kf_44, kf_62, ld_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_12 * if0_12[k]
                   - f_13 * if1_24[k]
                   + pa_y[k] * kf_44[k];

        t_118[k] = f_8 * kd_31[k]
                   + pb_z[k] * ld_48[k];

        t_119[k] = f_12 * if0_27[k]
                   - f_13 * if1_46[k]
                   + pa_x[k] * kf_62[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_x, pb_y, pb_z, if0_28, if1_47, kd_32, kd_36, \
                         kf_63, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * kd_32[k]
                   + pb_z[k] * ld_49[k];

        t_121[k] = f_5 * kd_36[k]
                   + pb_y[k] * ld_50[k];

        t_122[k] = f_12 * if0_28[k]
                   - f_13 * if1_47[k]
                   + pa_x[k] * kf_63[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pa_x, pa_y, pb_z, if0_13, if0_29, if1_25, \
                         if1_48, kd_34, kf_47, kf_65, ld_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_6 * if0_13[k]
                   - f_7 * if1_25[k]
                   + pa_y[k] * kf_47[k];

        t_124[k] = f_5 * kd_34[k]
                   + pb_z[k] * ld_51[k];

        t_125[k] = f_12 * if0_29[k]
                   - f_13 * if1_48[k]
                   + pa_x[k] * kf_65[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_x, pa_y, pb_y, pb_z, if0_30, if1_49, \
                         kd_35, kd_38, kf_48, kf_66, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_5 * kd_35[k]
                   + pb_z[k] * ld_52[k];

        t_127[k] = f_8 * kd_38[k]
                   + pb_y[k] * ld_53[k];

        t_128[k] = f_12 * if0_30[k]
                   - f_13 * if1_49[k]
                   + pa_x[k] * kf_66[k];

        t_129[k] = pa_y[k] * kf_48[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, kd_37, kd_40, \
                         kd_41, kf_50, kf_52, kf_54, ld_54, ld_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_y[k] * kf_50[k];

        t_131[k] = f_5 * kd_40[k]
                   + pa_y[k] * kf_52[k];

        t_132[k] = f_19 * kd_37[k]
                   + pb_z[k] * ld_54[k];

        t_133[k] = f_3 * kd_41[k]
                   + pb_y[k] * ld_55[k];

        t_134[k] = pa_y[k] * kf_54[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_x, pb_y, pb_z, if0_13, if1_25, \
                         kd_39, kd_58, kf_48, ld_56, ld_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_15 * if0_13[k]
                   - f_16 * if1_25[k]
                   + pa_z[k] * kf_48[k];

        t_136[k] = pb_y[k] * ld_56[k];

        t_137[k] = f_14 * kd_39[k]
                   + pb_z[k] * ld_56[k];

        t_138[k] = f_5 * kd_58[k]
                   + pb_x[k] * ld_58[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_y, pb_z, if0_31, if1_53, kd_40, \
                         kf_74, lp0_10, lp1_10, ld_57, ld_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_1 * lp0_10[k]
                   - f_2 * lp1_10[k]
                   + pb_y[k] * ld_57[k];

        t_140[k] = f_14 * kd_40[k]
                   + pb_z[k] * ld_57[k];

        t_141[k] = pb_y[k] * ld_58[k];

        t_142[k] = f_12 * if0_31[k]
                   - f_13 * if1_53[k]
                   + pa_x[k] * kf_74[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_y, pb_x, pb_y, pb_z, if0_16, if1_30, \
                         kd_42, kd_60, kf_55, ld_59, ld_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_10 * if0_16[k]
                   - f_11 * if1_30[k]
                   + pa_y[k] * kf_55[k];

        t_144[k] = f_9 * kd_42[k]
                   + pb_y[k] * ld_59[k];

        t_145[k] = pb_z[k] * ld_59[k];

        t_146[k] = f_8 * kd_60[k]
                   + pb_x[k] * ld_60[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_x, pb_y, pb_z, if0_32, if1_57, kd_44, \
                         kf_76, lp0_11, lp1_11, ld_60, ld_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_6 * if0_32[k]
                   - f_7 * if1_57[k]
                   + pa_x[k] * kf_76[k];

        t_148[k] = pb_z[k] * ld_60[k];

        t_149[k] = f_9 * kd_44[k]
                   + pb_y[k] * ld_61[k];

        t_150[k] = f_1 * lp0_11[k]
                   - f_2 * lp1_11[k]
                   + pb_z[k] * ld_61[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pa_z, pb_y, pb_z, kd_42, kd_43, \
                         kd_47, kf_55, kf_58, ld_62, ld_63, ld_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = pa_z[k] * kf_55[k];

        t_152[k] = f_3 * kd_42[k]
                   + pb_z[k] * ld_62[k];

        t_153[k] = pa_z[k] * kf_58[k];

        t_154[k] = f_3 * kd_43[k]
                   + pb_z[k] * ld_63[k];

        t_155[k] = f_14 * kd_47[k]
                   + pb_y[k] * ld_64[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_y, pa_z, pb_z, if0_19, if1_34, kd_44, kd_45, \
                         kf_60, kf_61, ld_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_5 * kd_44[k]
                   + pa_z[k] * kf_60[k];

        t_157[k] = f_17 * if0_19[k]
                   - f_18 * if1_34[k]
                   + pa_y[k] * kf_61[k];

        t_158[k] = f_8 * kd_45[k]
                   + pb_z[k] * ld_65[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_x, pb_y, pb_z, if0_34, if1_66, kd_46, kd_50, \
                         kf_77, ld_66, ld_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_6 * if0_34[k]
                   - f_7 * if1_66[k]
                   + pa_x[k] * kf_77[k];

        t_160[k] = f_8 * kd_46[k]
                   + pb_z[k] * ld_66[k];

        t_161[k] = f_19 * kd_50[k]
                   + pb_y[k] * ld_67[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pa_x, pa_y, pb_z, if0_22, if0_36, if1_37, \
                         if1_69, kd_48, kf_64, kf_78, ld_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_6 * if0_36[k]
                   - f_7 * if1_69[k]
                   + pa_x[k] * kf_78[k];

        t_163[k] = f_12 * if0_22[k]
                   - f_13 * if1_37[k]
                   + pa_y[k] * kf_64[k];

        t_164[k] = f_5 * kd_48[k]
                   + pb_z[k] * ld_68[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_x, pb_y, pb_z, if0_37, if1_71, kd_49, kd_53, \
                         kf_79, ld_69, ld_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_6 * if0_37[k]
                   - f_7 * if1_71[k]
                   + pa_x[k] * kf_79[k];

        t_166[k] = f_5 * kd_49[k]
                   + pb_z[k] * ld_69[k];

        t_167[k] = f_5 * kd_53[k]
                   + pb_y[k] * ld_70[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_x, pa_y, pb_z, if0_23, if0_39, if1_38, \
                         if1_74, kd_51, kf_67, kf_80, ld_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_6 * if0_39[k]
                   - f_7 * if1_74[k]
                   + pa_x[k] * kf_80[k];

        t_169[k] = f_6 * if0_23[k]
                   - f_7 * if1_38[k]
                   + pa_y[k] * kf_67[k];

        t_170[k] = f_19 * kd_51[k]
                   + pb_z[k] * ld_71[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_x, pb_y, pb_z, if0_40, if1_76, kd_52, kd_55, \
                         kf_81, ld_72, ld_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_6 * if0_40[k]
                   - f_7 * if1_76[k]
                   + pa_x[k] * kf_81[k];

        t_172[k] = f_19 * kd_52[k]
                   + pb_z[k] * ld_72[k];

        t_173[k] = f_8 * kd_55[k]
                   + pb_y[k] * ld_73[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_x, pa_y, if0_42, if1_79, kd_57, kf_68, \
                         kf_70, kf_72, kf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * if0_42[k]
                   - f_7 * if1_79[k]
                   + pa_x[k] * kf_82[k];

        t_175[k] = pa_y[k] * kf_68[k];

        t_176[k] = pa_y[k] * kf_70[k];

        t_177[k] = f_5 * kd_57[k]
                   + pa_y[k] * kf_72[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_y, pa_z, pb_y, pb_z, if0_23, if1_38, \
                         kd_54, kd_58, kf_68, kf_74, ld_74, ld_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * kd_54[k]
                   + pb_z[k] * ld_74[k];

        t_179[k] = f_3 * kd_58[k]
                   + pb_y[k] * ld_75[k];

        t_180[k] = pa_y[k] * kf_74[k];

        t_181[k] = f_10 * if0_23[k]
                   - f_11 * if1_38[k]
                   + pa_z[k] * kf_68[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, pb_x, pb_y, pb_z, kd_56, kd_57, \
                         kd_66, lp0_12, lp1_12, ld_76, ld_77, ld_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = pb_y[k] * ld_76[k];

        t_183[k] = f_9 * kd_56[k]
                   + pb_z[k] * ld_76[k];

        t_184[k] = f_8 * kd_66[k]
                   + pb_x[k] * ld_78[k];

        t_185[k] = f_1 * lp0_12[k]
                   - f_2 * lp1_12[k]
                   + pb_y[k] * ld_77[k];

        t_186[k] = f_9 * kd_57[k]
                   + pb_z[k] * ld_77[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pa_x, pb_y, if0_44, if1_91, kd_59, kd_67, \
                         kf_85, kf_86, ld_78, ld_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = pb_y[k] * ld_78[k];

        t_188[k] = f_6 * if0_44[k]
                   - f_7 * if1_91[k]
                   + pa_x[k] * kf_85[k];

        t_189[k] = f_5 * kd_67[k]
                   + pa_x[k] * kf_86[k];

        t_190[k] = f_4 * kd_59[k]
                   + pb_y[k] * ld_79[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pa_x, pa_z, pb_x, kd_68, kf_75, \
                         kf_89, kf_91, kf_92, ld_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_3 * kd_68[k]
                   + pb_x[k] * ld_80[k];

        t_192[k] = pa_x[k] * kf_89[k];

        t_193[k] = pa_x[k] * kf_91[k];

        t_194[k] = pa_x[k] * kf_92[k];

        t_195[k] = pa_z[k] * kf_75[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_x, pb_z, kd_59, kd_73, kf_94, \
                         kf_95, kf_96, kf_97, ld_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_3 * kd_59[k]
                   + pb_z[k] * ld_81[k];

        t_197[k] = pa_x[k] * kf_94[k];

        t_198[k] = pa_x[k] * kf_95[k];

        t_199[k] = pa_x[k] * kf_96[k];

        t_200[k] = f_5 * kd_73[k]
                   + pa_x[k] * kf_97[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pa_x, pb_z, kd_61, kd_76, \
                         kf_100, kf_101, kf_102, kf_103, kf_104, \
                         ld_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_8 * kd_61[k]
                   + pb_z[k] * ld_82[k];

        t_202[k] = pa_x[k] * kf_100[k];

        t_203[k] = pa_x[k] * kf_101[k];

        t_204[k] = pa_x[k] * kf_102[k];

        t_205[k] = pa_x[k] * kf_103[k];

        t_206[k] = f_5 * kd_76[k]
                   + pa_x[k] * kf_104[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, pa_x, pb_z, kd_62, kd_79, \
                         kf_107, kf_108, kf_109, kf_110, kf_111, \
                         ld_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_5 * kd_62[k]
                   + pb_z[k] * ld_83[k];

        t_208[k] = pa_x[k] * kf_107[k];

        t_209[k] = pa_x[k] * kf_108[k];

        t_210[k] = pa_x[k] * kf_109[k];

        t_211[k] = pa_x[k] * kf_110[k];

        t_212[k] = f_5 * kd_79[k]
                   + pa_x[k] * kf_111[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, pa_x, pb_z, kd_63, kd_82, \
                         kf_114, kf_115, kf_116, kf_117, kf_118, \
                         ld_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_19 * kd_63[k]
                   + pb_z[k] * ld_84[k];

        t_214[k] = pa_x[k] * kf_114[k];

        t_215[k] = pa_x[k] * kf_115[k];

        t_216[k] = pa_x[k] * kf_116[k];

        t_217[k] = pa_x[k] * kf_117[k];

        t_218[k] = f_5 * kd_82[k]
                   + pa_x[k] * kf_118[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, pa_x, pa_y, pb_z, kd_64, \
                         kf_83, kf_121, kf_122, kf_123, kf_124, ld_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_14 * kd_64[k]
                   + pb_z[k] * ld_85[k];

        t_220[k] = pa_x[k] * kf_121[k];

        t_221[k] = pa_x[k] * kf_122[k];

        t_222[k] = pa_x[k] * kf_123[k];

        t_223[k] = pa_x[k] * kf_124[k];

        t_224[k] = pa_y[k] * kf_83[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, pa_x, pa_y, kd_87, kf_84, kf_125, \
                         kf_126, kf_127, kf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_y[k] * kf_84[k];

        t_226[k] = pa_x[k] * kf_125[k];

        t_227[k] = pa_x[k] * kf_126[k];

        t_228[k] = pa_x[k] * kf_127[k];

        t_229[k] = f_5 * kd_87[k]
                   + pa_x[k] * kf_129[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, pa_x, pb_x, pb_z, kd_65, kd_89, \
                         kf_133, kf_134, kf_136, ld_86, ld_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * kd_65[k]
                   + pb_z[k] * ld_86[k];

        t_231[k] = f_3 * kd_89[k]
                   + pb_x[k] * ld_87[k];

        t_232[k] = pa_x[k] * kf_133[k];

        t_233[k] = pa_x[k] * kf_134[k];

        t_234[k] = pa_x[k] * kf_136[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pb_x, pb_y, kd_67, kd_68, lp0_13, \
                         lp0_14, lp1_13, lp1_14, ld_88, ld_89, ld_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_1 * lp0_13[k]
                   - f_2 * lp1_13[k]
                   + pb_x[k] * ld_88[k];

        t_236[k] = f_0 * kd_67[k]
                   + pb_y[k] * ld_88[k];

        t_237[k] = pb_x[k] * ld_89[k];

        t_238[k] = pb_x[k] * ld_90[k];

        t_239[k] = f_0 * kd_68[k]
                   + f_1 * lp0_14[k]
                   - f_2 * lp1_14[k]
                   + pb_y[k] * ld_89[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pa_z, pb_y, pb_z, kd_67, kd_69, \
                         kf_86, lp0_15, lp1_15, ld_89, ld_90, ld_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pb_z[k] * ld_89[k];

        t_241[k] = f_0 * kd_69[k]
                   + pb_y[k] * ld_90[k];

        t_242[k] = f_1 * lp0_15[k]
                   - f_2 * lp1_15[k]
                   + pb_z[k] * ld_90[k];

        t_243[k] = pa_z[k] * kf_86[k];

        t_244[k] = f_3 * kd_67[k]
                   + pb_z[k] * ld_91[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_z, pb_y, pb_z, kd_68, kd_69, kd_72, \
                         kf_89, kf_92, ld_92, ld_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pa_z[k] * kf_89[k];

        t_246[k] = f_3 * kd_68[k]
                   + pb_z[k] * ld_92[k];

        t_247[k] = f_4 * kd_72[k]
                   + pb_y[k] * ld_93[k];

        t_248[k] = f_5 * kd_69[k]
                   + pa_z[k] * kf_92[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pb_x, pb_z, kd_70, lp0_16, lp1_16, ld_94, \
                         ld_95, ld_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_1 * lp0_16[k]
                   - f_2 * lp1_16[k]
                   + pb_x[k] * ld_94[k];

        t_250[k] = f_8 * kd_70[k]
                   + pb_z[k] * ld_94[k];

        t_251[k] = pb_x[k] * ld_95[k];

        t_252[k] = pb_x[k] * ld_96[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, pa_z, pb_y, pb_z, if0_32, if1_57, kd_71, kd_75, \
                         kf_93, ld_95, ld_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_6 * if0_32[k]
                   - f_7 * if1_57[k]
                   + pa_z[k] * kf_93[k];

        t_254[k] = f_8 * kd_71[k]
                   + pb_z[k] * ld_95[k];

        t_255[k] = f_9 * kd_75[k]
                   + pb_y[k] * ld_96[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_y, pb_x, pb_z, if0_36, if1_69, kd_73, \
                         kf_103, lp0_17, lp1_17, ld_97, ld_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_10 * if0_36[k]
                   - f_11 * if1_69[k]
                   + pa_y[k] * kf_103[k];

        t_257[k] = f_1 * lp0_17[k]
                   - f_2 * lp1_17[k]
                   + pb_x[k] * ld_97[k];

        t_258[k] = f_5 * kd_73[k]
                   + pb_z[k] * ld_97[k];

        t_259[k] = pb_x[k] * ld_98[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_z, pb_x, pb_y, pb_z, if0_33, if1_61, \
                         kd_74, kd_78, kf_100, ld_98, ld_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pb_x[k] * ld_99[k];

        t_261[k] = f_12 * if0_33[k]
                   - f_13 * if1_61[k]
                   + pa_z[k] * kf_100[k];

        t_262[k] = f_5 * kd_74[k]
                   + pb_z[k] * ld_98[k];

        t_263[k] = f_14 * kd_78[k]
                   + pb_y[k] * ld_99[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_y, pb_x, pb_z, if0_39, if1_74, kd_76, \
                         kf_110, lp0_18, lp1_18, ld_100, ld_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_15 * if0_39[k]
                   - f_16 * if1_74[k]
                   + pa_y[k] * kf_110[k];

        t_265[k] = f_1 * lp0_18[k]
                   - f_2 * lp1_18[k]
                   + pb_x[k] * ld_100[k];

        t_266[k] = f_19 * kd_76[k]
                   + pb_z[k] * ld_100[k];

        t_267[k] = pb_x[k] * ld_101[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_z, pb_x, pb_y, pb_z, if0_34, if1_66, \
                         kd_77, kd_81, kf_107, ld_101, ld_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = pb_x[k] * ld_102[k];

        t_269[k] = f_17 * if0_34[k]
                   - f_18 * if1_66[k]
                   + pa_z[k] * kf_107[k];

        t_270[k] = f_19 * kd_77[k]
                   + pb_z[k] * ld_101[k];

        t_271[k] = f_19 * kd_81[k]
                   + pb_y[k] * ld_102[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_y, pb_x, pb_z, if0_42, if1_79, kd_79, \
                         kf_117, lp0_19, lp1_19, ld_103, ld_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * if0_42[k]
                   - f_18 * if1_79[k]
                   + pa_y[k] * kf_117[k];

        t_273[k] = f_1 * lp0_19[k]
                   - f_2 * lp1_19[k]
                   + pb_x[k] * ld_103[k];

        t_274[k] = f_14 * kd_79[k]
                   + pb_z[k] * ld_103[k];

        t_275[k] = pb_x[k] * ld_104[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_z, pb_x, pb_y, pb_z, if0_37, if1_71, \
                         kd_80, kd_84, kf_114, ld_104, ld_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_x[k] * ld_105[k];

        t_277[k] = f_15 * if0_37[k]
                   - f_16 * if1_71[k]
                   + pa_z[k] * kf_114[k];

        t_278[k] = f_14 * kd_80[k]
                   + pb_z[k] * ld_104[k];

        t_279[k] = f_5 * kd_84[k]
                   + pb_y[k] * ld_105[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_z, if0_43, if1_83, kd_82, \
                         kf_124, lp0_20, lp1_20, ld_106, ld_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_12 * if0_43[k]
                   - f_13 * if1_83[k]
                   + pa_y[k] * kf_124[k];

        t_281[k] = f_1 * lp0_20[k]
                   - f_2 * lp1_20[k]
                   + pb_x[k] * ld_106[k];

        t_282[k] = f_9 * kd_82[k]
                   + pb_z[k] * ld_106[k];

        t_283[k] = pb_x[k] * ld_107[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pa_z, pb_x, pb_y, pb_z, if0_40, if1_76, \
                         kd_83, kd_86, kf_121, ld_107, ld_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pb_x[k] * ld_108[k];

        t_285[k] = f_10 * if0_40[k]
                   - f_11 * if1_76[k]
                   + pa_z[k] * kf_121[k];

        t_286[k] = f_9 * kd_83[k]
                   + pb_z[k] * ld_107[k];

        t_287[k] = f_8 * kd_86[k]
                   + pb_y[k] * ld_108[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pa_y, pb_z, if0_44, if1_91, kd_85, \
                         kd_88, kf_128, kf_129, kf_130, kf_133, \
                         ld_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * if0_44[k]
                   - f_7 * if1_91[k]
                   + pa_y[k] * kf_128[k];

        t_289[k] = pa_y[k] * kf_129[k];

        t_290[k] = pa_y[k] * kf_130[k];

        t_291[k] = f_5 * kd_88[k]
                   + pa_y[k] * kf_133[k];

        t_292[k] = f_4 * kd_85[k]
                   + pb_z[k] * ld_109[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_y, pb_x, pb_y, pb_z, kd_87, kd_89, \
                         kf_136, lp0_21, lp1_21, ld_110, ld_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_3 * kd_89[k]
                   + pb_y[k] * ld_110[k];

        t_294[k] = pa_y[k] * kf_136[k];

        t_295[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_111[k];

        t_296[k] = f_0 * kd_87[k]
                   + pb_z[k] * ld_111[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, pb_x, pb_y, pb_z, kd_88, lp0_22, \
                         lp1_22, ld_112, ld_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pb_x[k] * ld_112[k];

        t_298[k] = pb_x[k] * ld_113[k];

        t_299[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_112[k];

        t_300[k] = f_0 * kd_88[k]
                   + pb_z[k] * ld_112[k];

        t_301[k] = pb_y[k] * ld_113[k];
    }

#pragma omp simd aligned(t_302, pb_z, kd_89, lp0_23, lp1_23, ld_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * kd_89[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_113[k];
    }
}

auto
compute_prim_lf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_1 = buffer.data(if0 + 1);
    const auto *if0_2 = buffer.data(if0 + 2);
    const auto *if0_3 = buffer.data(if0 + 3);
    const auto *if0_5 = buffer.data(if0 + 5);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_9 = buffer.data(if0 + 9);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_12 = buffer.data(if0 + 12);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_15 = buffer.data(if0 + 15);
    const auto *if0_17 = buffer.data(if0 + 17);
    const auto *if0_18 = buffer.data(if0 + 18);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_24 = buffer.data(if0 + 24);
    const auto *if0_25 = buffer.data(if0 + 25);
    const auto *if0_27 = buffer.data(if0 + 27);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_30 = buffer.data(if0 + 30);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_33 = buffer.data(if0 + 33);
    const auto *if0_34 = buffer.data(if0 + 34);
    const auto *if0_35 = buffer.data(if0 + 35);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_1 = buffer.data(if1 + 1);
    const auto *if1_2 = buffer.data(if1 + 2);
    const auto *if1_3 = buffer.data(if1 + 3);
    const auto *if1_5 = buffer.data(if1 + 5);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_9 = buffer.data(if1 + 9);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_12 = buffer.data(if1 + 12);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_15 = buffer.data(if1 + 15);
    const auto *if1_17 = buffer.data(if1 + 17);
    const auto *if1_18 = buffer.data(if1 + 18);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_21 = buffer.data(if1 + 21);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_24 = buffer.data(if1 + 24);
    const auto *if1_25 = buffer.data(if1 + 25);
    const auto *if1_27 = buffer.data(if1 + 27);
    const auto *if1_28 = buffer.data(if1 + 28);
    const auto *if1_30 = buffer.data(if1 + 30);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_33 = buffer.data(if1 + 33);
    const auto *if1_34 = buffer.data(if1 + 34);
    const auto *if1_35 = buffer.data(if1 + 35);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_47 = buffer.data(kd + 47);

    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_4, kf_6, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_6[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_4[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_5, if1_5, kf_11, lp0_3, lp1_3, ld_4, \
                         ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_5[k]
                 - f_7 * if1_5[k]
                 + pa_x[k] * kf_11[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_8, kf_7, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_7[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_8[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_1, if0_8, if1_1, \
                         if1_8, kf_8, kf_19, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_8[k]
                  - f_7 * if1_8[k]
                  + pa_x[k] * kf_19[k];

        t_18[k] = f_8 * if0_1[k]
                  - f_9 * if1_1[k]
                  + pa_y[k] * kf_8[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_11, if1_11, kd_10, \
                         kf_23, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_10[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_11[k]
                  - f_12 * if1_11[k]
                  + pa_x[k] * kf_23[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_2, if1_2, kd_14, kf_14, \
                         lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_2[k]
                  - f_9 * if1_2[k]
                  + pa_z[k] * kf_14[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_14[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_3, if0_14, if1_3, \
                         if1_14, kf_20, kf_31, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_14[k]
                  - f_12 * if1_14[k]
                  + pa_x[k] * kf_31[k];

        t_30[k] = f_13 * if0_3[k]
                  - f_14 * if1_3[k]
                  + pa_y[k] * kf_20[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_17, if1_17, kd_16, \
                         kf_35, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_16[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_17[k]
                  - f_14 * if1_17[k]
                  + pa_x[k] * kf_35[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_6, if1_6, kd_20, kf_26, \
                         lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_6[k]
                  - f_14 * if1_6[k]
                  + pa_z[k] * kf_26[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_20[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_9, if0_20, if1_9, \
                         if1_20, kf_32, kf_43, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_20[k]
                  - f_14 * if1_20[k]
                  + pa_x[k] * kf_43[k];

        t_42[k] = f_11 * if0_9[k]
                  - f_12 * if1_9[k]
                  + pa_y[k] * kf_32[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_21, if1_21, kd_22, \
                         kf_47, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_22[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_21[k]
                  - f_9 * if1_21[k]
                  + pa_x[k] * kf_47[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_12, if1_12, kd_26, \
                         kf_38, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_12[k]
                  - f_12 * if1_12[k]
                  + pa_z[k] * kf_38[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_26[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_15, if0_22, \
                         if1_15, if1_22, kf_44, kf_55, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_22[k]
                  - f_9 * if1_22[k]
                  + pa_x[k] * kf_55[k];

        t_54[k] = f_6 * if0_15[k]
                  - f_7 * if1_15[k]
                  + pa_y[k] * kf_44[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_23, if1_23, kd_27, \
                         kf_56, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_27[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_23[k]
                  - f_4 * if1_23[k]
                  + pa_x[k] * kf_56[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_18, if1_18, kd_28, \
                         kf_50, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_18[k]
                  - f_7 * if1_18[k]
                  + pa_z[k] * kf_50[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_28[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_35, if1_35, kf_57, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_35[k]
                  - f_4 * if1_35[k]
                  + pa_x[k] * kf_57[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_30, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_30[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_23, if1_23, kf_64, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_23[k]
                  - f_4 * if1_23[k]
                  + pa_z[k] * kf_64[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_27, if1_27, kd_34, \
                         kf_70, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_34[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_27[k]
                  - f_7 * if1_27[k]
                  + pa_y[k] * kf_70[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_24, if0_30, \
                         if1_24, if1_30, kd_37, kf_68, kf_76, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_24[k]
                  - f_9 * if1_24[k]
                  + pa_z[k] * kf_68[k];

        t_82[k] = f_10 * kd_37[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_30[k]
                  - f_12 * if1_30[k]
                  + pa_y[k] * kf_76[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_25, if1_25, kf_74, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_25[k]
                  - f_14 * if1_25[k]
                  + pa_z[k] * kf_74[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_33, if1_33, kd_40, \
                         kf_82, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_40[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_33[k]
                  - f_14 * if1_33[k]
                  + pa_y[k] * kf_82[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_28, if0_34, \
                         if1_28, if1_34, kd_43, kf_80, kf_88, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_28[k]
                  - f_12 * if1_28[k]
                  + pa_z[k] * kf_80[k];

        t_94[k] = f_16 * kd_43[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_34[k]
                  - f_9 * if1_34[k]
                  + pa_y[k] * kf_88[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_31, if1_31, kf_86, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_31[k]
                  - f_7 * if1_31[k]
                  + pa_z[k] * kf_86[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_35, if1_35, kd_44, \
                         kf_89, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_44[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_35[k]
                   - f_4 * if1_35[k]
                   + pa_y[k] * kf_89[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_47, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_47[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_1 = buffer.data(if0 + 1);
    const auto *if0_2 = buffer.data(if0 + 2);
    const auto *if0_3 = buffer.data(if0 + 3);
    const auto *if0_5 = buffer.data(if0 + 5);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_9 = buffer.data(if0 + 9);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_12 = buffer.data(if0 + 12);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_15 = buffer.data(if0 + 15);
    const auto *if0_17 = buffer.data(if0 + 17);
    const auto *if0_18 = buffer.data(if0 + 18);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_24 = buffer.data(if0 + 24);
    const auto *if0_25 = buffer.data(if0 + 25);
    const auto *if0_27 = buffer.data(if0 + 27);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_30 = buffer.data(if0 + 30);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_33 = buffer.data(if0 + 33);
    const auto *if0_34 = buffer.data(if0 + 34);
    const auto *if0_35 = buffer.data(if0 + 35);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_10 = buffer.data(if1 + 10);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_17 = buffer.data(if1 + 17);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_27 = buffer.data(if1 + 27);
    const auto *if1_28 = buffer.data(if1 + 28);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_40 = buffer.data(if1 + 40);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_46 = buffer.data(if1 + 46);
    const auto *if1_49 = buffer.data(if1 + 49);
    const auto *if1_61 = buffer.data(if1 + 61);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_69 = buffer.data(if1 + 69);
    const auto *if1_79 = buffer.data(if1 + 79);
    const auto *if1_84 = buffer.data(if1 + 84);
    const auto *if1_89 = buffer.data(if1 + 89);
    const auto *if1_94 = buffer.data(if1 + 94);
    const auto *if1_96 = buffer.data(if1 + 96);
    const auto *if1_100 = buffer.data(if1 + 100);
    const auto *if1_102 = buffer.data(if1 + 102);
    const auto *if1_106 = buffer.data(if1 + 106);
    const auto *if1_108 = buffer.data(if1 + 108);
    const auto *if1_112 = buffer.data(if1 + 112);
    const auto *if1_119 = buffer.data(if1 + 119);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_136 = buffer.data(kf + 136);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_6, kf_7, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_7[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_6[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_5, if1_17, kf_17, lp0_3, lp1_3, \
                         ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_5[k]
                 - f_7 * if1_17[k]
                 + pa_x[k] * kf_17[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_10, kf_10, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_10[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_10[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_1, if0_8, if1_7, \
                         if1_27, kf_14, kf_27, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_8[k]
                  - f_7 * if1_27[k]
                  + pa_x[k] * kf_27[k];

        t_18[k] = f_8 * if0_1[k]
                  - f_9 * if1_7[k]
                  + pa_y[k] * kf_14[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_11, if1_31, kd_12, \
                         kf_31, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_12[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_11[k]
                  - f_12 * if1_31[k]
                  + pa_x[k] * kf_31[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_2, if1_10, kd_16, \
                         kf_22, lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_2[k]
                  - f_9 * if1_10[k]
                  + pa_z[k] * kf_22[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_16[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_3, if0_14, \
                         if1_14, if1_45, kf_28, kf_44, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_14[k]
                  - f_12 * if1_45[k]
                  + pa_x[k] * kf_44[k];

        t_30[k] = f_13 * if0_3[k]
                  - f_14 * if1_14[k]
                  + pa_y[k] * kf_28[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_17, if1_49, kd_18, \
                         kf_48, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_18[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_17[k]
                  - f_14 * if1_49[k]
                  + pa_x[k] * kf_48[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_6, if1_22, kd_22, \
                         kf_39, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_6[k]
                  - f_14 * if1_22[k]
                  + pa_z[k] * kf_39[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_22[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_9, if0_20, \
                         if1_28, if1_66, kf_45, kf_64, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_20[k]
                  - f_14 * if1_66[k]
                  + pa_x[k] * kf_64[k];

        t_42[k] = f_11 * if0_9[k]
                  - f_12 * if1_28[k]
                  + pa_y[k] * kf_45[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_21, if1_69, kd_24, \
                         kf_68, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_24[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_21[k]
                  - f_9 * if1_69[k]
                  + pa_x[k] * kf_68[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_12, if1_40, kd_28, \
                         kf_59, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_12[k]
                  - f_12 * if1_40[k]
                  + pa_z[k] * kf_59[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_28[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_15, if0_22, \
                         if1_46, if1_79, kf_65, kf_87, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_22[k]
                  - f_9 * if1_79[k]
                  + pa_x[k] * kf_87[k];

        t_54[k] = f_6 * if0_15[k]
                  - f_7 * if1_46[k]
                  + pa_y[k] * kf_65[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_23, if1_84, kd_29, \
                         kf_90, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_29[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_23[k]
                  - f_4 * if1_84[k]
                  + pa_x[k] * kf_90[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_18, if1_61, kd_30, \
                         kf_82, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_18[k]
                  - f_7 * if1_61[k]
                  + pa_z[k] * kf_82[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_35, if1_119, kf_97, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_35[k]
                  - f_4 * if1_119[k]
                  + pa_x[k] * kf_97[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_32, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_32[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_23, if1_84, kf_107, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_23[k]
                  - f_4 * if1_84[k]
                  + pa_z[k] * kf_107[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_27, if1_96, kd_37, \
                         kf_114, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_37[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_27[k]
                  - f_7 * if1_96[k]
                  + pa_y[k] * kf_114[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_24, if0_30, \
                         if1_89, if1_102, kd_40, kf_112, kf_120, \
                         ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_24[k]
                  - f_9 * if1_89[k]
                  + pa_z[k] * kf_112[k];

        t_82[k] = f_10 * kd_40[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_30[k]
                  - f_12 * if1_102[k]
                  + pa_y[k] * kf_120[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_25, if1_94, kf_118, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_25[k]
                  - f_14 * if1_94[k]
                  + pa_z[k] * kf_118[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_33, if1_108, kd_43, \
                         kf_126, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_43[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_33[k]
                  - f_14 * if1_108[k]
                  + pa_y[k] * kf_126[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_28, if0_34, \
                         if1_100, if1_112, kd_46, kf_124, kf_132, \
                         ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_28[k]
                  - f_12 * if1_100[k]
                  + pa_z[k] * kf_124[k];

        t_94[k] = f_16 * kd_46[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_34[k]
                  - f_9 * if1_112[k]
                  + pa_y[k] * kf_132[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_31, if1_106, kf_130, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_31[k]
                  - f_7 * if1_106[k]
                  + pa_z[k] * kf_130[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_35, if1_119, kd_47, \
                         kf_136, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_47[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_35[k]
                   - f_4 * if1_119[k]
                   + pa_y[k] * kf_136[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_50, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_50[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 3.0 / p;
    const auto f_7 = 2.5 / alpha;
    const auto f_8 = 2.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / alpha;
    const auto f_13 = 2.0 * beta / (alpha * p);
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 2.0 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_10 = buffer.data(if0 + 10);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_17 = buffer.data(if0 + 17);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_27 = buffer.data(if0 + 27);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_37 = buffer.data(if0 + 37);
    const auto *if0_40 = buffer.data(if0 + 40);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_46 = buffer.data(if0 + 46);
    const auto *if0_49 = buffer.data(if0 + 49);
    const auto *if0_55 = buffer.data(if0 + 55);
    const auto *if0_56 = buffer.data(if0 + 56);
    const auto *if0_57 = buffer.data(if0 + 57);
    const auto *if0_58 = buffer.data(if0 + 58);
    const auto *if0_61 = buffer.data(if0 + 61);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_69 = buffer.data(if0 + 69);
    const auto *if0_72 = buffer.data(if0 + 72);
    const auto *if0_73 = buffer.data(if0 + 73);
    const auto *if0_75 = buffer.data(if0 + 75);
    const auto *if0_76 = buffer.data(if0 + 76);
    const auto *if0_79 = buffer.data(if0 + 79);
    const auto *if0_84 = buffer.data(if0 + 84);
    const auto *if0_89 = buffer.data(if0 + 89);
    const auto *if0_94 = buffer.data(if0 + 94);
    const auto *if0_96 = buffer.data(if0 + 96);
    const auto *if0_100 = buffer.data(if0 + 100);
    const auto *if0_102 = buffer.data(if0 + 102);
    const auto *if0_106 = buffer.data(if0 + 106);
    const auto *if0_108 = buffer.data(if0 + 108);
    const auto *if0_112 = buffer.data(if0 + 112);
    const auto *if0_119 = buffer.data(if0 + 119);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_9 = buffer.data(if1 + 9);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_17 = buffer.data(if1 + 17);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_29 = buffer.data(if1 + 29);
    const auto *if1_30 = buffer.data(if1 + 30);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_36 = buffer.data(if1 + 36);
    const auto *if1_39 = buffer.data(if1 + 39);
    const auto *if1_42 = buffer.data(if1 + 42);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_44 = buffer.data(if1 + 44);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_46 = buffer.data(if1 + 46);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_55 = buffer.data(if1 + 55);
    const auto *if1_56 = buffer.data(if1 + 56);
    const auto *if1_57 = buffer.data(if1 + 57);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_61 = buffer.data(if1 + 61);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_69 = buffer.data(if1 + 69);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_76 = buffer.data(if1 + 76);
    const auto *if1_80 = buffer.data(if1 + 80);
    const auto *if1_82 = buffer.data(if1 + 82);
    const auto *if1_86 = buffer.data(if1 + 86);
    const auto *if1_88 = buffer.data(if1 + 88);
    const auto *if1_91 = buffer.data(if1 + 91);
    const auto *if1_98 = buffer.data(if1 + 98);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_119 = buffer.data(kf + 119);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pa_z, pb_z, kd_1, kf_0, kf_3, \
                         kf_5, lp0_2, lp1_2, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = pa_y[k] * kf_0[k];

        t_7[k] = f_3 * kd_1[k]
                 + pa_y[k] * kf_3[k];

        t_8[k] = pa_y[k] * kf_5[k];

        t_9[k] = pa_z[k] * kf_0[k];

        t_10[k] = pa_z[k] * kf_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_y, pa_z, pb_x, pb_z, if0_0, if1_0, kd_2, \
                         kd_7, kf_5, kf_6, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * kd_2[k]
                  + pa_z[k] * kf_5[k];

        t_12[k] = f_4 * if0_0[k]
                  - f_5 * if1_0[k]
                  + pa_y[k] * kf_6[k];

        t_13[k] = pb_z[k] * ld_3[k];

        t_14[k] = f_6 * kd_7[k]
                  + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_z, pb_z, if0_17, if1_14, kf_7, \
                         kf_13, lp0_3, lp1_3, ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * if0_17[k]
                  - f_8 * if1_14[k]
                  + pa_x[k] * kf_13[k];

        t_16[k] = pb_z[k] * ld_4[k];

        t_17[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];

        t_18[k] = pa_z[k] * kf_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_x, pb_y, if0_0, if1_0, kd_11, \
                         kf_8, kf_9, ld_6, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * kf_9[k];

        t_20[k] = f_4 * if0_0[k]
                  - f_5 * if1_0[k]
                  + pa_z[k] * kf_8[k];

        t_21[k] = pb_y[k] * ld_6[k];

        t_22[k] = f_6 * kd_11[k]
                  + pb_x[k] * ld_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_y, if0_27, if1_22, kf_21, lp0_4, lp1_4, \
                         ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];

        t_24[k] = pb_y[k] * ld_8[k];

        t_25[k] = f_7 * if0_27[k]
                  - f_8 * if1_22[k]
                  + pa_x[k] * kf_21[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, pb_z, if0_7, if1_7, kd_13, kf_10, ld_9, \
                         ld_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_9 * if0_7[k]
                  - f_10 * if1_7[k]
                  + pa_y[k] * kf_10[k];

        t_27[k] = pb_z[k] * ld_9[k];

        t_28[k] = f_11 * kd_13[k]
                  + pb_x[k] * ld_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, if0_31, if1_26, kf_10, \
                         kf_25, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_12 * if0_31[k]
                  - f_13 * if1_26[k]
                  + pa_x[k] * kf_25[k];

        t_30[k] = pb_z[k] * ld_10[k];

        t_31[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];

        t_32[k] = pa_z[k] * kf_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, if0_10, if1_9, kd_8, kd_10, \
                         kf_13, kf_15, kf_16, kf_19, kf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * kf_13[k];

        t_34[k] = f_3 * kd_8[k]
                  + pa_z[k] * kf_15[k];

        t_35[k] = f_3 * kd_10[k]
                  + pa_y[k] * kf_19[k];

        t_36[k] = pa_y[k] * kf_21[k];

        t_37[k] = f_9 * if0_10[k]
                  - f_10 * if1_9[k]
                  + pa_z[k] * kf_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, kd_17, lp0_6, lp1_6, ld_12, \
                         ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * ld_12[k];

        t_39[k] = f_11 * kd_17[k]
                  + pb_x[k] * ld_14[k];

        t_40[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];

        t_41[k] = pb_y[k] * ld_14[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pa_y, pb_z, if0_14, if0_45, if1_11, if1_35, \
                         kf_22, kf_34, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_12 * if0_45[k]
                  - f_13 * if1_35[k]
                  + pa_x[k] * kf_34[k];

        t_43[k] = f_14 * if0_14[k]
                  - f_15 * if1_11[k]
                  + pa_y[k] * kf_22[k];

        t_44[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, if0_49, if1_39, kd_19, \
                         kf_38, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_16 * kd_19[k]
                  + pb_x[k] * ld_16[k];

        t_46[k] = f_14 * if0_49[k]
                  - f_15 * if1_39[k]
                  + pa_x[k] * kf_38[k];

        t_47[k] = pb_z[k] * ld_16[k];

        t_48[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pa_z, if0_22, if1_17, kd_14, kf_22, \
                         kf_25, kf_27, kf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_z[k] * kf_22[k];

        t_50[k] = pa_z[k] * kf_25[k];

        t_51[k] = f_3 * kd_14[k]
                  + pa_z[k] * kf_27[k];

        t_52[k] = f_4 * if0_22[k]
                  - f_5 * if1_17[k]
                  + pa_y[k] * kf_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pa_y, if0_56, if0_57, if1_43, if1_44, \
                         kd_16, kf_32, kf_34, kf_42, kf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_14 * if0_56[k]
                  - f_15 * if1_43[k]
                  + pa_x[k] * kf_42[k];

        t_54[k] = f_14 * if0_57[k]
                  - f_15 * if1_44[k]
                  + pa_x[k] * kf_43[k];

        t_55[k] = f_3 * kd_16[k]
                  + pa_y[k] * kf_32[k];

        t_56[k] = pa_y[k] * kf_34[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_z, pb_x, pb_y, if0_22, if1_17, kd_23, \
                         kf_29, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_14 * if0_22[k]
                  - f_15 * if1_17[k]
                  + pa_z[k] * kf_29[k];

        t_58[k] = pb_y[k] * ld_18[k];

        t_59[k] = f_16 * kd_23[k]
                  + pb_x[k] * ld_20[k];

        t_60[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pa_y, pb_y, pb_z, if0_28, if0_66, \
                         if1_23, if1_51, kf_35, kf_50, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_y[k] * ld_20[k];

        t_62[k] = f_14 * if0_66[k]
                  - f_15 * if1_51[k]
                  + pa_x[k] * kf_50[k];

        t_63[k] = f_12 * if0_28[k]
                  - f_13 * if1_23[k]
                  + pa_y[k] * kf_35[k];

        t_64[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_x, pb_z, if0_69, if1_54, kd_25, \
                         kf_54, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * kd_25[k]
                  + pb_x[k] * ld_22[k];

        t_66[k] = f_9 * if0_69[k]
                  - f_10 * if1_54[k]
                  + pa_x[k] * kf_54[k];

        t_67[k] = pb_z[k] * ld_22[k];

        t_68[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pa_z, if0_37, if1_29, kd_20, kf_35, \
                         kf_38, kf_40, kf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_z[k] * kf_35[k];

        t_70[k] = pa_z[k] * kf_38[k];

        t_71[k] = f_3 * kd_20[k]
                  + pa_z[k] * kf_40[k];

        t_72[k] = f_9 * if0_37[k]
                  - f_10 * if1_29[k]
                  + pa_y[k] * kf_41[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_y, if0_40, if0_72, if0_73, if1_30, if1_55, \
                         if1_56, kf_44, kf_58, kf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_9 * if0_72[k]
                  - f_10 * if1_55[k]
                  + pa_x[k] * kf_58[k];

        t_74[k] = f_9 * if0_73[k]
                  - f_10 * if1_56[k]
                  + pa_x[k] * kf_59[k];

        t_75[k] = f_4 * if0_40[k]
                  - f_5 * if1_30[k]
                  + pa_y[k] * kf_44[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_x, pa_y, if0_75, if0_76, if1_57, if1_58, \
                         kd_22, kf_48, kf_50, kf_61, kf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_9 * if0_75[k]
                  - f_10 * if1_57[k]
                  + pa_x[k] * kf_61[k];

        t_77[k] = f_9 * if0_76[k]
                  - f_10 * if1_58[k]
                  + pa_x[k] * kf_62[k];

        t_78[k] = f_3 * kd_22[k]
                  + pa_y[k] * kf_48[k];

        t_79[k] = pa_y[k] * kf_50[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_z, pb_x, pb_y, if0_40, if1_30, kd_29, \
                         kf_45, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * if0_40[k]
                  - f_13 * if1_30[k]
                  + pa_z[k] * kf_45[k];

        t_81[k] = pb_y[k] * ld_24[k];

        t_82[k] = f_3 * kd_29[k]
                  + pb_x[k] * ld_26[k];

        t_83[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pa_y, pb_y, pb_z, if0_46, if0_79, \
                         if1_36, if1_61, kf_51, kf_69, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_y[k] * ld_26[k];

        t_85[k] = f_9 * if0_79[k]
                  - f_10 * if1_61[k]
                  + pa_x[k] * kf_69[k];

        t_86[k] = f_7 * if0_46[k]
                  - f_8 * if1_36[k]
                  + pa_y[k] * kf_51[k];

        t_87[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pb_x, pb_z, if0_84, if1_66, kd_30, \
                         kf_71, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_28[k];

        t_89[k] = f_4 * if0_84[k]
                  - f_5 * if1_66[k]
                  + pa_x[k] * kf_71[k];

        t_90[k] = pb_z[k] * ld_28[k];

        t_91[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, if0_55, if1_42, kd_26, kf_51, \
                         kf_54, kf_56, kf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * kf_51[k];

        t_93[k] = pa_z[k] * kf_54[k];

        t_94[k] = f_3 * kd_26[k]
                  + pa_z[k] * kf_56[k];

        t_95[k] = f_14 * if0_55[k]
                  - f_15 * if1_42[k]
                  + pa_y[k] * kf_57[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_x, pa_y, if0_58, if0_94, if0_96, if1_45, if1_74, \
                         if1_76, kf_60, kf_72, kf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * if0_94[k]
                  - f_5 * if1_74[k]
                  + pa_x[k] * kf_72[k];

        t_97[k] = f_4 * if0_96[k]
                  - f_5 * if1_76[k]
                  + pa_x[k] * kf_73[k];

        t_98[k] = f_9 * if0_58[k]
                  - f_10 * if1_45[k]
                  + pa_y[k] * kf_60[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_x, pa_y, if0_61, if0_100, if0_102, if1_46, \
                         if1_80, if1_82, kf_63, kf_74, kf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_4 * if0_100[k]
                  - f_5 * if1_80[k]
                  + pa_x[k] * kf_74[k];

        t_100[k] = f_4 * if0_102[k]
                   - f_5 * if1_82[k]
                   + pa_x[k] * kf_75[k];

        t_101[k] = f_4 * if0_61[k]
                   - f_5 * if1_46[k]
                   + pa_y[k] * kf_63[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_x, pa_y, if0_106, if0_108, if1_86, \
                         if1_88, kd_28, kf_67, kf_69, kf_76, kf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * if0_106[k]
                   - f_5 * if1_86[k]
                   + pa_x[k] * kf_76[k];

        t_103[k] = f_4 * if0_108[k]
                   - f_5 * if1_88[k]
                   + pa_x[k] * kf_77[k];

        t_104[k] = f_3 * kd_28[k]
                   + pa_y[k] * kf_67[k];

        t_105[k] = pa_y[k] * kf_69[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_z, pb_x, pb_y, if0_61, if1_46, kd_31, \
                         kf_64, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_7 * if0_61[k]
                   - f_8 * if1_46[k]
                   + pa_z[k] * kf_64[k];

        t_107[k] = pb_y[k] * ld_30[k];

        t_108[k] = f_17 * kd_31[k]
                   + pb_x[k] * ld_32[k];

        t_109[k] = f_1 * lp0_12[k]
                   - f_2 * lp1_12[k]
                   + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_x, pa_z, pb_y, if0_119, if1_98, \
                         kd_32, kf_70, kf_79, kf_80, kf_83, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_y[k] * ld_32[k];

        t_111[k] = f_4 * if0_119[k]
                   - f_5 * if1_98[k]
                   + pa_x[k] * kf_79[k];

        t_112[k] = f_3 * kd_32[k]
                   + pa_x[k] * kf_80[k];

        t_113[k] = pa_x[k] * kf_83[k];

        t_114[k] = pa_z[k] * kf_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_x, kd_37, kd_40, kd_43, kd_46, \
                         kd_51, kf_88, kf_94, kf_100, kf_106, kf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * kd_37[k]
                   + pa_x[k] * kf_88[k];

        t_116[k] = f_3 * kd_40[k]
                   + pa_x[k] * kf_94[k];

        t_117[k] = f_3 * kd_43[k]
                   + pa_x[k] * kf_100[k];

        t_118[k] = f_3 * kd_46[k]
                   + pa_x[k] * kf_106[k];

        t_119[k] = f_3 * kd_51[k]
                   + pa_x[k] * kf_114[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, pb_x, kf_119, lp0_13, lp1_13, \
                         ld_33, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * kf_119[k];

        t_121[k] = f_1 * lp0_13[k]
                   - f_2 * lp1_13[k]
                   + pb_x[k] * ld_33[k];

        t_122[k] = pb_x[k] * ld_34[k];

        t_123[k] = pb_x[k] * ld_35[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_z, pb_y, pb_z, kd_33, kf_80, lp0_14, \
                         lp0_15, lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_0 * kd_33[k]
                   + f_1 * lp0_14[k]
                   - f_2 * lp1_14[k]
                   + pb_y[k] * ld_34[k];

        t_125[k] = pb_z[k] * ld_34[k];

        t_126[k] = f_1 * lp0_15[k]
                   - f_2 * lp1_15[k]
                   + pb_z[k] * ld_35[k];

        t_127[k] = pa_z[k] * kf_80[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_z, pb_x, kd_34, kf_83, kf_85, \
                         lp0_16, lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pa_z[k] * kf_83[k];

        t_129[k] = f_3 * kd_34[k]
                   + pa_z[k] * kf_85[k];

        t_130[k] = f_1 * lp0_16[k]
                   - f_2 * lp1_16[k]
                   + pb_x[k] * ld_36[k];

        t_131[k] = pb_x[k] * ld_37[k];

        t_132[k] = pb_x[k] * ld_38[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_y, pa_z, pb_y, if0_84, if0_96, if1_66, \
                         if1_76, kd_39, kf_86, kf_93, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_4 * if0_84[k]
                   - f_5 * if1_66[k]
                   + pa_z[k] * kf_86[k];

        t_134[k] = f_6 * kd_39[k]
                   + pb_y[k] * ld_38[k];

        t_135[k] = f_7 * if0_96[k]
                   - f_8 * if1_76[k]
                   + pa_y[k] * kf_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, if0_89, if1_69, kf_91, \
                         lp0_17, lp1_17, ld_39, ld_40, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_1 * lp0_17[k]
                   - f_2 * lp1_17[k]
                   + pb_x[k] * ld_39[k];

        t_137[k] = pb_x[k] * ld_40[k];

        t_138[k] = pb_x[k] * ld_41[k];

        t_139[k] = f_9 * if0_89[k]
                   - f_10 * if1_69[k]
                   + pa_z[k] * kf_91[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_y, pb_x, pb_y, if0_102, if1_82, kd_42, \
                         kf_99, lp0_18, lp1_18, ld_41, ld_42, ld_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_11 * kd_42[k]
                   + pb_y[k] * ld_41[k];

        t_141[k] = f_12 * if0_102[k]
                   - f_13 * if1_82[k]
                   + pa_y[k] * kf_99[k];

        t_142[k] = f_1 * lp0_18[k]
                   - f_2 * lp1_18[k]
                   + pb_x[k] * ld_42[k];

        t_143[k] = pb_x[k] * ld_43[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_y, pa_z, pb_x, pb_y, if0_94, if0_108, \
                         if1_74, if1_88, kd_45, kf_97, kf_105, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_x[k] * ld_44[k];

        t_145[k] = f_14 * if0_94[k]
                   - f_15 * if1_74[k]
                   + pa_z[k] * kf_97[k];

        t_146[k] = f_16 * kd_45[k]
                   + pb_y[k] * ld_44[k];

        t_147[k] = f_14 * if0_108[k]
                   - f_15 * if1_88[k]
                   + pa_y[k] * kf_105[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_z, pb_x, if0_100, if1_80, kf_103, \
                         lp0_19, lp1_19, ld_45, ld_46, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_1 * lp0_19[k]
                   - f_2 * lp1_19[k]
                   + pb_x[k] * ld_45[k];

        t_149[k] = pb_x[k] * ld_46[k];

        t_150[k] = pb_x[k] * ld_47[k];

        t_151[k] = f_12 * if0_100[k]
                   - f_13 * if1_80[k]
                   + pa_z[k] * kf_103[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_y, pb_x, pb_y, if0_112, if1_91, kd_48, \
                         kf_111, lp0_20, lp1_20, ld_47, ld_48, ld_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * kd_48[k]
                   + pb_y[k] * ld_47[k];

        t_153[k] = f_9 * if0_112[k]
                   - f_10 * if1_91[k]
                   + pa_y[k] * kf_111[k];

        t_154[k] = f_1 * lp0_20[k]
                   - f_2 * lp1_20[k]
                   + pb_x[k] * ld_48[k];

        t_155[k] = pb_x[k] * ld_49[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pa_y, pa_z, pb_x, pb_y, if0_106, if0_119, \
                         if1_86, if1_98, kd_50, kf_109, kf_113, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_x[k] * ld_50[k];

        t_157[k] = f_7 * if0_106[k]
                   - f_8 * if1_86[k]
                   + pa_z[k] * kf_109[k];

        t_158[k] = f_17 * kd_50[k]
                   + pb_y[k] * ld_50[k];

        t_159[k] = f_4 * if0_119[k]
                   - f_5 * if1_98[k]
                   + pa_y[k] * kf_113[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_y, pb_x, kd_52, kf_117, kf_119, \
                         lp0_21, lp1_21, ld_51, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_3 * kd_52[k]
                   + pa_y[k] * kf_117[k];

        t_161[k] = pa_y[k] * kf_119[k];

        t_162[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_163[k] = pb_x[k] * ld_52[k];

        t_164[k] = pb_x[k] * ld_53[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pb_y, pb_z, kd_53, lp0_22, lp0_23, lp1_22, \
                         lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_166[k] = pb_y[k] * ld_53[k];

        t_167[k] = f_0 * kd_53[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_1 = buffer.data(if0 + 1);
    const auto *if0_2 = buffer.data(if0 + 2);
    const auto *if0_3 = buffer.data(if0 + 3);
    const auto *if0_5 = buffer.data(if0 + 5);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_9 = buffer.data(if0 + 9);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_12 = buffer.data(if0 + 12);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_15 = buffer.data(if0 + 15);
    const auto *if0_17 = buffer.data(if0 + 17);
    const auto *if0_18 = buffer.data(if0 + 18);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_24 = buffer.data(if0 + 24);
    const auto *if0_25 = buffer.data(if0 + 25);
    const auto *if0_27 = buffer.data(if0 + 27);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_30 = buffer.data(if0 + 30);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_33 = buffer.data(if0 + 33);
    const auto *if0_34 = buffer.data(if0 + 34);
    const auto *if0_35 = buffer.data(if0 + 35);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_72 = buffer.data(if1 + 72);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_80 = buffer.data(if1 + 80);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_92 = buffer.data(kf + 92);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_6, kf_6, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_6[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_6[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_5, if1_11, kf_11, lp0_3, lp1_3, \
                         ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_5[k]
                 - f_7 * if1_11[k]
                 + pa_x[k] * kf_11[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_10, kf_7, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_7[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_10[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_1, if0_8, if1_6, \
                         if1_19, kf_8, kf_19, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_8[k]
                  - f_7 * if1_19[k]
                  + pa_x[k] * kf_19[k];

        t_18[k] = f_8 * if0_1[k]
                  - f_9 * if1_6[k]
                  + pa_y[k] * kf_8[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_11, if1_23, kd_12, \
                         kf_23, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_12[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_11[k]
                  - f_12 * if1_23[k]
                  + pa_x[k] * kf_23[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_2, if1_7, kd_16, kf_14, \
                         lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_2[k]
                  - f_9 * if1_7[k]
                  + pa_z[k] * kf_14[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_16[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_3, if0_14, if1_8, \
                         if1_31, kf_20, kf_31, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_14[k]
                  - f_12 * if1_31[k]
                  + pa_x[k] * kf_31[k];

        t_30[k] = f_13 * if0_3[k]
                  - f_14 * if1_8[k]
                  + pa_y[k] * kf_20[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_17, if1_35, kd_18, \
                         kf_35, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_18[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_17[k]
                  - f_14 * if1_35[k]
                  + pa_x[k] * kf_35[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_6, if1_14, kd_22, \
                         kf_26, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_6[k]
                  - f_14 * if1_14[k]
                  + pa_z[k] * kf_26[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_22[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_9, if0_20, \
                         if1_20, if1_43, kf_32, kf_43, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_20[k]
                  - f_14 * if1_43[k]
                  + pa_x[k] * kf_43[k];

        t_42[k] = f_11 * if0_9[k]
                  - f_12 * if1_20[k]
                  + pa_y[k] * kf_32[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_21, if1_45, kd_24, \
                         kf_47, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_24[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_21[k]
                  - f_9 * if1_45[k]
                  + pa_x[k] * kf_47[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_12, if1_26, kd_28, \
                         kf_38, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_12[k]
                  - f_12 * if1_26[k]
                  + pa_z[k] * kf_38[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_28[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_15, if0_22, \
                         if1_32, if1_47, kf_44, kf_55, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_22[k]
                  - f_9 * if1_47[k]
                  + pa_x[k] * kf_55[k];

        t_54[k] = f_6 * if0_15[k]
                  - f_7 * if1_32[k]
                  + pa_y[k] * kf_44[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_23, if1_51, kd_29, \
                         kf_57, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_29[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_23[k]
                  - f_4 * if1_51[k]
                  + pa_x[k] * kf_57[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_18, if1_38, kd_30, \
                         kf_50, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_18[k]
                  - f_7 * if1_38[k]
                  + pa_z[k] * kf_50[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_35, if1_80, kf_59, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_35[k]
                  - f_4 * if1_80[k]
                  + pa_x[k] * kf_59[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_32, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_32[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_23, if1_51, kf_66, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_23[k]
                  - f_4 * if1_51[k]
                  + pa_z[k] * kf_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_27, if1_60, kd_37, \
                         kf_72, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_37[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_27[k]
                  - f_7 * if1_60[k]
                  + pa_y[k] * kf_72[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_24, if0_30, \
                         if1_54, if1_66, kd_40, kf_70, kf_78, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_24[k]
                  - f_9 * if1_54[k]
                  + pa_z[k] * kf_70[k];

        t_82[k] = f_10 * kd_40[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_30[k]
                  - f_12 * if1_66[k]
                  + pa_y[k] * kf_78[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_25, if1_58, kf_76, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_25[k]
                  - f_14 * if1_58[k]
                  + pa_z[k] * kf_76[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_33, if1_72, kd_43, \
                         kf_84, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_43[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_33[k]
                  - f_14 * if1_72[k]
                  + pa_y[k] * kf_84[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_28, if0_34, \
                         if1_64, if1_74, kd_46, kf_82, kf_90, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_28[k]
                  - f_12 * if1_64[k]
                  + pa_z[k] * kf_82[k];

        t_94[k] = f_16 * kd_46[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_34[k]
                  - f_9 * if1_74[k]
                  + pa_y[k] * kf_90[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_31, if1_70, kf_88, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_31[k]
                  - f_7 * if1_70[k]
                  + pa_z[k] * kf_88[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_35, if1_80, kd_47, \
                         kf_92, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_47[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_35[k]
                   - f_4 * if1_80[k]
                   + pa_y[k] * kf_92[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_50, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_50[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_19 = buffer.data(if0 + 19);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_26 = buffer.data(if0 + 26);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_32 = buffer.data(if0 + 32);
    const auto *if0_35 = buffer.data(if0 + 35);
    const auto *if0_38 = buffer.data(if0 + 38);
    const auto *if0_43 = buffer.data(if0 + 43);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_47 = buffer.data(if0 + 47);
    const auto *if0_51 = buffer.data(if0 + 51);
    const auto *if0_54 = buffer.data(if0 + 54);
    const auto *if0_58 = buffer.data(if0 + 58);
    const auto *if0_60 = buffer.data(if0 + 60);
    const auto *if0_64 = buffer.data(if0 + 64);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_70 = buffer.data(if0 + 70);
    const auto *if0_72 = buffer.data(if0 + 72);
    const auto *if0_74 = buffer.data(if0 + 74);
    const auto *if0_80 = buffer.data(if0 + 80);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_10 = buffer.data(if1 + 10);
    const auto *if1_13 = buffer.data(if1 + 13);
    const auto *if1_16 = buffer.data(if1 + 16);
    const auto *if1_21 = buffer.data(if1 + 21);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_25 = buffer.data(if1 + 25);
    const auto *if1_28 = buffer.data(if1 + 28);
    const auto *if1_33 = buffer.data(if1 + 33);
    const auto *if1_34 = buffer.data(if1 + 34);
    const auto *if1_37 = buffer.data(if1 + 37);
    const auto *if1_40 = buffer.data(if1 + 40);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_49 = buffer.data(if1 + 49);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_57 = buffer.data(if1 + 57);
    const auto *if1_62 = buffer.data(if1 + 62);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_68 = buffer.data(if1 + 68);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_76 = buffer.data(if1 + 76);
    const auto *if1_79 = buffer.data(if1 + 79);
    const auto *if1_86 = buffer.data(if1 + 86);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_95 = buffer.data(kf + 95);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_6, kf_6, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_6[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_6[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_11, if1_13, kf_12, lp0_3, lp1_3, \
                         ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_11[k]
                 - f_7 * if1_13[k]
                 + pa_x[k] * kf_12[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_10, kf_7, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_7[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_10[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_6, if0_19, if1_7, \
                         if1_21, kf_9, kf_20, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_19[k]
                  - f_7 * if1_21[k]
                  + pa_x[k] * kf_20[k];

        t_18[k] = f_8 * if0_6[k]
                  - f_9 * if1_7[k]
                  + pa_y[k] * kf_9[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_23, if1_25, kd_12, \
                         kf_24, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_12[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_23[k]
                  - f_12 * if1_25[k]
                  + pa_x[k] * kf_24[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_7, if1_8, kd_16, kf_15, \
                         lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_7[k]
                  - f_9 * if1_8[k]
                  + pa_z[k] * kf_15[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_16[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_8, if0_31, \
                         if1_10, if1_33, kf_21, kf_32, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_31[k]
                  - f_12 * if1_33[k]
                  + pa_x[k] * kf_32[k];

        t_30[k] = f_13 * if0_8[k]
                  - f_14 * if1_10[k]
                  + pa_y[k] * kf_21[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_35, if1_37, kd_18, \
                         kf_36, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_18[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_35[k]
                  - f_14 * if1_37[k]
                  + pa_x[k] * kf_36[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_14, if1_16, kd_22, \
                         kf_27, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_14[k]
                  - f_14 * if1_16[k]
                  + pa_z[k] * kf_27[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_22[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_20, if0_43, \
                         if1_22, if1_45, kf_33, kf_44, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_43[k]
                  - f_14 * if1_45[k]
                  + pa_x[k] * kf_44[k];

        t_42[k] = f_11 * if0_20[k]
                  - f_12 * if1_22[k]
                  + pa_y[k] * kf_33[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_45, if1_47, kd_24, \
                         kf_48, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_24[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_45[k]
                  - f_9 * if1_47[k]
                  + pa_x[k] * kf_48[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_26, if1_28, kd_28, \
                         kf_39, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_26[k]
                  - f_12 * if1_28[k]
                  + pa_z[k] * kf_39[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_28[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_32, if0_47, \
                         if1_34, if1_49, kf_45, kf_56, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_47[k]
                  - f_9 * if1_49[k]
                  + pa_x[k] * kf_56[k];

        t_54[k] = f_6 * if0_32[k]
                  - f_7 * if1_34[k]
                  + pa_y[k] * kf_45[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_51, if1_54, kd_29, \
                         kf_58, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_29[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_51[k]
                  - f_4 * if1_54[k]
                  + pa_x[k] * kf_58[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_38, if1_40, kd_30, \
                         kf_51, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_38[k]
                  - f_7 * if1_40[k]
                  + pa_z[k] * kf_51[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_80, if1_86, kf_60, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_80[k]
                  - f_4 * if1_86[k]
                  + pa_x[k] * kf_60[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_32, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_32[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_51, if1_54, kf_67, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_51[k]
                  - f_4 * if1_54[k]
                  + pa_z[k] * kf_67[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_60, if1_64, kd_37, \
                         kf_74, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_37[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_60[k]
                  - f_7 * if1_64[k]
                  + pa_y[k] * kf_74[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_54, if0_66, \
                         if1_57, if1_70, kd_40, kf_72, kf_80, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_54[k]
                  - f_9 * if1_57[k]
                  + pa_z[k] * kf_72[k];

        t_82[k] = f_10 * kd_40[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_66[k]
                  - f_12 * if1_70[k]
                  + pa_y[k] * kf_80[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_58, if1_62, kf_78, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_58[k]
                  - f_14 * if1_62[k]
                  + pa_z[k] * kf_78[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_72, if1_76, kd_43, \
                         kf_86, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_43[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_72[k]
                  - f_14 * if1_76[k]
                  + pa_y[k] * kf_86[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_64, if0_74, \
                         if1_68, if1_79, kd_46, kf_84, kf_92, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_64[k]
                  - f_12 * if1_68[k]
                  + pa_z[k] * kf_84[k];

        t_94[k] = f_16 * kd_46[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_74[k]
                  - f_9 * if1_79[k]
                  + pa_y[k] * kf_92[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_70, if1_74, kf_90, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_70[k]
                  - f_7 * if1_74[k]
                  + pa_z[k] * kf_90[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_80, if1_86, kd_47, \
                         kf_95, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_47[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_80[k]
                   - f_4 * if1_86[k]
                   + pa_y[k] * kf_95[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_50, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_50[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_10 = buffer.data(if0 + 10);
    const auto *if0_13 = buffer.data(if0 + 13);
    const auto *if0_16 = buffer.data(if0 + 16);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_25 = buffer.data(if0 + 25);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_33 = buffer.data(if0 + 33);
    const auto *if0_34 = buffer.data(if0 + 34);
    const auto *if0_37 = buffer.data(if0 + 37);
    const auto *if0_40 = buffer.data(if0 + 40);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_47 = buffer.data(if0 + 47);
    const auto *if0_49 = buffer.data(if0 + 49);
    const auto *if0_54 = buffer.data(if0 + 54);
    const auto *if0_57 = buffer.data(if0 + 57);
    const auto *if0_62 = buffer.data(if0 + 62);
    const auto *if0_64 = buffer.data(if0 + 64);
    const auto *if0_68 = buffer.data(if0 + 68);
    const auto *if0_70 = buffer.data(if0 + 70);
    const auto *if0_74 = buffer.data(if0 + 74);
    const auto *if0_76 = buffer.data(if0 + 76);
    const auto *if0_79 = buffer.data(if0 + 79);
    const auto *if0_86 = buffer.data(if0 + 86);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_72 = buffer.data(if1 + 72);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_80 = buffer.data(if1 + 80);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_6, kf_6, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_6[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_6[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_13, if1_11, kf_11, lp0_3, lp1_3, \
                         ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_13[k]
                 - f_7 * if1_11[k]
                 + pa_x[k] * kf_11[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_10, kf_7, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_7[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_10[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_7, if0_21, if1_6, \
                         if1_19, kf_8, kf_19, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_21[k]
                  - f_7 * if1_19[k]
                  + pa_x[k] * kf_19[k];

        t_18[k] = f_8 * if0_7[k]
                  - f_9 * if1_6[k]
                  + pa_y[k] * kf_8[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_25, if1_23, kd_12, \
                         kf_23, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_12[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_25[k]
                  - f_12 * if1_23[k]
                  + pa_x[k] * kf_23[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_8, if1_7, kd_16, kf_14, \
                         lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_8[k]
                  - f_9 * if1_7[k]
                  + pa_z[k] * kf_14[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_16[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_10, if0_33, \
                         if1_8, if1_31, kf_20, kf_31, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_33[k]
                  - f_12 * if1_31[k]
                  + pa_x[k] * kf_31[k];

        t_30[k] = f_13 * if0_10[k]
                  - f_14 * if1_8[k]
                  + pa_y[k] * kf_20[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_37, if1_35, kd_18, \
                         kf_35, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_18[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_37[k]
                  - f_14 * if1_35[k]
                  + pa_x[k] * kf_35[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_16, if1_14, kd_22, \
                         kf_26, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_16[k]
                  - f_14 * if1_14[k]
                  + pa_z[k] * kf_26[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_22[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_22, if0_45, \
                         if1_20, if1_43, kf_32, kf_43, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_45[k]
                  - f_14 * if1_43[k]
                  + pa_x[k] * kf_43[k];

        t_42[k] = f_11 * if0_22[k]
                  - f_12 * if1_20[k]
                  + pa_y[k] * kf_32[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_47, if1_45, kd_24, \
                         kf_47, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_24[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_47[k]
                  - f_9 * if1_45[k]
                  + pa_x[k] * kf_47[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_28, if1_26, kd_28, \
                         kf_38, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_28[k]
                  - f_12 * if1_26[k]
                  + pa_z[k] * kf_38[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_28[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_34, if0_49, \
                         if1_32, if1_47, kf_44, kf_55, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_49[k]
                  - f_9 * if1_47[k]
                  + pa_x[k] * kf_55[k];

        t_54[k] = f_6 * if0_34[k]
                  - f_7 * if1_32[k]
                  + pa_y[k] * kf_44[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_54, if1_51, kd_29, \
                         kf_56, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_29[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_54[k]
                  - f_4 * if1_51[k]
                  + pa_x[k] * kf_56[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_40, if1_38, kd_30, \
                         kf_50, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_40[k]
                  - f_7 * if1_38[k]
                  + pa_z[k] * kf_50[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_86, if1_80, kf_57, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_86[k]
                  - f_4 * if1_80[k]
                  + pa_x[k] * kf_57[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_32, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_32[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_54, if1_51, kf_64, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_54[k]
                  - f_4 * if1_51[k]
                  + pa_z[k] * kf_64[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_64, if1_60, kd_37, \
                         kf_70, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_37[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_64[k]
                  - f_7 * if1_60[k]
                  + pa_y[k] * kf_70[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_57, if0_70, \
                         if1_54, if1_66, kd_40, kf_68, kf_76, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_57[k]
                  - f_9 * if1_54[k]
                  + pa_z[k] * kf_68[k];

        t_82[k] = f_10 * kd_40[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_70[k]
                  - f_12 * if1_66[k]
                  + pa_y[k] * kf_76[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_62, if1_58, kf_74, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_62[k]
                  - f_14 * if1_58[k]
                  + pa_z[k] * kf_74[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_76, if1_72, kd_43, \
                         kf_82, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_43[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_76[k]
                  - f_14 * if1_72[k]
                  + pa_y[k] * kf_82[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_68, if0_79, \
                         if1_64, if1_74, kd_46, kf_80, kf_88, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_68[k]
                  - f_12 * if1_64[k]
                  + pa_z[k] * kf_80[k];

        t_94[k] = f_16 * kd_46[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_79[k]
                  - f_9 * if1_74[k]
                  + pa_y[k] * kf_88[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_74, if1_70, kf_86, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_74[k]
                  - f_7 * if1_70[k]
                  + pa_z[k] * kf_86[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_86, if1_80, kd_47, \
                         kf_89, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_47[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_86[k]
                   - f_4 * if1_80[k]
                   + pa_y[k] * kf_89[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_50, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_50[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_19 = buffer.data(if0 + 19);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_26 = buffer.data(if0 + 26);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_32 = buffer.data(if0 + 32);
    const auto *if0_35 = buffer.data(if0 + 35);
    const auto *if0_38 = buffer.data(if0 + 38);
    const auto *if0_43 = buffer.data(if0 + 43);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_47 = buffer.data(if0 + 47);
    const auto *if0_51 = buffer.data(if0 + 51);
    const auto *if0_54 = buffer.data(if0 + 54);
    const auto *if0_58 = buffer.data(if0 + 58);
    const auto *if0_60 = buffer.data(if0 + 60);
    const auto *if0_64 = buffer.data(if0 + 64);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_70 = buffer.data(if0 + 70);
    const auto *if0_72 = buffer.data(if0 + 72);
    const auto *if0_74 = buffer.data(if0 + 74);
    const auto *if0_80 = buffer.data(if0 + 80);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_72 = buffer.data(if1 + 72);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_80 = buffer.data(if1 + 80);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_92 = buffer.data(kf + 92);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_6, kf_6, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_6[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_6[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_11, if1_11, kf_11, lp0_3, lp1_3, \
                         ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_11[k]
                 - f_7 * if1_11[k]
                 + pa_x[k] * kf_11[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_10, kf_7, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_7[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_10[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_6, if0_19, if1_6, \
                         if1_19, kf_8, kf_19, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_19[k]
                  - f_7 * if1_19[k]
                  + pa_x[k] * kf_19[k];

        t_18[k] = f_8 * if0_6[k]
                  - f_9 * if1_6[k]
                  + pa_y[k] * kf_8[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_23, if1_23, kd_12, \
                         kf_23, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_12[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_23[k]
                  - f_12 * if1_23[k]
                  + pa_x[k] * kf_23[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_7, if1_7, kd_16, kf_14, \
                         lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_7[k]
                  - f_9 * if1_7[k]
                  + pa_z[k] * kf_14[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_16[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_8, if0_31, if1_8, \
                         if1_31, kf_20, kf_31, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_31[k]
                  - f_12 * if1_31[k]
                  + pa_x[k] * kf_31[k];

        t_30[k] = f_13 * if0_8[k]
                  - f_14 * if1_8[k]
                  + pa_y[k] * kf_20[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_35, if1_35, kd_18, \
                         kf_35, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_18[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_35[k]
                  - f_14 * if1_35[k]
                  + pa_x[k] * kf_35[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_14, if1_14, kd_22, \
                         kf_26, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_14[k]
                  - f_14 * if1_14[k]
                  + pa_z[k] * kf_26[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_22[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_20, if0_43, \
                         if1_20, if1_43, kf_32, kf_43, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_43[k]
                  - f_14 * if1_43[k]
                  + pa_x[k] * kf_43[k];

        t_42[k] = f_11 * if0_20[k]
                  - f_12 * if1_20[k]
                  + pa_y[k] * kf_32[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_45, if1_45, kd_24, \
                         kf_47, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_24[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_45[k]
                  - f_9 * if1_45[k]
                  + pa_x[k] * kf_47[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_26, if1_26, kd_28, \
                         kf_38, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_26[k]
                  - f_12 * if1_26[k]
                  + pa_z[k] * kf_38[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_28[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_32, if0_47, \
                         if1_32, if1_47, kf_44, kf_55, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_47[k]
                  - f_9 * if1_47[k]
                  + pa_x[k] * kf_55[k];

        t_54[k] = f_6 * if0_32[k]
                  - f_7 * if1_32[k]
                  + pa_y[k] * kf_44[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_51, if1_51, kd_29, \
                         kf_57, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_29[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_51[k]
                  - f_4 * if1_51[k]
                  + pa_x[k] * kf_57[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_38, if1_38, kd_30, \
                         kf_50, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_38[k]
                  - f_7 * if1_38[k]
                  + pa_z[k] * kf_50[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_80, if1_80, kf_59, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_80[k]
                  - f_4 * if1_80[k]
                  + pa_x[k] * kf_59[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_32, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_32[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_51, if1_51, kf_66, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_51[k]
                  - f_4 * if1_51[k]
                  + pa_z[k] * kf_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_60, if1_60, kd_37, \
                         kf_72, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_37[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_60[k]
                  - f_7 * if1_60[k]
                  + pa_y[k] * kf_72[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_54, if0_66, \
                         if1_54, if1_66, kd_40, kf_70, kf_78, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_54[k]
                  - f_9 * if1_54[k]
                  + pa_z[k] * kf_70[k];

        t_82[k] = f_10 * kd_40[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_66[k]
                  - f_12 * if1_66[k]
                  + pa_y[k] * kf_78[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_58, if1_58, kf_76, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_58[k]
                  - f_14 * if1_58[k]
                  + pa_z[k] * kf_76[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_72, if1_72, kd_43, \
                         kf_84, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_43[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_72[k]
                  - f_14 * if1_72[k]
                  + pa_y[k] * kf_84[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_64, if0_74, \
                         if1_64, if1_74, kd_46, kf_82, kf_90, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_64[k]
                  - f_12 * if1_64[k]
                  + pa_z[k] * kf_82[k];

        t_94[k] = f_16 * kd_46[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_74[k]
                  - f_9 * if1_74[k]
                  + pa_y[k] * kf_90[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_70, if1_70, kf_88, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_70[k]
                  - f_7 * if1_70[k]
                  + pa_z[k] * kf_88[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_80, if1_80, kd_47, \
                         kf_92, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_47[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_80[k]
                   - f_4 * if1_80[k]
                   + pa_y[k] * kf_92[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_50, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_50[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

auto
compute_prim_lf_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / alpha;
    const auto f_12 = 2.0 * beta / (alpha * p);
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 2.0 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_19 = buffer.data(if0 + 19);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_26 = buffer.data(if0 + 26);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_32 = buffer.data(if0 + 32);
    const auto *if0_35 = buffer.data(if0 + 35);
    const auto *if0_38 = buffer.data(if0 + 38);
    const auto *if0_43 = buffer.data(if0 + 43);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_47 = buffer.data(if0 + 47);
    const auto *if0_51 = buffer.data(if0 + 51);
    const auto *if0_54 = buffer.data(if0 + 54);
    const auto *if0_58 = buffer.data(if0 + 58);
    const auto *if0_60 = buffer.data(if0 + 60);
    const auto *if0_64 = buffer.data(if0 + 64);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_70 = buffer.data(if0 + 70);
    const auto *if0_72 = buffer.data(if0 + 72);
    const auto *if0_74 = buffer.data(if0 + 74);
    const auto *if0_80 = buffer.data(if0 + 80);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_72 = buffer.data(if1 + 72);
    const auto *if1_74 = buffer.data(if1 + 74);
    const auto *if1_80 = buffer.data(if1 + 80);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);

    const auto *lp0_0 = buffer.data(lp0 + 0);
    const auto *lp0_1 = buffer.data(lp0 + 1);
    const auto *lp0_2 = buffer.data(lp0 + 2);
    const auto *lp0_3 = buffer.data(lp0 + 3);
    const auto *lp0_4 = buffer.data(lp0 + 4);
    const auto *lp0_5 = buffer.data(lp0 + 5);
    const auto *lp0_6 = buffer.data(lp0 + 6);
    const auto *lp0_7 = buffer.data(lp0 + 7);
    const auto *lp0_8 = buffer.data(lp0 + 8);
    const auto *lp0_9 = buffer.data(lp0 + 9);
    const auto *lp0_10 = buffer.data(lp0 + 10);
    const auto *lp0_11 = buffer.data(lp0 + 11);
    const auto *lp0_12 = buffer.data(lp0 + 12);
    const auto *lp0_13 = buffer.data(lp0 + 13);
    const auto *lp0_14 = buffer.data(lp0 + 14);
    const auto *lp0_15 = buffer.data(lp0 + 15);
    const auto *lp0_16 = buffer.data(lp0 + 16);
    const auto *lp0_17 = buffer.data(lp0 + 17);
    const auto *lp0_18 = buffer.data(lp0 + 18);
    const auto *lp0_19 = buffer.data(lp0 + 19);
    const auto *lp0_20 = buffer.data(lp0 + 20);
    const auto *lp0_21 = buffer.data(lp0 + 21);
    const auto *lp0_22 = buffer.data(lp0 + 22);
    const auto *lp0_23 = buffer.data(lp0 + 23);

    const auto *lp1_0 = buffer.data(lp1 + 0);
    const auto *lp1_1 = buffer.data(lp1 + 1);
    const auto *lp1_2 = buffer.data(lp1 + 2);
    const auto *lp1_3 = buffer.data(lp1 + 3);
    const auto *lp1_4 = buffer.data(lp1 + 4);
    const auto *lp1_5 = buffer.data(lp1 + 5);
    const auto *lp1_6 = buffer.data(lp1 + 6);
    const auto *lp1_7 = buffer.data(lp1 + 7);
    const auto *lp1_8 = buffer.data(lp1 + 8);
    const auto *lp1_9 = buffer.data(lp1 + 9);
    const auto *lp1_10 = buffer.data(lp1 + 10);
    const auto *lp1_11 = buffer.data(lp1 + 11);
    const auto *lp1_12 = buffer.data(lp1 + 12);
    const auto *lp1_13 = buffer.data(lp1 + 13);
    const auto *lp1_14 = buffer.data(lp1 + 14);
    const auto *lp1_15 = buffer.data(lp1 + 15);
    const auto *lp1_16 = buffer.data(lp1 + 16);
    const auto *lp1_17 = buffer.data(lp1 + 17);
    const auto *lp1_18 = buffer.data(lp1 + 18);
    const auto *lp1_19 = buffer.data(lp1 + 19);
    const auto *lp1_20 = buffer.data(lp1 + 20);
    const auto *lp1_21 = buffer.data(lp1 + 21);
    const auto *lp1_22 = buffer.data(lp1 + 22);
    const auto *lp1_23 = buffer.data(lp1 + 23);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kd_0, lp0_0, lp0_1, lp1_0, \
                         lp1_1, ld_0, ld_1, ld_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_0[k]
                 + f_1 * lp0_0[k]
                 - f_2 * lp1_0[k]
                 + pb_x[k] * ld_0[k];

        t_1[k] = pb_y[k] * ld_0[k];

        t_2[k] = pb_z[k] * ld_0[k];

        t_3[k] = f_1 * lp0_1[k]
                 - f_2 * lp1_1[k]
                 + pb_y[k] * ld_1[k];

        t_4[k] = pb_y[k] * ld_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, if0_0, if1_0, kd_6, kf_6, \
                         lp0_2, lp1_2, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * lp0_2[k]
                 - f_2 * lp1_2[k]
                 + pb_z[k] * ld_2[k];

        t_6[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pa_y[k] * kf_6[k];

        t_7[k] = pb_z[k] * ld_3[k];

        t_8[k] = f_5 * kd_6[k]
                 + pb_x[k] * ld_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, if0_11, if1_11, kf_11, lp0_3, lp1_3, \
                         ld_4, ld_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * if0_11[k]
                 - f_7 * if1_11[k]
                 + pa_x[k] * kf_11[k];

        t_10[k] = pb_z[k] * ld_4[k];

        t_11[k] = f_1 * lp0_3[k]
                  - f_2 * lp1_3[k]
                  + pb_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, if0_0, if1_0, kd_10, kf_7, \
                         lp0_4, lp1_4, ld_6, ld_7, ld_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_0[k]
                  - f_4 * if1_0[k]
                  + pa_z[k] * kf_7[k];

        t_13[k] = pb_y[k] * ld_6[k];

        t_14[k] = f_5 * kd_10[k]
                  + pb_x[k] * ld_8[k];

        t_15[k] = f_1 * lp0_4[k]
                  - f_2 * lp1_4[k]
                  + pb_y[k] * ld_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, if0_6, if0_19, if1_6, \
                         if1_19, kf_8, kf_19, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ld_8[k];

        t_17[k] = f_6 * if0_19[k]
                  - f_7 * if1_19[k]
                  + pa_x[k] * kf_19[k];

        t_18[k] = f_8 * if0_6[k]
                  - f_9 * if1_6[k]
                  + pa_y[k] * kf_8[k];

        t_19[k] = pb_z[k] * ld_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, if0_23, if1_23, kd_12, \
                         kf_23, lp0_5, lp1_5, ld_10, ld_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * kd_12[k]
                  + pb_x[k] * ld_10[k];

        t_21[k] = f_11 * if0_23[k]
                  - f_12 * if1_23[k]
                  + pa_x[k] * kf_23[k];

        t_22[k] = pb_z[k] * ld_10[k];

        t_23[k] = f_1 * lp0_5[k]
                  - f_2 * lp1_5[k]
                  + pb_z[k] * ld_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, if0_7, if1_7, kd_16, kf_14, \
                         lp0_6, lp1_6, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * if0_7[k]
                  - f_9 * if1_7[k]
                  + pa_z[k] * kf_14[k];

        t_25[k] = pb_y[k] * ld_12[k];

        t_26[k] = f_10 * kd_16[k]
                  + pb_x[k] * ld_14[k];

        t_27[k] = f_1 * lp0_6[k]
                  - f_2 * lp1_6[k]
                  + pb_y[k] * ld_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, if0_8, if0_31, if1_8, \
                         if1_31, kf_20, kf_31, ld_14, ld_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ld_14[k];

        t_29[k] = f_11 * if0_31[k]
                  - f_12 * if1_31[k]
                  + pa_x[k] * kf_31[k];

        t_30[k] = f_13 * if0_8[k]
                  - f_14 * if1_8[k]
                  + pa_y[k] * kf_20[k];

        t_31[k] = pb_z[k] * ld_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, if0_35, if1_35, kd_18, \
                         kf_35, lp0_7, lp1_7, ld_16, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_15 * kd_18[k]
                  + pb_x[k] * ld_16[k];

        t_33[k] = f_13 * if0_35[k]
                  - f_14 * if1_35[k]
                  + pa_x[k] * kf_35[k];

        t_34[k] = pb_z[k] * ld_16[k];

        t_35[k] = f_1 * lp0_7[k]
                  - f_2 * lp1_7[k]
                  + pb_z[k] * ld_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, if0_14, if1_14, kd_22, \
                         kf_26, lp0_8, lp1_8, ld_18, ld_19, ld_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * if0_14[k]
                  - f_14 * if1_14[k]
                  + pa_z[k] * kf_26[k];

        t_37[k] = pb_y[k] * ld_18[k];

        t_38[k] = f_15 * kd_22[k]
                  + pb_x[k] * ld_20[k];

        t_39[k] = f_1 * lp0_8[k]
                  - f_2 * lp1_8[k]
                  + pb_y[k] * ld_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, if0_20, if0_43, \
                         if1_20, if1_43, kf_32, kf_43, ld_20, ld_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ld_20[k];

        t_41[k] = f_13 * if0_43[k]
                  - f_14 * if1_43[k]
                  + pa_x[k] * kf_43[k];

        t_42[k] = f_11 * if0_20[k]
                  - f_12 * if1_20[k]
                  + pa_y[k] * kf_32[k];

        t_43[k] = pb_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, if0_45, if1_45, kd_24, \
                         kf_47, lp0_9, lp1_9, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_16 * kd_24[k]
                  + pb_x[k] * ld_22[k];

        t_45[k] = f_8 * if0_45[k]
                  - f_9 * if1_45[k]
                  + pa_x[k] * kf_47[k];

        t_46[k] = pb_z[k] * ld_22[k];

        t_47[k] = f_1 * lp0_9[k]
                  - f_2 * lp1_9[k]
                  + pb_z[k] * ld_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, if0_26, if1_26, kd_28, \
                         kf_38, lp0_10, lp1_10, ld_24, ld_25, ld_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if0_26[k]
                  - f_12 * if1_26[k]
                  + pa_z[k] * kf_38[k];

        t_49[k] = pb_y[k] * ld_24[k];

        t_50[k] = f_16 * kd_28[k]
                  + pb_x[k] * ld_26[k];

        t_51[k] = f_1 * lp0_10[k]
                  - f_2 * lp1_10[k]
                  + pb_y[k] * ld_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, pb_y, pb_z, if0_32, if0_47, \
                         if1_32, if1_47, kf_44, kf_55, ld_26, ld_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ld_26[k];

        t_53[k] = f_8 * if0_47[k]
                  - f_9 * if1_47[k]
                  + pa_x[k] * kf_55[k];

        t_54[k] = f_6 * if0_32[k]
                  - f_7 * if1_32[k]
                  + pa_y[k] * kf_44[k];

        t_55[k] = pb_z[k] * ld_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, if0_51, if1_51, kd_29, \
                         kf_56, lp0_11, lp1_11, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_17 * kd_29[k]
                  + pb_x[k] * ld_28[k];

        t_57[k] = f_3 * if0_51[k]
                  - f_4 * if1_51[k]
                  + pa_x[k] * kf_56[k];

        t_58[k] = pb_z[k] * ld_28[k];

        t_59[k] = f_1 * lp0_11[k]
                  - f_2 * lp1_11[k]
                  + pb_z[k] * ld_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, pb_y, if0_38, if1_38, kd_30, \
                         kf_50, lp0_12, lp1_12, ld_30, ld_31, ld_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_6 * if0_38[k]
                  - f_7 * if1_38[k]
                  + pa_z[k] * kf_50[k];

        t_61[k] = pb_y[k] * ld_30[k];

        t_62[k] = f_17 * kd_30[k]
                  + pb_x[k] * ld_32[k];

        t_63[k] = f_1 * lp0_12[k]
                  - f_2 * lp1_12[k]
                  + pb_y[k] * ld_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_y, if0_80, if1_80, kf_57, \
                         lp0_13, lp1_13, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * ld_32[k];

        t_65[k] = f_3 * if0_80[k]
                  - f_4 * if1_80[k]
                  + pa_x[k] * kf_57[k];

        t_66[k] = f_1 * lp0_13[k]
                  - f_2 * lp1_13[k]
                  + pb_x[k] * ld_33[k];

        t_67[k] = pb_x[k] * ld_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, kd_32, lp0_14, lp0_15, \
                         lp1_14, lp1_15, ld_34, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * ld_35[k];

        t_69[k] = f_0 * kd_32[k]
                  + f_1 * lp0_14[k]
                  - f_2 * lp1_14[k]
                  + pb_y[k] * ld_34[k];

        t_70[k] = pb_z[k] * ld_34[k];

        t_71[k] = f_1 * lp0_15[k]
                  - f_2 * lp1_15[k]
                  + pb_z[k] * ld_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, if0_51, if1_51, kf_64, lp0_16, \
                         lp1_16, ld_36, ld_37, ld_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * lp0_16[k]
                  - f_2 * lp1_16[k]
                  + pb_x[k] * ld_36[k];

        t_73[k] = pb_x[k] * ld_37[k];

        t_74[k] = pb_x[k] * ld_38[k];

        t_75[k] = f_3 * if0_51[k]
                  - f_4 * if1_51[k]
                  + pa_z[k] * kf_64[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, if0_60, if1_60, kd_37, \
                         kf_70, lp0_17, lp1_17, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * kd_37[k]
                  + pb_y[k] * ld_38[k];

        t_77[k] = f_6 * if0_60[k]
                  - f_7 * if1_60[k]
                  + pa_y[k] * kf_70[k];

        t_78[k] = f_1 * lp0_17[k]
                  - f_2 * lp1_17[k]
                  + pb_x[k] * ld_39[k];

        t_79[k] = pb_x[k] * ld_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, if0_54, if0_66, \
                         if1_54, if1_66, kd_40, kf_68, kf_76, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * ld_41[k];

        t_81[k] = f_8 * if0_54[k]
                  - f_9 * if1_54[k]
                  + pa_z[k] * kf_68[k];

        t_82[k] = f_10 * kd_40[k]
                  + pb_y[k] * ld_41[k];

        t_83[k] = f_11 * if0_66[k]
                  - f_12 * if1_66[k]
                  + pa_y[k] * kf_76[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_z, pb_x, if0_58, if1_58, kf_74, lp0_18, \
                         lp1_18, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * lp0_18[k]
                  - f_2 * lp1_18[k]
                  + pb_x[k] * ld_42[k];

        t_85[k] = pb_x[k] * ld_43[k];

        t_86[k] = pb_x[k] * ld_44[k];

        t_87[k] = f_13 * if0_58[k]
                  - f_14 * if1_58[k]
                  + pa_z[k] * kf_74[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pb_x, pb_y, if0_72, if1_72, kd_43, \
                         kf_82, lp0_19, lp1_19, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * kd_43[k]
                  + pb_y[k] * ld_44[k];

        t_89[k] = f_13 * if0_72[k]
                  - f_14 * if1_72[k]
                  + pa_y[k] * kf_82[k];

        t_90[k] = f_1 * lp0_19[k]
                  - f_2 * lp1_19[k]
                  + pb_x[k] * ld_45[k];

        t_91[k] = pb_x[k] * ld_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, if0_64, if0_74, \
                         if1_64, if1_74, kd_46, kf_80, kf_88, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * ld_47[k];

        t_93[k] = f_11 * if0_64[k]
                  - f_12 * if1_64[k]
                  + pa_z[k] * kf_80[k];

        t_94[k] = f_16 * kd_46[k]
                  + pb_y[k] * ld_47[k];

        t_95[k] = f_8 * if0_74[k]
                  - f_9 * if1_74[k]
                  + pa_y[k] * kf_88[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, if0_70, if1_70, kf_86, lp0_20, \
                         lp1_20, ld_48, ld_49, ld_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * lp0_20[k]
                  - f_2 * lp1_20[k]
                  + pb_x[k] * ld_48[k];

        t_97[k] = pb_x[k] * ld_49[k];

        t_98[k] = pb_x[k] * ld_50[k];

        t_99[k] = f_6 * if0_70[k]
                  - f_7 * if1_70[k]
                  + pa_z[k] * kf_86[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, if0_80, if1_80, kd_47, \
                         kf_89, lp0_21, lp1_21, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kd_47[k]
                   + pb_y[k] * ld_50[k];

        t_101[k] = f_3 * if0_80[k]
                   - f_4 * if1_80[k]
                   + pa_y[k] * kf_89[k];

        t_102[k] = f_1 * lp0_21[k]
                   - f_2 * lp1_21[k]
                   + pb_x[k] * ld_51[k];

        t_103[k] = pb_x[k] * ld_52[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, kd_50, lp0_22, lp0_23, \
                         lp1_22, lp1_23, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * ld_53[k];

        t_105[k] = f_1 * lp0_22[k]
                   - f_2 * lp1_22[k]
                   + pb_y[k] * ld_52[k];

        t_106[k] = pb_y[k] * ld_53[k];

        t_107[k] = f_0 * kd_50[k]
                   + f_1 * lp0_23[k]
                   - f_2 * lp1_23[k]
                   + pb_z[k] * ld_53[k];
    }
}

}  // namespace simdt2ceri
