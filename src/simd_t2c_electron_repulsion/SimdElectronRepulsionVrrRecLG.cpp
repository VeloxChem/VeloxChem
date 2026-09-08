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


#include "SimdElectronRepulsionVrrRecLG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_lg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ig0, const size_t ig1,
                                     const size_t kf, const size_t kg, const size_t ld0,
                                     const size_t ld1, const size_t lf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 3.5 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 3.0 / p;
    const auto f_12 = 2.5 / alpha;
    const auto f_13 = 2.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 2.0 / alpha;
    const auto f_19 = 2.0 * beta / (alpha * p);
    const auto f_20 = 1.5 / alpha;
    const auto f_21 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_15 = buffer.data(ig0 + 15);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_45 = buffer.data(ig0 + 45);
    const auto *ig0_48 = buffer.data(ig0 + 48);
    const auto *ig0_55 = buffer.data(ig0 + 55);
    const auto *ig0_75 = buffer.data(ig0 + 75);
    const auto *ig0_80 = buffer.data(ig0 + 80);
    const auto *ig0_89 = buffer.data(ig0 + 89);
    const auto *ig0_90 = buffer.data(ig0 + 90);
    const auto *ig0_93 = buffer.data(ig0 + 93);
    const auto *ig0_100 = buffer.data(ig0 + 100);
    const auto *ig0_108 = buffer.data(ig0 + 108);
    const auto *ig0_120 = buffer.data(ig0 + 120);
    const auto *ig0_125 = buffer.data(ig0 + 125);
    const auto *ig0_135 = buffer.data(ig0 + 135);
    const auto *ig0_140 = buffer.data(ig0 + 140);
    const auto *ig0_149 = buffer.data(ig0 + 149);
    const auto *ig0_150 = buffer.data(ig0 + 150);
    const auto *ig0_153 = buffer.data(ig0 + 153);
    const auto *ig0_160 = buffer.data(ig0 + 160);
    const auto *ig0_168 = buffer.data(ig0 + 168);
    const auto *ig0_180 = buffer.data(ig0 + 180);
    const auto *ig0_183 = buffer.data(ig0 + 183);
    const auto *ig0_185 = buffer.data(ig0 + 185);
    const auto *ig0_190 = buffer.data(ig0 + 190);
    const auto *ig0_192 = buffer.data(ig0 + 192);
    const auto *ig0_194 = buffer.data(ig0 + 194);
    const auto *ig0_195 = buffer.data(ig0 + 195);
    const auto *ig0_200 = buffer.data(ig0 + 200);
    const auto *ig0_210 = buffer.data(ig0 + 210);
    const auto *ig0_215 = buffer.data(ig0 + 215);
    const auto *ig0_224 = buffer.data(ig0 + 224);
    const auto *ig0_235 = buffer.data(ig0 + 235);
    const auto *ig0_265 = buffer.data(ig0 + 265);
    const auto *ig0_267 = buffer.data(ig0 + 267);
    const auto *ig0_269 = buffer.data(ig0 + 269);
    const auto *ig0_280 = buffer.data(ig0 + 280);
    const auto *ig0_282 = buffer.data(ig0 + 282);
    const auto *ig0_284 = buffer.data(ig0 + 284);
    const auto *ig0_314 = buffer.data(ig0 + 314);
    const auto *ig0_325 = buffer.data(ig0 + 325);
    const auto *ig0_340 = buffer.data(ig0 + 340);
    const auto *ig0_355 = buffer.data(ig0 + 355);
    const auto *ig0_357 = buffer.data(ig0 + 357);
    const auto *ig0_359 = buffer.data(ig0 + 359);
    const auto *ig0_370 = buffer.data(ig0 + 370);
    const auto *ig0_372 = buffer.data(ig0 + 372);
    const auto *ig0_374 = buffer.data(ig0 + 374);
    const auto *ig0_385 = buffer.data(ig0 + 385);
    const auto *ig0_387 = buffer.data(ig0 + 387);
    const auto *ig0_389 = buffer.data(ig0 + 389);
    const auto *ig0_404 = buffer.data(ig0 + 404);
    const auto *ig0_419 = buffer.data(ig0 + 419);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_15 = buffer.data(ig1 + 15);
    const auto *ig1_30 = buffer.data(ig1 + 30);
    const auto *ig1_45 = buffer.data(ig1 + 45);
    const auto *ig1_48 = buffer.data(ig1 + 48);
    const auto *ig1_55 = buffer.data(ig1 + 55);
    const auto *ig1_75 = buffer.data(ig1 + 75);
    const auto *ig1_80 = buffer.data(ig1 + 80);
    const auto *ig1_89 = buffer.data(ig1 + 89);
    const auto *ig1_90 = buffer.data(ig1 + 90);
    const auto *ig1_93 = buffer.data(ig1 + 93);
    const auto *ig1_100 = buffer.data(ig1 + 100);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_120 = buffer.data(ig1 + 120);
    const auto *ig1_125 = buffer.data(ig1 + 125);
    const auto *ig1_135 = buffer.data(ig1 + 135);
    const auto *ig1_140 = buffer.data(ig1 + 140);
    const auto *ig1_149 = buffer.data(ig1 + 149);
    const auto *ig1_150 = buffer.data(ig1 + 150);
    const auto *ig1_153 = buffer.data(ig1 + 153);
    const auto *ig1_160 = buffer.data(ig1 + 160);
    const auto *ig1_168 = buffer.data(ig1 + 168);
    const auto *ig1_180 = buffer.data(ig1 + 180);
    const auto *ig1_183 = buffer.data(ig1 + 183);
    const auto *ig1_185 = buffer.data(ig1 + 185);
    const auto *ig1_190 = buffer.data(ig1 + 190);
    const auto *ig1_192 = buffer.data(ig1 + 192);
    const auto *ig1_194 = buffer.data(ig1 + 194);
    const auto *ig1_195 = buffer.data(ig1 + 195);
    const auto *ig1_200 = buffer.data(ig1 + 200);
    const auto *ig1_210 = buffer.data(ig1 + 210);
    const auto *ig1_215 = buffer.data(ig1 + 215);
    const auto *ig1_224 = buffer.data(ig1 + 224);
    const auto *ig1_235 = buffer.data(ig1 + 235);
    const auto *ig1_265 = buffer.data(ig1 + 265);
    const auto *ig1_267 = buffer.data(ig1 + 267);
    const auto *ig1_269 = buffer.data(ig1 + 269);
    const auto *ig1_280 = buffer.data(ig1 + 280);
    const auto *ig1_282 = buffer.data(ig1 + 282);
    const auto *ig1_284 = buffer.data(ig1 + 284);
    const auto *ig1_314 = buffer.data(ig1 + 314);
    const auto *ig1_325 = buffer.data(ig1 + 325);
    const auto *ig1_340 = buffer.data(ig1 + 340);
    const auto *ig1_355 = buffer.data(ig1 + 355);
    const auto *ig1_357 = buffer.data(ig1 + 357);
    const auto *ig1_359 = buffer.data(ig1 + 359);
    const auto *ig1_370 = buffer.data(ig1 + 370);
    const auto *ig1_372 = buffer.data(ig1 + 372);
    const auto *ig1_374 = buffer.data(ig1 + 374);
    const auto *ig1_385 = buffer.data(ig1 + 385);
    const auto *ig1_387 = buffer.data(ig1 + 387);
    const auto *ig1_389 = buffer.data(ig1 + 389);
    const auto *ig1_404 = buffer.data(ig1 + 404);
    const auto *ig1_419 = buffer.data(ig1 + 419);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_282 = buffer.data(kf + 282);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_351 = buffer.data(kf + 351);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_3 = buffer.data(kg + 3);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_6 = buffer.data(kg + 6);
    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_18 = buffer.data(kg + 18);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_45 = buffer.data(kg + 45);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_84 = buffer.data(kg + 84);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_93 = buffer.data(kg + 93);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_149 = buffer.data(kg + 149);
    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_153 = buffer.data(kg + 153);
    const auto *kg_155 = buffer.data(kg + 155);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_160 = buffer.data(kg + 160);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_164 = buffer.data(kg + 164);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_185 = buffer.data(kg + 185);
    const auto *kg_190 = buffer.data(kg + 190);
    const auto *kg_192 = buffer.data(kg + 192);
    const auto *kg_194 = buffer.data(kg + 194);
    const auto *kg_195 = buffer.data(kg + 195);
    const auto *kg_200 = buffer.data(kg + 200);
    const auto *kg_210 = buffer.data(kg + 210);
    const auto *kg_212 = buffer.data(kg + 212);
    const auto *kg_213 = buffer.data(kg + 213);
    const auto *kg_215 = buffer.data(kg + 215);
    const auto *kg_219 = buffer.data(kg + 219);
    const auto *kg_220 = buffer.data(kg + 220);
    const auto *kg_222 = buffer.data(kg + 222);
    const auto *kg_224 = buffer.data(kg + 224);
    const auto *kg_225 = buffer.data(kg + 225);
    const auto *kg_226 = buffer.data(kg + 226);
    const auto *kg_228 = buffer.data(kg + 228);
    const auto *kg_230 = buffer.data(kg + 230);
    const auto *kg_231 = buffer.data(kg + 231);
    const auto *kg_235 = buffer.data(kg + 235);
    const auto *kg_237 = buffer.data(kg + 237);
    const auto *kg_239 = buffer.data(kg + 239);
    const auto *kg_243 = buffer.data(kg + 243);
    const auto *kg_255 = buffer.data(kg + 255);
    const auto *kg_258 = buffer.data(kg + 258);
    const auto *kg_260 = buffer.data(kg + 260);
    const auto *kg_265 = buffer.data(kg + 265);
    const auto *kg_267 = buffer.data(kg + 267);
    const auto *kg_269 = buffer.data(kg + 269);
    const auto *kg_270 = buffer.data(kg + 270);
    const auto *kg_273 = buffer.data(kg + 273);
    const auto *kg_275 = buffer.data(kg + 275);
    const auto *kg_280 = buffer.data(kg + 280);
    const auto *kg_282 = buffer.data(kg + 282);
    const auto *kg_284 = buffer.data(kg + 284);
    const auto *kg_285 = buffer.data(kg + 285);
    const auto *kg_290 = buffer.data(kg + 290);
    const auto *kg_300 = buffer.data(kg + 300);
    const auto *kg_302 = buffer.data(kg + 302);
    const auto *kg_303 = buffer.data(kg + 303);
    const auto *kg_305 = buffer.data(kg + 305);
    const auto *kg_309 = buffer.data(kg + 309);
    const auto *kg_310 = buffer.data(kg + 310);
    const auto *kg_312 = buffer.data(kg + 312);
    const auto *kg_314 = buffer.data(kg + 314);
    const auto *kg_315 = buffer.data(kg + 315);
    const auto *kg_316 = buffer.data(kg + 316);
    const auto *kg_318 = buffer.data(kg + 318);
    const auto *kg_321 = buffer.data(kg + 321);
    const auto *kg_325 = buffer.data(kg + 325);
    const auto *kg_355 = buffer.data(kg + 355);
    const auto *kg_357 = buffer.data(kg + 357);
    const auto *kg_359 = buffer.data(kg + 359);
    const auto *kg_370 = buffer.data(kg + 370);
    const auto *kg_372 = buffer.data(kg + 372);
    const auto *kg_374 = buffer.data(kg + 374);
    const auto *kg_385 = buffer.data(kg + 385);
    const auto *kg_387 = buffer.data(kg + 387);
    const auto *kg_389 = buffer.data(kg + 389);
    const auto *kg_405 = buffer.data(kg + 405);
    const auto *kg_407 = buffer.data(kg + 407);
    const auto *kg_410 = buffer.data(kg + 410);
    const auto *kg_414 = buffer.data(kg + 414);
    const auto *kg_419 = buffer.data(kg + 419);
    const auto *kg_420 = buffer.data(kg + 420);
    const auto *kg_421 = buffer.data(kg + 421);
    const auto *kg_423 = buffer.data(kg + 423);
    const auto *kg_425 = buffer.data(kg + 425);
    const auto *kg_430 = buffer.data(kg + 430);
    const auto *kg_432 = buffer.data(kg + 432);
    const auto *kg_433 = buffer.data(kg + 433);
    const auto *kg_434 = buffer.data(kg + 434);
    const auto *kg_440 = buffer.data(kg + 440);
    const auto *kg_445 = buffer.data(kg + 445);
    const auto *kg_446 = buffer.data(kg + 446);
    const auto *kg_447 = buffer.data(kg + 447);
    const auto *kg_448 = buffer.data(kg + 448);
    const auto *kg_449 = buffer.data(kg + 449);
    const auto *kg_450 = buffer.data(kg + 450);
    const auto *kg_453 = buffer.data(kg + 453);
    const auto *kg_455 = buffer.data(kg + 455);
    const auto *kg_460 = buffer.data(kg + 460);
    const auto *kg_461 = buffer.data(kg + 461);
    const auto *kg_462 = buffer.data(kg + 462);
    const auto *kg_463 = buffer.data(kg + 463);
    const auto *kg_464 = buffer.data(kg + 464);
    const auto *kg_465 = buffer.data(kg + 465);
    const auto *kg_468 = buffer.data(kg + 468);
    const auto *kg_470 = buffer.data(kg + 470);
    const auto *kg_475 = buffer.data(kg + 475);
    const auto *kg_476 = buffer.data(kg + 476);
    const auto *kg_477 = buffer.data(kg + 477);
    const auto *kg_478 = buffer.data(kg + 478);
    const auto *kg_479 = buffer.data(kg + 479);
    const auto *kg_480 = buffer.data(kg + 480);
    const auto *kg_483 = buffer.data(kg + 483);
    const auto *kg_485 = buffer.data(kg + 485);
    const auto *kg_490 = buffer.data(kg + 490);
    const auto *kg_491 = buffer.data(kg + 491);
    const auto *kg_492 = buffer.data(kg + 492);
    const auto *kg_493 = buffer.data(kg + 493);
    const auto *kg_494 = buffer.data(kg + 494);
    const auto *kg_495 = buffer.data(kg + 495);
    const auto *kg_498 = buffer.data(kg + 498);
    const auto *kg_500 = buffer.data(kg + 500);
    const auto *kg_505 = buffer.data(kg + 505);
    const auto *kg_506 = buffer.data(kg + 506);
    const auto *kg_507 = buffer.data(kg + 507);
    const auto *kg_508 = buffer.data(kg + 508);
    const auto *kg_509 = buffer.data(kg + 509);
    const auto *kg_513 = buffer.data(kg + 513);
    const auto *kg_520 = buffer.data(kg + 520);
    const auto *kg_521 = buffer.data(kg + 521);
    const auto *kg_522 = buffer.data(kg + 522);
    const auto *kg_523 = buffer.data(kg + 523);
    const auto *kg_524 = buffer.data(kg + 524);
    const auto *kg_525 = buffer.data(kg + 525);
    const auto *kg_527 = buffer.data(kg + 527);
    const auto *kg_528 = buffer.data(kg + 528);
    const auto *kg_530 = buffer.data(kg + 530);
    const auto *kg_535 = buffer.data(kg + 535);
    const auto *kg_536 = buffer.data(kg + 536);
    const auto *kg_537 = buffer.data(kg + 537);
    const auto *kg_539 = buffer.data(kg + 539);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_54 = buffer.data(ld0 + 54);
    const auto *ld0_57 = buffer.data(ld0 + 57);
    const auto *ld0_59 = buffer.data(ld0 + 59);
    const auto *ld0_60 = buffer.data(ld0 + 60);
    const auto *ld0_63 = buffer.data(ld0 + 63);
    const auto *ld0_65 = buffer.data(ld0 + 65);
    const auto *ld0_84 = buffer.data(ld0 + 84);
    const auto *ld0_87 = buffer.data(ld0 + 87);
    const auto *ld0_89 = buffer.data(ld0 + 89);
    const auto *ld0_90 = buffer.data(ld0 + 90);
    const auto *ld0_93 = buffer.data(ld0 + 93);
    const auto *ld0_95 = buffer.data(ld0 + 95);
    const auto *ld0_120 = buffer.data(ld0 + 120);
    const auto *ld0_123 = buffer.data(ld0 + 123);
    const auto *ld0_125 = buffer.data(ld0 + 125);
    const auto *ld0_126 = buffer.data(ld0 + 126);
    const auto *ld0_129 = buffer.data(ld0 + 129);
    const auto *ld0_131 = buffer.data(ld0 + 131);
    const auto *ld0_162 = buffer.data(ld0 + 162);
    const auto *ld0_165 = buffer.data(ld0 + 165);
    const auto *ld0_167 = buffer.data(ld0 + 167);
    const auto *ld0_216 = buffer.data(ld0 + 216);
    const auto *ld0_219 = buffer.data(ld0 + 219);
    const auto *ld0_221 = buffer.data(ld0 + 221);
    const auto *ld0_228 = buffer.data(ld0 + 228);
    const auto *ld0_231 = buffer.data(ld0 + 231);
    const auto *ld0_233 = buffer.data(ld0 + 233);
    const auto *ld0_234 = buffer.data(ld0 + 234);
    const auto *ld0_237 = buffer.data(ld0 + 237);
    const auto *ld0_239 = buffer.data(ld0 + 239);
    const auto *ld0_240 = buffer.data(ld0 + 240);
    const auto *ld0_243 = buffer.data(ld0 + 243);
    const auto *ld0_245 = buffer.data(ld0 + 245);
    const auto *ld0_246 = buffer.data(ld0 + 246);
    const auto *ld0_249 = buffer.data(ld0 + 249);
    const auto *ld0_251 = buffer.data(ld0 + 251);
    const auto *ld0_252 = buffer.data(ld0 + 252);
    const auto *ld0_255 = buffer.data(ld0 + 255);
    const auto *ld0_257 = buffer.data(ld0 + 257);
    const auto *ld0_264 = buffer.data(ld0 + 264);
    const auto *ld0_267 = buffer.data(ld0 + 267);
    const auto *ld0_269 = buffer.data(ld0 + 269);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_54 = buffer.data(ld1 + 54);
    const auto *ld1_57 = buffer.data(ld1 + 57);
    const auto *ld1_59 = buffer.data(ld1 + 59);
    const auto *ld1_60 = buffer.data(ld1 + 60);
    const auto *ld1_63 = buffer.data(ld1 + 63);
    const auto *ld1_65 = buffer.data(ld1 + 65);
    const auto *ld1_84 = buffer.data(ld1 + 84);
    const auto *ld1_87 = buffer.data(ld1 + 87);
    const auto *ld1_89 = buffer.data(ld1 + 89);
    const auto *ld1_90 = buffer.data(ld1 + 90);
    const auto *ld1_93 = buffer.data(ld1 + 93);
    const auto *ld1_95 = buffer.data(ld1 + 95);
    const auto *ld1_120 = buffer.data(ld1 + 120);
    const auto *ld1_123 = buffer.data(ld1 + 123);
    const auto *ld1_125 = buffer.data(ld1 + 125);
    const auto *ld1_126 = buffer.data(ld1 + 126);
    const auto *ld1_129 = buffer.data(ld1 + 129);
    const auto *ld1_131 = buffer.data(ld1 + 131);
    const auto *ld1_162 = buffer.data(ld1 + 162);
    const auto *ld1_165 = buffer.data(ld1 + 165);
    const auto *ld1_167 = buffer.data(ld1 + 167);
    const auto *ld1_216 = buffer.data(ld1 + 216);
    const auto *ld1_219 = buffer.data(ld1 + 219);
    const auto *ld1_221 = buffer.data(ld1 + 221);
    const auto *ld1_228 = buffer.data(ld1 + 228);
    const auto *ld1_231 = buffer.data(ld1 + 231);
    const auto *ld1_233 = buffer.data(ld1 + 233);
    const auto *ld1_234 = buffer.data(ld1 + 234);
    const auto *ld1_237 = buffer.data(ld1 + 237);
    const auto *ld1_239 = buffer.data(ld1 + 239);
    const auto *ld1_240 = buffer.data(ld1 + 240);
    const auto *ld1_243 = buffer.data(ld1 + 243);
    const auto *ld1_245 = buffer.data(ld1 + 245);
    const auto *ld1_246 = buffer.data(ld1 + 246);
    const auto *ld1_249 = buffer.data(ld1 + 249);
    const auto *ld1_251 = buffer.data(ld1 + 251);
    const auto *ld1_252 = buffer.data(ld1 + 252);
    const auto *ld1_255 = buffer.data(ld1 + 255);
    const auto *ld1_257 = buffer.data(ld1 + 257);
    const auto *ld1_264 = buffer.data(ld1 + 264);
    const auto *ld1_267 = buffer.data(ld1 + 267);
    const auto *ld1_269 = buffer.data(ld1 + 269);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_201 = buffer.data(lf + 201);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_207 = buffer.data(lf + 207);
    const auto *lf_208 = buffer.data(lf + 208);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_211 = buffer.data(lf + 211);
    const auto *lf_212 = buffer.data(lf + 212);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_217 = buffer.data(lf + 217);
    const auto *lf_218 = buffer.data(lf + 218);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_226 = buffer.data(lf + 226);
    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_271 = buffer.data(lf + 271);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_277 = buffer.data(lf + 277);
    const auto *lf_278 = buffer.data(lf + 278);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_281 = buffer.data(lf + 281);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_288 = buffer.data(lf + 288);
    const auto *lf_289 = buffer.data(lf + 289);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_302 = buffer.data(lf + 302);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_307 = buffer.data(lf + 307);
    const auto *lf_308 = buffer.data(lf + 308);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_310 = buffer.data(lf + 310);
    const auto *lf_312 = buffer.data(lf + 312);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_317 = buffer.data(lf + 317);
    const auto *lf_318 = buffer.data(lf + 318);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_320 = buffer.data(lf + 320);
    const auto *lf_322 = buffer.data(lf + 322);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_327 = buffer.data(lf + 327);
    const auto *lf_328 = buffer.data(lf + 328);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_330 = buffer.data(lf + 330);
    const auto *lf_332 = buffer.data(lf + 332);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_337 = buffer.data(lf + 337);
    const auto *lf_338 = buffer.data(lf + 338);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_340 = buffer.data(lf + 340);
    const auto *lf_342 = buffer.data(lf + 342);
    const auto *lf_346 = buffer.data(lf + 346);
    const auto *lf_347 = buffer.data(lf + 347);
    const auto *lf_348 = buffer.data(lf + 348);
    const auto *lf_350 = buffer.data(lf + 350);
    const auto *lf_352 = buffer.data(lf + 352);
    const auto *lf_355 = buffer.data(lf + 355);
    const auto *lf_356 = buffer.data(lf + 356);
    const auto *lf_357 = buffer.data(lf + 357);
    const auto *lf_359 = buffer.data(lf + 359);
    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_361 = buffer.data(lf + 361);
    const auto *lf_363 = buffer.data(lf + 363);
    const auto *lf_365 = buffer.data(lf + 365);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_367 = buffer.data(lf + 367);
    const auto *lf_368 = buffer.data(lf + 368);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_370 = buffer.data(lf + 370);
    const auto *lf_372 = buffer.data(lf + 372);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_377 = buffer.data(lf + 377);
    const auto *lf_378 = buffer.data(lf + 378);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_382 = buffer.data(lf + 382);
    const auto *lf_383 = buffer.data(lf + 383);
    const auto *lf_385 = buffer.data(lf + 385);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_387 = buffer.data(lf + 387);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_392 = buffer.data(lf + 392);
    const auto *lf_393 = buffer.data(lf + 393);
    const auto *lf_395 = buffer.data(lf + 395);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_397 = buffer.data(lf + 397);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_402 = buffer.data(lf + 402);
    const auto *lf_403 = buffer.data(lf + 403);
    const auto *lf_405 = buffer.data(lf + 405);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_407 = buffer.data(lf + 407);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_412 = buffer.data(lf + 412);
    const auto *lf_413 = buffer.data(lf + 413);
    const auto *lf_415 = buffer.data(lf + 415);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_417 = buffer.data(lf + 417);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);
    const auto *lf_422 = buffer.data(lf + 422);
    const auto *lf_423 = buffer.data(lf + 423);
    const auto *lf_425 = buffer.data(lf + 425);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_427 = buffer.data(lf + 427);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_430 = buffer.data(lf + 430);
    const auto *lf_432 = buffer.data(lf + 432);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_437 = buffer.data(lf + 437);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_439 = buffer.data(lf + 439);
    const auto *lf_440 = buffer.data(lf + 440);
    const auto *lf_442 = buffer.data(lf + 442);
    const auto *lf_443 = buffer.data(lf + 443);
    const auto *lf_445 = buffer.data(lf + 445);
    const auto *lf_446 = buffer.data(lf + 446);
    const auto *lf_447 = buffer.data(lf + 447);
    const auto *lf_448 = buffer.data(lf + 448);
    const auto *lf_449 = buffer.data(lf + 449);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, \
                         lf_0, lf_1, lf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kf_0[k]
                 + f_1 * ld0_0[k]
                 - f_2 * ld1_0[k]
                 + pb_x[k] * lf_0[k];

        t_1[k] = pb_y[k] * lf_0[k];

        t_2[k] = pb_z[k] * lf_0[k];

        t_3[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_y[k] * lf_1[k];

        t_4[k] = pb_y[k] * lf_2[k];

        t_5[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, kf_6, kf_9, ld0_3, ld1_3, \
                         lf_3, lf_5, lf_6, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * kf_6[k]
                 + pb_x[k] * lf_6[k];

        t_7[k] = pb_z[k] * lf_3[k];

        t_8[k] = pb_y[k] * lf_5[k];

        t_9[k] = f_0 * kf_9[k]
                 + pb_x[k] * lf_9[k];

        t_10[k] = f_1 * ld0_3[k]
                  - f_2 * ld1_3[k]
                  + pb_y[k] * lf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, kg_0, ld0_5, ld1_5, \
                         lf_6, lf_8, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * lf_6[k];

        t_12[k] = f_3 * ld0_5[k]
                  - f_4 * ld1_5[k]
                  + pb_y[k] * lf_8[k];

        t_13[k] = pb_y[k] * lf_9[k];

        t_14[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_9[k];

        t_15[k] = pa_y[k] * kg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, kf_0, kf_1, kg_3, \
                         kg_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * kf_0[k]
                  + pb_y[k] * lf_10[k];

        t_17[k] = pb_z[k] * lf_10[k];

        t_18[k] = f_6 * kf_1[k]
                  + pa_y[k] * kg_3[k];

        t_19[k] = pb_z[k] * lf_11[k];

        t_20[k] = pa_y[k] * kg_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, kf_6, kf_16, kf_18, \
                         kg_9, kg_10, lf_13, lf_16, lf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * kf_16[k]
                  + pb_x[k] * lf_16[k];

        t_22[k] = pb_z[k] * lf_13[k];

        t_23[k] = f_7 * kf_18[k]
                  + pb_x[k] * lf_18[k];

        t_24[k] = pa_y[k] * kg_9[k];

        t_25[k] = f_8 * kf_6[k]
                  + pa_y[k] * kg_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, kf_8, kf_9, \
                         kg_0, kg_12, kg_14, lf_16, lf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * lf_16[k];

        t_27[k] = f_6 * kf_8[k]
                  + pa_y[k] * kg_12[k];

        t_28[k] = f_5 * kf_9[k]
                  + pb_y[k] * lf_19[k];

        t_29[k] = pa_y[k] * kg_14[k];

        t_30[k] = pa_z[k] * kg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, kf_0, kf_2, \
                         kg_3, kg_5, kg_6, lf_20, lf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * lf_20[k];

        t_32[k] = f_5 * kf_0[k]
                  + pb_z[k] * lf_20[k];

        t_33[k] = pa_z[k] * kg_3[k];

        t_34[k] = pb_y[k] * lf_22[k];

        t_35[k] = f_6 * kf_2[k]
                  + pa_z[k] * kg_5[k];

        t_36[k] = pa_z[k] * kg_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, kf_27, kf_29, kg_10, lf_25, \
                         lf_27, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * kf_27[k]
                  + pb_x[k] * lf_27[k];

        t_38[k] = pb_y[k] * lf_25[k];

        t_39[k] = f_7 * kf_29[k]
                  + pb_x[k] * lf_29[k];

        t_40[k] = pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, kf_6, kf_7, kf_9, kg_12, \
                         kg_14, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * kf_6[k]
                  + pb_z[k] * lf_26[k];

        t_42[k] = f_6 * kf_7[k]
                  + pa_z[k] * kg_12[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_8 * kf_9[k]
                  + pa_z[k] * kg_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, ig0_0, ig1_0, kf_10, kg_15, \
                         lf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * ig0_0[k]
                  - f_10 * ig1_0[k]
                  + pa_y[k] * kg_15[k];

        t_46[k] = f_6 * kf_10[k]
                  + pb_y[k] * lf_30[k];

        t_47[k] = pb_z[k] * lf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, kf_33, kf_36, ld0_18, ld0_21, \
                         ld1_18, ld1_21, lf_31, lf_32, lf_33, lf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * kf_33[k]
                  + f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_x[k] * lf_33[k];

        t_49[k] = pb_z[k] * lf_31[k];

        t_50[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_z[k] * lf_32[k];

        t_51[k] = f_11 * kf_36[k]
                  + pb_x[k] * lf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, ig0_55, ig1_55, kf_38, \
                         kf_39, kg_55, lf_33, lf_38, lf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * lf_33[k];

        t_53[k] = f_11 * kf_38[k]
                  + pb_x[k] * lf_38[k];

        t_54[k] = f_11 * kf_39[k]
                  + pb_x[k] * lf_39[k];

        t_55[k] = f_12 * ig0_55[k]
                  - f_13 * ig1_55[k]
                  + pa_x[k] * kg_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, kf_19, ld0_21, ld0_23, ld1_21, \
                         ld1_23, lf_36, lf_37, lf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * lf_36[k];

        t_57[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_37[k];

        t_58[k] = f_6 * kf_19[k]
                  + pb_y[k] * lf_39[k];

        t_59[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, kf_22, kg_16, \
                         kg_18, kg_30, kg_32, kg_35, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * kg_30[k];

        t_61[k] = pa_z[k] * kg_16[k];

        t_62[k] = pa_y[k] * kg_32[k];

        t_63[k] = pa_z[k] * kg_18[k];

        t_64[k] = f_5 * kf_22[k]
                  + pb_y[k] * lf_42[k];

        t_65[k] = pa_y[k] * kg_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, kf_47, kf_48, kg_21, \
                         kg_25, kg_39, lf_47, lf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * kg_21[k];

        t_67[k] = f_11 * kf_47[k]
                  + pb_x[k] * lf_47[k];

        t_68[k] = f_11 * kf_48[k]
                  + pb_x[k] * lf_48[k];

        t_69[k] = pa_y[k] * kg_39[k];

        t_70[k] = pa_z[k] * kg_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, kf_16, kf_28, kf_29, kg_42, \
                         kg_44, lf_46, lf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * kf_16[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_6 * kf_28[k]
                  + pa_y[k] * kg_42[k];

        t_73[k] = f_5 * kf_29[k]
                  + pb_y[k] * lf_49[k];

        t_74[k] = pa_y[k] * kg_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, ig0_0, ig1_0, kf_20, kg_30, \
                         ld0_30, ld1_30, lf_50, lf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * ig0_0[k]
                  - f_10 * ig1_0[k]
                  + pa_z[k] * kg_30[k];

        t_76[k] = pb_y[k] * lf_50[k];

        t_77[k] = f_6 * kf_20[k]
                  + pb_z[k] * lf_50[k];

        t_78[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_51[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, kf_55, kf_56, kf_57, \
                         ld0_35, ld1_35, lf_52, lf_55, lf_56, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * lf_52[k];

        t_80[k] = f_11 * kf_55[k]
                  + f_3 * ld0_35[k]
                  - f_4 * ld1_35[k]
                  + pb_x[k] * lf_55[k];

        t_81[k] = f_11 * kf_56[k]
                  + pb_x[k] * lf_56[k];

        t_82[k] = f_11 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_83[k] = pb_y[k] * lf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, kf_26, kf_59, ld0_33, \
                         ld0_35, ld1_33, ld1_35, lf_56, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * kf_59[k]
                  + pb_x[k] * lf_59[k];

        t_85[k] = f_1 * ld0_33[k]
                  - f_2 * ld1_33[k]
                  + pb_y[k] * lf_56[k];

        t_86[k] = f_6 * kf_26[k]
                  + pb_z[k] * lf_56[k];

        t_87[k] = f_3 * ld0_35[k]
                  - f_4 * ld1_35[k]
                  + pb_y[k] * lf_58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, ig0_15, ig0_89, ig1_15, \
                         ig1_89, kf_30, kg_45, kg_89, lf_59, lf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * lf_59[k];

        t_89[k] = f_12 * ig0_89[k]
                  - f_13 * ig1_89[k]
                  + pa_x[k] * kg_89[k];

        t_90[k] = f_14 * ig0_15[k]
                  - f_15 * ig1_15[k]
                  + pa_y[k] * kg_45[k];

        t_91[k] = f_16 * kf_30[k]
                  + pb_y[k] * lf_60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, kf_63, ld0_36, ld0_39, ld1_36, \
                         ld1_39, lf_60, lf_61, lf_62, lf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * lf_60[k];

        t_93[k] = f_17 * kf_63[k]
                  + f_3 * ld0_39[k]
                  - f_4 * ld1_39[k]
                  + pb_x[k] * lf_63[k];

        t_94[k] = pb_z[k] * lf_61[k];

        t_95[k] = f_3 * ld0_36[k]
                  - f_4 * ld1_36[k]
                  + pb_z[k] * lf_62[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, kf_66, kf_68, kf_69, lf_63, \
                         lf_66, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_17 * kf_66[k]
                  + pb_x[k] * lf_66[k];

        t_97[k] = pb_z[k] * lf_63[k];

        t_98[k] = f_17 * kf_68[k]
                  + pb_x[k] * lf_68[k];

        t_99[k] = f_17 * kf_69[k]
                  + pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, ig0_100, ig1_100, \
                         kf_39, kg_100, ld0_39, ld1_39, lf_66, lf_67, \
                         lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_18 * ig0_100[k]
                   - f_19 * ig1_100[k]
                   + pa_x[k] * kg_100[k];

        t_101[k] = pb_z[k] * lf_66[k];

        t_102[k] = f_3 * ld0_39[k]
                   - f_4 * ld1_39[k]
                   + pb_z[k] * lf_67[k];

        t_103[k] = f_16 * kf_39[k]
                   + pb_y[k] * lf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, kf_30, kg_45, kg_46, \
                         kg_48, ld0_41, ld1_41, lf_69, lf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * ld0_41[k]
                   - f_2 * ld1_41[k]
                   + pb_z[k] * lf_69[k];

        t_105[k] = pa_z[k] * kg_45[k];

        t_106[k] = pa_z[k] * kg_46[k];

        t_107[k] = f_5 * kf_30[k]
                   + pb_z[k] * lf_70[k];

        t_108[k] = pa_z[k] * kg_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, kf_32, kf_42, kf_77, \
                         kg_50, kg_51, lf_72, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * kf_42[k]
                   + pb_y[k] * lf_72[k];

        t_110[k] = f_6 * kf_32[k]
                   + pa_z[k] * kg_50[k];

        t_111[k] = pa_z[k] * kg_51[k];

        t_112[k] = f_17 * kf_77[k]
                   + pb_x[k] * lf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, kf_36, kf_78, kf_79, \
                         kg_55, lf_76, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_17 * kf_78[k]
                   + pb_x[k] * lf_78[k];

        t_114[k] = f_17 * kf_79[k]
                   + pb_x[k] * lf_79[k];

        t_115[k] = pa_z[k] * kg_55[k];

        t_116[k] = f_5 * kf_36[k]
                   + pb_z[k] * lf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, kf_37, kf_39, kf_49, \
                         kg_57, kg_59, kg_75, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * kf_37[k]
                   + pa_z[k] * kg_57[k];

        t_118[k] = f_6 * kf_49[k]
                   + pb_y[k] * lf_79[k];

        t_119[k] = f_8 * kf_39[k]
                   + pa_z[k] * kg_59[k];

        t_120[k] = pa_y[k] * kg_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, kf_50, kf_51, kf_52, \
                         kg_77, kg_78, kg_80, lf_80, lf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * kf_50[k]
                   + pb_y[k] * lf_80[k];

        t_122[k] = pa_y[k] * kg_77[k];

        t_123[k] = f_6 * kf_51[k]
                   + pa_y[k] * kg_78[k];

        t_124[k] = f_5 * kf_52[k]
                   + pb_y[k] * lf_82[k];

        t_125[k] = pa_y[k] * kg_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, kf_56, kf_86, kf_87, \
                         kf_88, kg_84, kg_85, lf_86, lf_87, lf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_17 * kf_86[k]
                   + pb_x[k] * lf_86[k];

        t_127[k] = f_17 * kf_87[k]
                   + pb_x[k] * lf_87[k];

        t_128[k] = f_17 * kf_88[k]
                   + pb_x[k] * lf_88[k];

        t_129[k] = pa_y[k] * kg_84[k];

        t_130[k] = f_8 * kf_56[k]
                   + pa_y[k] * kg_85[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, kf_46, kf_58, kf_59, \
                         kg_87, kg_89, lf_86, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * kf_46[k]
                   + pb_z[k] * lf_86[k];

        t_132[k] = f_6 * kf_58[k]
                   + pa_y[k] * kg_87[k];

        t_133[k] = f_5 * kf_59[k]
                   + pb_y[k] * lf_89[k];

        t_134[k] = pa_y[k] * kg_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, ig0_30, ig1_30, kf_50, \
                         kg_75, ld0_54, ld1_54, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_14 * ig0_30[k]
                   - f_15 * ig1_30[k]
                   + pa_z[k] * kg_75[k];

        t_136[k] = pb_y[k] * lf_90[k];

        t_137[k] = f_16 * kf_50[k]
                   + pb_z[k] * lf_90[k];

        t_138[k] = f_3 * ld0_54[k]
                   - f_4 * ld1_54[k]
                   + pb_y[k] * lf_91[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, kf_95, kf_96, kf_97, \
                         ld0_59, ld1_59, lf_92, lf_95, lf_96, lf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * lf_92[k];

        t_140[k] = f_17 * kf_95[k]
                   + f_3 * ld0_59[k]
                   - f_4 * ld1_59[k]
                   + pb_x[k] * lf_95[k];

        t_141[k] = f_17 * kf_96[k]
                   + pb_x[k] * lf_96[k];

        t_142[k] = f_17 * kf_97[k]
                   + pb_x[k] * lf_97[k];

        t_143[k] = pb_y[k] * lf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, kf_56, kf_99, ld0_57, \
                         ld0_59, ld1_57, ld1_59, lf_96, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_17 * kf_99[k]
                   + pb_x[k] * lf_99[k];

        t_145[k] = f_1 * ld0_57[k]
                   - f_2 * ld1_57[k]
                   + pb_y[k] * lf_96[k];

        t_146[k] = f_16 * kf_56[k]
                   + pb_z[k] * lf_96[k];

        t_147[k] = f_3 * ld0_59[k]
                   - f_4 * ld1_59[k]
                   + pb_y[k] * lf_98[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pa_y, pb_y, ig0_45, ig0_149, \
                         ig1_45, ig1_149, kf_60, kg_90, kg_149, lf_99, \
                         lf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * lf_99[k];

        t_149[k] = f_18 * ig0_149[k]
                   - f_19 * ig1_149[k]
                   + pa_x[k] * kg_149[k];

        t_150[k] = f_20 * ig0_45[k]
                   - f_21 * ig1_45[k]
                   + pa_y[k] * kg_90[k];

        t_151[k] = f_8 * kf_60[k]
                   + pb_y[k] * lf_100[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_z, kf_103, ld0_60, ld0_63, \
                         ld1_60, ld1_63, lf_100, lf_101, lf_102, \
                         lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_z[k] * lf_100[k];

        t_153[k] = f_8 * kf_103[k]
                   + f_3 * ld0_63[k]
                   - f_4 * ld1_63[k]
                   + pb_x[k] * lf_103[k];

        t_154[k] = pb_z[k] * lf_101[k];

        t_155[k] = f_3 * ld0_60[k]
                   - f_4 * ld1_60[k]
                   + pb_z[k] * lf_102[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pb_z, kf_106, kf_108, kf_109, \
                         lf_103, lf_106, lf_108, lf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_8 * kf_106[k]
                   + pb_x[k] * lf_106[k];

        t_157[k] = pb_z[k] * lf_103[k];

        t_158[k] = f_8 * kf_108[k]
                   + pb_x[k] * lf_108[k];

        t_159[k] = f_8 * kf_109[k]
                   + pb_x[k] * lf_109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_y, pb_z, ig0_160, ig1_160, \
                         kf_69, kg_160, ld0_63, ld1_63, lf_106, lf_107, \
                         lf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_20 * ig0_160[k]
                   - f_21 * ig1_160[k]
                   + pa_x[k] * kg_160[k];

        t_161[k] = pb_z[k] * lf_106[k];

        t_162[k] = f_3 * ld0_63[k]
                   - f_4 * ld1_63[k]
                   + pb_z[k] * lf_107[k];

        t_163[k] = f_8 * kf_69[k]
                   + pb_y[k] * lf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_z, pb_z, kf_60, kg_90, kg_91, \
                         kg_93, ld0_65, ld1_65, lf_109, lf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * ld0_65[k]
                   - f_2 * ld1_65[k]
                   + pb_z[k] * lf_109[k];

        t_165[k] = pa_z[k] * kg_90[k];

        t_166[k] = pa_z[k] * kg_91[k];

        t_167[k] = f_5 * kf_60[k]
                   + pb_z[k] * lf_110[k];

        t_168[k] = pa_z[k] * kg_93[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_z, pb_x, pb_y, kf_62, kf_72, kf_117, \
                         kg_95, kg_96, lf_112, lf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_16 * kf_72[k]
                   + pb_y[k] * lf_112[k];

        t_170[k] = f_6 * kf_62[k]
                   + pa_z[k] * kg_95[k];

        t_171[k] = pa_z[k] * kg_96[k];

        t_172[k] = f_8 * kf_117[k]
                   + pb_x[k] * lf_117[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pb_x, pb_z, kf_66, kf_118, kf_119, \
                         kg_100, lf_116, lf_118, lf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_8 * kf_118[k]
                   + pb_x[k] * lf_118[k];

        t_174[k] = f_8 * kf_119[k]
                   + pb_x[k] * lf_119[k];

        t_175[k] = pa_z[k] * kg_100[k];

        t_176[k] = f_5 * kf_66[k]
                   + pb_z[k] * lf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pa_z, pb_y, ig0_75, ig1_75, kf_67, \
                         kf_69, kf_79, kg_102, kg_104, kg_120, lf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * kf_67[k]
                   + pa_z[k] * kg_102[k];

        t_178[k] = f_16 * kf_79[k]
                   + pb_y[k] * lf_119[k];

        t_179[k] = f_8 * kf_69[k]
                   + pa_z[k] * kg_104[k];

        t_180[k] = f_9 * ig0_75[k]
                   - f_10 * ig1_75[k]
                   + pa_y[k] * kg_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_z, pb_y, pb_z, ig0_48, ig1_48, kf_70, \
                         kf_80, kf_82, kg_108, lf_120, lf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * kf_80[k]
                   + pb_y[k] * lf_120[k];

        t_182[k] = f_6 * kf_70[k]
                   + pb_z[k] * lf_120[k];

        t_183[k] = f_9 * ig0_48[k]
                   - f_10 * ig1_48[k]
                   + pa_z[k] * kg_108[k];

        t_184[k] = f_6 * kf_82[k]
                   + pb_y[k] * lf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pb_x, ig0_80, ig1_80, kf_126, \
                         kf_127, kf_128, kg_125, lf_126, lf_127, \
                         lf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_9 * ig0_80[k]
                   - f_10 * ig1_80[k]
                   + pa_y[k] * kg_125[k];

        t_186[k] = f_8 * kf_126[k]
                   + pb_x[k] * lf_126[k];

        t_187[k] = f_8 * kf_127[k]
                   + pb_x[k] * lf_127[k];

        t_188[k] = f_8 * kf_128[k]
                   + pb_x[k] * lf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_x, pb_z, ig0_190, ig1_190, kf_76, \
                         kf_129, kg_190, lf_126, lf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_8 * kf_129[k]
                   + pb_x[k] * lf_129[k];

        t_190[k] = f_20 * ig0_190[k]
                   - f_21 * ig1_190[k]
                   + pa_x[k] * kg_190[k];

        t_191[k] = f_6 * kf_76[k]
                   + pb_z[k] * lf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_x, pa_y, pb_y, ig0_192, ig0_194, \
                         ig1_192, ig1_194, kf_89, kg_135, kg_192, kg_194, \
                         lf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_20 * ig0_192[k]
                   - f_21 * ig1_192[k]
                   + pa_x[k] * kg_192[k];

        t_193[k] = f_6 * kf_89[k]
                   + pb_y[k] * lf_129[k];

        t_194[k] = f_20 * ig0_194[k]
                   - f_21 * ig1_194[k]
                   + pa_x[k] * kg_194[k];

        t_195[k] = pa_y[k] * kg_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, kf_90, kf_91, kf_92, \
                         kg_137, kg_138, kg_140, lf_130, lf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * kf_90[k]
                   + pb_y[k] * lf_130[k];

        t_197[k] = pa_y[k] * kg_137[k];

        t_198[k] = f_6 * kf_91[k]
                   + pa_y[k] * kg_138[k];

        t_199[k] = f_5 * kf_92[k]
                   + pb_y[k] * lf_132[k];

        t_200[k] = pa_y[k] * kg_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pa_y, pb_x, kf_96, kf_136, kf_137, \
                         kf_138, kg_144, kg_145, lf_136, lf_137, \
                         lf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_8 * kf_136[k]
                   + pb_x[k] * lf_136[k];

        t_202[k] = f_8 * kf_137[k]
                   + pb_x[k] * lf_137[k];

        t_203[k] = f_8 * kf_138[k]
                   + pb_x[k] * lf_138[k];

        t_204[k] = pa_y[k] * kg_144[k];

        t_205[k] = f_8 * kf_96[k]
                   + pa_y[k] * kg_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_y, pb_y, pb_z, kf_86, kf_98, kf_99, \
                         kg_147, kg_149, lf_136, lf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_16 * kf_86[k]
                   + pb_z[k] * lf_136[k];

        t_207[k] = f_6 * kf_98[k]
                   + pa_y[k] * kg_147[k];

        t_208[k] = f_5 * kf_99[k]
                   + pb_y[k] * lf_139[k];

        t_209[k] = pa_y[k] * kg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_z, pb_y, pb_z, ig0_75, ig1_75, kf_90, \
                         kg_135, ld0_84, ld1_84, lf_140, lf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_20 * ig0_75[k]
                   - f_21 * ig1_75[k]
                   + pa_z[k] * kg_135[k];

        t_211[k] = pb_y[k] * lf_140[k];

        t_212[k] = f_8 * kf_90[k]
                   + pb_z[k] * lf_140[k];

        t_213[k] = f_3 * ld0_84[k]
                   - f_4 * ld1_84[k]
                   + pb_y[k] * lf_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, kf_145, kf_146, \
                         kf_147, ld0_89, ld1_89, lf_142, lf_145, lf_146, \
                         lf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * lf_142[k];

        t_215[k] = f_8 * kf_145[k]
                   + f_3 * ld0_89[k]
                   - f_4 * ld1_89[k]
                   + pb_x[k] * lf_145[k];

        t_216[k] = f_8 * kf_146[k]
                   + pb_x[k] * lf_146[k];

        t_217[k] = f_8 * kf_147[k]
                   + pb_x[k] * lf_147[k];

        t_218[k] = pb_y[k] * lf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pb_y, pb_z, kf_96, kf_149, ld0_87, \
                         ld0_89, ld1_87, ld1_89, lf_146, lf_148, \
                         lf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_8 * kf_149[k]
                   + pb_x[k] * lf_149[k];

        t_220[k] = f_1 * ld0_87[k]
                   - f_2 * ld1_87[k]
                   + pb_y[k] * lf_146[k];

        t_221[k] = f_8 * kf_96[k]
                   + pb_z[k] * lf_146[k];

        t_222[k] = f_3 * ld0_89[k]
                   - f_4 * ld1_89[k]
                   + pb_y[k] * lf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pa_y, pb_y, ig0_90, ig0_224, \
                         ig1_90, ig1_224, kf_100, kg_150, kg_224, lf_149, \
                         lf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_y[k] * lf_149[k];

        t_224[k] = f_20 * ig0_224[k]
                   - f_21 * ig1_224[k]
                   + pa_x[k] * kg_224[k];

        t_225[k] = f_18 * ig0_90[k]
                   - f_19 * ig1_90[k]
                   + pa_y[k] * kg_150[k];

        t_226[k] = f_17 * kf_100[k]
                   + pb_y[k] * lf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pb_x, pb_z, kf_153, ld0_90, ld0_93, \
                         ld1_90, ld1_93, lf_150, lf_151, lf_152, \
                         lf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_z[k] * lf_150[k];

        t_228[k] = f_16 * kf_153[k]
                   + f_3 * ld0_93[k]
                   - f_4 * ld1_93[k]
                   + pb_x[k] * lf_153[k];

        t_229[k] = pb_z[k] * lf_151[k];

        t_230[k] = f_3 * ld0_90[k]
                   - f_4 * ld1_90[k]
                   + pb_z[k] * lf_152[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pb_z, kf_156, kf_158, kf_159, \
                         lf_153, lf_156, lf_158, lf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_16 * kf_156[k]
                   + pb_x[k] * lf_156[k];

        t_232[k] = pb_z[k] * lf_153[k];

        t_233[k] = f_16 * kf_158[k]
                   + pb_x[k] * lf_158[k];

        t_234[k] = f_16 * kf_159[k]
                   + pb_x[k] * lf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pa_x, pb_y, pb_z, ig0_235, ig1_235, \
                         kf_109, kg_235, ld0_93, ld1_93, lf_156, lf_157, \
                         lf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_14 * ig0_235[k]
                   - f_15 * ig1_235[k]
                   + pa_x[k] * kg_235[k];

        t_236[k] = pb_z[k] * lf_156[k];

        t_237[k] = f_3 * ld0_93[k]
                   - f_4 * ld1_93[k]
                   + pb_z[k] * lf_157[k];

        t_238[k] = f_17 * kf_109[k]
                   + pb_y[k] * lf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_z, pb_z, kf_100, kg_150, \
                         kg_151, kg_153, ld0_95, ld1_95, lf_159, \
                         lf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * ld0_95[k]
                   - f_2 * ld1_95[k]
                   + pb_z[k] * lf_159[k];

        t_240[k] = pa_z[k] * kg_150[k];

        t_241[k] = pa_z[k] * kg_151[k];

        t_242[k] = f_5 * kf_100[k]
                   + pb_z[k] * lf_160[k];

        t_243[k] = pa_z[k] * kg_153[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pb_x, pb_y, kf_102, kf_112, kf_167, \
                         kg_155, kg_156, lf_162, lf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_8 * kf_112[k]
                   + pb_y[k] * lf_162[k];

        t_245[k] = f_6 * kf_102[k]
                   + pa_z[k] * kg_155[k];

        t_246[k] = pa_z[k] * kg_156[k];

        t_247[k] = f_16 * kf_167[k]
                   + pb_x[k] * lf_167[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_z, pb_x, pb_z, kf_106, kf_168, kf_169, \
                         kg_160, lf_166, lf_168, lf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_16 * kf_168[k]
                   + pb_x[k] * lf_168[k];

        t_249[k] = f_16 * kf_169[k]
                   + pb_x[k] * lf_169[k];

        t_250[k] = pa_z[k] * kg_160[k];

        t_251[k] = f_5 * kf_106[k]
                   + pb_z[k] * lf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_y, pa_z, pb_y, ig0_120, ig1_120, \
                         kf_107, kf_109, kf_119, kg_162, kg_164, kg_180, \
                         lf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_6 * kf_107[k]
                   + pa_z[k] * kg_162[k];

        t_253[k] = f_8 * kf_119[k]
                   + pb_y[k] * lf_169[k];

        t_254[k] = f_8 * kf_109[k]
                   + pa_z[k] * kg_164[k];

        t_255[k] = f_14 * ig0_120[k]
                   - f_15 * ig1_120[k]
                   + pa_y[k] * kg_180[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_z, pb_y, pb_z, ig0_93, ig1_93, kf_110, \
                         kf_120, kf_122, kg_168, lf_170, lf_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_16 * kf_120[k]
                   + pb_y[k] * lf_170[k];

        t_257[k] = f_6 * kf_110[k]
                   + pb_z[k] * lf_170[k];

        t_258[k] = f_9 * ig0_93[k]
                   - f_10 * ig1_93[k]
                   + pa_z[k] * kg_168[k];

        t_259[k] = f_16 * kf_122[k]
                   + pb_y[k] * lf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pb_x, ig0_125, ig1_125, kf_176, \
                         kf_177, kf_178, kg_185, lf_176, lf_177, \
                         lf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * ig0_125[k]
                   - f_15 * ig1_125[k]
                   + pa_y[k] * kg_185[k];

        t_261[k] = f_16 * kf_176[k]
                   + pb_x[k] * lf_176[k];

        t_262[k] = f_16 * kf_177[k]
                   + pb_x[k] * lf_177[k];

        t_263[k] = f_16 * kf_178[k]
                   + pb_x[k] * lf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_x, pb_x, pb_z, ig0_265, ig1_265, kf_116, \
                         kf_179, kg_265, lf_176, lf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_16 * kf_179[k]
                   + pb_x[k] * lf_179[k];

        t_265[k] = f_14 * ig0_265[k]
                   - f_15 * ig1_265[k]
                   + pa_x[k] * kg_265[k];

        t_266[k] = f_6 * kf_116[k]
                   + pb_z[k] * lf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_x, pb_y, ig0_267, ig0_269, ig1_267, ig1_269, \
                         kf_129, kg_267, kg_269, lf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * ig0_267[k]
                   - f_15 * ig1_267[k]
                   + pa_x[k] * kg_267[k];

        t_268[k] = f_16 * kf_129[k]
                   + pb_y[k] * lf_179[k];

        t_269[k] = f_14 * ig0_269[k]
                   - f_15 * ig1_269[k]
                   + pa_x[k] * kg_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_y, pb_y, pb_z, ig0_135, ig1_135, kf_120, \
                         kf_130, kg_195, lf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_9 * ig0_135[k]
                   - f_10 * ig1_135[k]
                   + pa_y[k] * kg_195[k];

        t_271[k] = f_6 * kf_130[k]
                   + pb_y[k] * lf_180[k];

        t_272[k] = f_16 * kf_120[k]
                   + pb_z[k] * lf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_y, pa_z, pb_y, ig0_108, ig0_140, ig1_108, \
                         ig1_140, kf_132, kg_183, kg_200, lf_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_14 * ig0_108[k]
                   - f_15 * ig1_108[k]
                   + pa_z[k] * kg_183[k];

        t_274[k] = f_6 * kf_132[k]
                   + pb_y[k] * lf_182[k];

        t_275[k] = f_9 * ig0_140[k]
                   - f_10 * ig1_140[k]
                   + pa_y[k] * kg_200[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, kf_186, kf_187, kf_188, kf_189, \
                         lf_186, lf_187, lf_188, lf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_16 * kf_186[k]
                   + pb_x[k] * lf_186[k];

        t_277[k] = f_16 * kf_187[k]
                   + pb_x[k] * lf_187[k];

        t_278[k] = f_16 * kf_188[k]
                   + pb_x[k] * lf_188[k];

        t_279[k] = f_16 * kf_189[k]
                   + pb_x[k] * lf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_x, pb_z, ig0_280, ig0_282, ig1_280, ig1_282, \
                         kf_126, kg_280, kg_282, lf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_14 * ig0_280[k]
                   - f_15 * ig1_280[k]
                   + pa_x[k] * kg_280[k];

        t_281[k] = f_16 * kf_126[k]
                   + pb_z[k] * lf_186[k];

        t_282[k] = f_14 * ig0_282[k]
                   - f_15 * ig1_282[k]
                   + pa_x[k] * kg_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pb_y, ig0_284, ig1_284, \
                         kf_139, kf_140, kg_210, kg_284, lf_189, \
                         lf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_6 * kf_139[k]
                   + pb_y[k] * lf_189[k];

        t_284[k] = f_14 * ig0_284[k]
                   - f_15 * ig1_284[k]
                   + pa_x[k] * kg_284[k];

        t_285[k] = pa_y[k] * kg_210[k];

        t_286[k] = f_5 * kf_140[k]
                   + pb_y[k] * lf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_y, pb_x, pb_y, kf_141, kf_142, \
                         kf_196, kg_212, kg_213, kg_215, lf_192, \
                         lf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_y[k] * kg_212[k];

        t_288[k] = f_6 * kf_141[k]
                   + pa_y[k] * kg_213[k];

        t_289[k] = f_5 * kf_142[k]
                   + pb_y[k] * lf_192[k];

        t_290[k] = pa_y[k] * kg_215[k];

        t_291[k] = f_16 * kf_196[k]
                   + pb_x[k] * lf_196[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_x, kf_146, kf_197, kf_198, \
                         kg_219, kg_220, lf_197, lf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_16 * kf_197[k]
                   + pb_x[k] * lf_197[k];

        t_293[k] = f_16 * kf_198[k]
                   + pb_x[k] * lf_198[k];

        t_294[k] = pa_y[k] * kg_219[k];

        t_295[k] = f_8 * kf_146[k]
                   + pa_y[k] * kg_220[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_y, pb_z, kf_136, kf_148, kf_149, \
                         kg_222, kg_224, lf_196, lf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_8 * kf_136[k]
                   + pb_z[k] * lf_196[k];

        t_297[k] = f_6 * kf_148[k]
                   + pa_y[k] * kg_222[k];

        t_298[k] = f_5 * kf_149[k]
                   + pb_y[k] * lf_199[k];

        t_299[k] = pa_y[k] * kg_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_z, pb_y, pb_z, ig0_135, ig1_135, \
                         kf_140, kg_210, ld0_120, ld1_120, lf_200, \
                         lf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_18 * ig0_135[k]
                   - f_19 * ig1_135[k]
                   + pa_z[k] * kg_210[k];

        t_301[k] = pb_y[k] * lf_200[k];

        t_302[k] = f_17 * kf_140[k]
                   + pb_z[k] * lf_200[k];

        t_303[k] = f_3 * ld0_120[k]
                   - f_4 * ld1_120[k]
                   + pb_y[k] * lf_201[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, kf_205, kf_206, \
                         kf_207, ld0_125, ld1_125, lf_202, lf_205, lf_206, \
                         lf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pb_y[k] * lf_202[k];

        t_305[k] = f_16 * kf_205[k]
                   + f_3 * ld0_125[k]
                   - f_4 * ld1_125[k]
                   + pb_x[k] * lf_205[k];

        t_306[k] = f_16 * kf_206[k]
                   + pb_x[k] * lf_206[k];

        t_307[k] = f_16 * kf_207[k]
                   + pb_x[k] * lf_207[k];

        t_308[k] = pb_y[k] * lf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pb_y, pb_z, kf_146, kf_209, \
                         ld0_123, ld0_125, ld1_123, ld1_125, lf_206, lf_208, \
                         lf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_16 * kf_209[k]
                   + pb_x[k] * lf_209[k];

        t_310[k] = f_1 * ld0_123[k]
                   - f_2 * ld1_123[k]
                   + pb_y[k] * lf_206[k];

        t_311[k] = f_17 * kf_146[k]
                   + pb_z[k] * lf_206[k];

        t_312[k] = f_3 * ld0_125[k]
                   - f_4 * ld1_125[k]
                   + pb_y[k] * lf_208[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_x, pa_y, pb_y, ig0_150, ig0_314, \
                         ig1_150, ig1_314, kf_150, kg_225, kg_314, lf_209, \
                         lf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * lf_209[k];

        t_314[k] = f_14 * ig0_314[k]
                   - f_15 * ig1_314[k]
                   + pa_x[k] * kg_314[k];

        t_315[k] = f_12 * ig0_150[k]
                   - f_13 * ig1_150[k]
                   + pa_y[k] * kg_225[k];

        t_316[k] = f_11 * kf_150[k]
                   + pb_y[k] * lf_210[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_x, pb_z, kf_213, ld0_126, ld0_129, \
                         ld1_126, ld1_129, lf_210, lf_211, lf_212, \
                         lf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pb_z[k] * lf_210[k];

        t_318[k] = f_6 * kf_213[k]
                   + f_3 * ld0_129[k]
                   - f_4 * ld1_129[k]
                   + pb_x[k] * lf_213[k];

        t_319[k] = pb_z[k] * lf_211[k];

        t_320[k] = f_3 * ld0_126[k]
                   - f_4 * ld1_126[k]
                   + pb_z[k] * lf_212[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_z, kf_216, kf_218, kf_219, \
                         lf_213, lf_216, lf_218, lf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_6 * kf_216[k]
                   + pb_x[k] * lf_216[k];

        t_322[k] = pb_z[k] * lf_213[k];

        t_323[k] = f_6 * kf_218[k]
                   + pb_x[k] * lf_218[k];

        t_324[k] = f_6 * kf_219[k]
                   + pb_x[k] * lf_219[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_x, pb_y, pb_z, ig0_325, ig1_325, \
                         kf_159, kg_325, ld0_129, ld1_129, lf_216, lf_217, \
                         lf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_9 * ig0_325[k]
                   - f_10 * ig1_325[k]
                   + pa_x[k] * kg_325[k];

        t_326[k] = pb_z[k] * lf_216[k];

        t_327[k] = f_3 * ld0_129[k]
                   - f_4 * ld1_129[k]
                   + pb_z[k] * lf_217[k];

        t_328[k] = f_11 * kf_159[k]
                   + pb_y[k] * lf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, pa_z, pb_z, kf_150, kg_225, \
                         kg_226, kg_228, ld0_131, ld1_131, lf_219, \
                         lf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * ld0_131[k]
                   - f_2 * ld1_131[k]
                   + pb_z[k] * lf_219[k];

        t_330[k] = pa_z[k] * kg_225[k];

        t_331[k] = pa_z[k] * kg_226[k];

        t_332[k] = f_5 * kf_150[k]
                   + pb_z[k] * lf_220[k];

        t_333[k] = pa_z[k] * kg_228[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_z, pb_x, pb_y, kf_152, kf_162, kf_227, \
                         kg_230, kg_231, lf_222, lf_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_17 * kf_162[k]
                   + pb_y[k] * lf_222[k];

        t_335[k] = f_6 * kf_152[k]
                   + pa_z[k] * kg_230[k];

        t_336[k] = pa_z[k] * kg_231[k];

        t_337[k] = f_6 * kf_227[k]
                   + pb_x[k] * lf_227[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_z, pb_x, pb_z, kf_156, kf_228, kf_229, \
                         kg_235, lf_226, lf_228, lf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_6 * kf_228[k]
                   + pb_x[k] * lf_228[k];

        t_339[k] = f_6 * kf_229[k]
                   + pb_x[k] * lf_229[k];

        t_340[k] = pa_z[k] * kg_235[k];

        t_341[k] = f_5 * kf_156[k]
                   + pb_z[k] * lf_226[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pa_y, pa_z, pb_y, ig0_180, ig1_180, \
                         kf_157, kf_159, kf_169, kg_237, kg_239, kg_255, \
                         lf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_6 * kf_157[k]
                   + pa_z[k] * kg_237[k];

        t_343[k] = f_17 * kf_169[k]
                   + pb_y[k] * lf_229[k];

        t_344[k] = f_8 * kf_159[k]
                   + pa_z[k] * kg_239[k];

        t_345[k] = f_20 * ig0_180[k]
                   - f_21 * ig1_180[k]
                   + pa_y[k] * kg_255[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_z, pb_y, pb_z, ig0_153, ig1_153, \
                         kf_160, kf_170, kf_172, kg_243, lf_230, \
                         lf_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_8 * kf_170[k]
                   + pb_y[k] * lf_230[k];

        t_347[k] = f_6 * kf_160[k]
                   + pb_z[k] * lf_230[k];

        t_348[k] = f_9 * ig0_153[k]
                   - f_10 * ig1_153[k]
                   + pa_z[k] * kg_243[k];

        t_349[k] = f_8 * kf_172[k]
                   + pb_y[k] * lf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_y, pb_x, ig0_185, ig1_185, kf_236, \
                         kf_237, kf_238, kg_260, lf_236, lf_237, \
                         lf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_20 * ig0_185[k]
                   - f_21 * ig1_185[k]
                   + pa_y[k] * kg_260[k];

        t_351[k] = f_6 * kf_236[k]
                   + pb_x[k] * lf_236[k];

        t_352[k] = f_6 * kf_237[k]
                   + pb_x[k] * lf_237[k];

        t_353[k] = f_6 * kf_238[k]
                   + pb_x[k] * lf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_x, pb_x, pb_z, ig0_355, ig1_355, kf_166, \
                         kf_239, kg_355, lf_236, lf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_6 * kf_239[k]
                   + pb_x[k] * lf_239[k];

        t_355[k] = f_9 * ig0_355[k]
                   - f_10 * ig1_355[k]
                   + pa_x[k] * kg_355[k];

        t_356[k] = f_6 * kf_166[k]
                   + pb_z[k] * lf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_x, pb_y, ig0_357, ig0_359, ig1_357, ig1_359, \
                         kf_179, kg_357, kg_359, lf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_9 * ig0_357[k]
                   - f_10 * ig1_357[k]
                   + pa_x[k] * kg_357[k];

        t_358[k] = f_8 * kf_179[k]
                   + pb_y[k] * lf_239[k];

        t_359[k] = f_9 * ig0_359[k]
                   - f_10 * ig1_359[k]
                   + pa_x[k] * kg_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pb_y, pb_z, ig0_195, ig1_195, kf_170, \
                         kf_180, kg_270, lf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_14 * ig0_195[k]
                   - f_15 * ig1_195[k]
                   + pa_y[k] * kg_270[k];

        t_361[k] = f_16 * kf_180[k]
                   + pb_y[k] * lf_240[k];

        t_362[k] = f_16 * kf_170[k]
                   + pb_z[k] * lf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_y, pa_z, pb_y, ig0_168, ig0_200, ig1_168, \
                         ig1_200, kf_182, kg_258, kg_275, lf_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * ig0_168[k]
                   - f_15 * ig1_168[k]
                   + pa_z[k] * kg_258[k];

        t_364[k] = f_16 * kf_182[k]
                   + pb_y[k] * lf_242[k];

        t_365[k] = f_14 * ig0_200[k]
                   - f_15 * ig1_200[k]
                   + pa_y[k] * kg_275[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, kf_246, kf_247, kf_248, kf_249, \
                         lf_246, lf_247, lf_248, lf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_6 * kf_246[k]
                   + pb_x[k] * lf_246[k];

        t_367[k] = f_6 * kf_247[k]
                   + pb_x[k] * lf_247[k];

        t_368[k] = f_6 * kf_248[k]
                   + pb_x[k] * lf_248[k];

        t_369[k] = f_6 * kf_249[k]
                   + pb_x[k] * lf_249[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pa_x, pb_z, ig0_370, ig0_372, ig1_370, ig1_372, \
                         kf_176, kg_370, kg_372, lf_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_9 * ig0_370[k]
                   - f_10 * ig1_370[k]
                   + pa_x[k] * kg_370[k];

        t_371[k] = f_16 * kf_176[k]
                   + pb_z[k] * lf_246[k];

        t_372[k] = f_9 * ig0_372[k]
                   - f_10 * ig1_372[k]
                   + pa_x[k] * kg_372[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pa_x, pa_y, pb_y, ig0_210, ig0_374, ig1_210, \
                         ig1_374, kf_189, kg_285, kg_374, lf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * kf_189[k]
                   + pb_y[k] * lf_249[k];

        t_374[k] = f_9 * ig0_374[k]
                   - f_10 * ig1_374[k]
                   + pa_x[k] * kg_374[k];

        t_375[k] = f_9 * ig0_210[k]
                   - f_10 * ig1_210[k]
                   + pa_y[k] * kg_285[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pb_y, pb_z, ig0_183, ig1_183, \
                         kf_180, kf_190, kf_192, kg_273, lf_250, \
                         lf_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_6 * kf_190[k]
                   + pb_y[k] * lf_250[k];

        t_377[k] = f_8 * kf_180[k]
                   + pb_z[k] * lf_250[k];

        t_378[k] = f_20 * ig0_183[k]
                   - f_21 * ig1_183[k]
                   + pa_z[k] * kg_273[k];

        t_379[k] = f_6 * kf_192[k]
                   + pb_y[k] * lf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pa_y, pb_x, ig0_215, ig1_215, kf_256, \
                         kf_257, kf_258, kg_290, lf_256, lf_257, \
                         lf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_9 * ig0_215[k]
                   - f_10 * ig1_215[k]
                   + pa_y[k] * kg_290[k];

        t_381[k] = f_6 * kf_256[k]
                   + pb_x[k] * lf_256[k];

        t_382[k] = f_6 * kf_257[k]
                   + pb_x[k] * lf_257[k];

        t_383[k] = f_6 * kf_258[k]
                   + pb_x[k] * lf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pa_x, pb_x, pb_z, ig0_385, ig1_385, kf_186, \
                         kf_259, kg_385, lf_256, lf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_6 * kf_259[k]
                   + pb_x[k] * lf_259[k];

        t_385[k] = f_9 * ig0_385[k]
                   - f_10 * ig1_385[k]
                   + pa_x[k] * kg_385[k];

        t_386[k] = f_8 * kf_186[k]
                   + pb_z[k] * lf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_x, pa_y, pb_y, ig0_387, ig0_389, \
                         ig1_387, ig1_389, kf_199, kg_300, kg_387, kg_389, \
                         lf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_9 * ig0_387[k]
                   - f_10 * ig1_387[k]
                   + pa_x[k] * kg_387[k];

        t_388[k] = f_6 * kf_199[k]
                   + pb_y[k] * lf_259[k];

        t_389[k] = f_9 * ig0_389[k]
                   - f_10 * ig1_389[k]
                   + pa_x[k] * kg_389[k];

        t_390[k] = pa_y[k] * kg_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, pa_y, pb_y, kf_200, kf_201, \
                         kf_202, kg_302, kg_303, kg_305, lf_260, \
                         lf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_5 * kf_200[k]
                   + pb_y[k] * lf_260[k];

        t_392[k] = pa_y[k] * kg_302[k];

        t_393[k] = f_6 * kf_201[k]
                   + pa_y[k] * kg_303[k];

        t_394[k] = f_5 * kf_202[k]
                   + pb_y[k] * lf_262[k];

        t_395[k] = pa_y[k] * kg_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, pa_y, pb_x, kf_206, kf_266, \
                         kf_267, kf_268, kg_309, kg_310, lf_266, lf_267, \
                         lf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_6 * kf_266[k]
                   + pb_x[k] * lf_266[k];

        t_397[k] = f_6 * kf_267[k]
                   + pb_x[k] * lf_267[k];

        t_398[k] = f_6 * kf_268[k]
                   + pb_x[k] * lf_268[k];

        t_399[k] = pa_y[k] * kg_309[k];

        t_400[k] = f_8 * kf_206[k]
                   + pa_y[k] * kg_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pb_y, pb_z, kf_196, kf_208, kf_209, \
                         kg_312, kg_314, lf_266, lf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_17 * kf_196[k]
                   + pb_z[k] * lf_266[k];

        t_402[k] = f_6 * kf_208[k]
                   + pa_y[k] * kg_312[k];

        t_403[k] = f_5 * kf_209[k]
                   + pb_y[k] * lf_269[k];

        t_404[k] = pa_y[k] * kg_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_z, pb_y, pb_z, ig0_210, ig1_210, \
                         kf_200, kg_300, ld0_162, ld1_162, lf_270, \
                         lf_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_12 * ig0_210[k]
                   - f_13 * ig1_210[k]
                   + pa_z[k] * kg_300[k];

        t_406[k] = pb_y[k] * lf_270[k];

        t_407[k] = f_11 * kf_200[k]
                   + pb_z[k] * lf_270[k];

        t_408[k] = f_3 * ld0_162[k]
                   - f_4 * ld1_162[k]
                   + pb_y[k] * lf_271[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pb_x, pb_y, kf_275, kf_276, \
                         kf_277, ld0_167, ld1_167, lf_272, lf_275, lf_276, \
                         lf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * lf_272[k];

        t_410[k] = f_6 * kf_275[k]
                   + f_3 * ld0_167[k]
                   - f_4 * ld1_167[k]
                   + pb_x[k] * lf_275[k];

        t_411[k] = f_6 * kf_276[k]
                   + pb_x[k] * lf_276[k];

        t_412[k] = f_6 * kf_277[k]
                   + pb_x[k] * lf_277[k];

        t_413[k] = pb_y[k] * lf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pb_y, pb_z, kf_206, kf_279, \
                         ld0_165, ld0_167, ld1_165, ld1_167, lf_276, lf_278, \
                         lf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_6 * kf_279[k]
                   + pb_x[k] * lf_279[k];

        t_415[k] = f_1 * ld0_165[k]
                   - f_2 * ld1_165[k]
                   + pb_y[k] * lf_276[k];

        t_416[k] = f_11 * kf_206[k]
                   + pb_z[k] * lf_276[k];

        t_417[k] = f_3 * ld0_167[k]
                   - f_4 * ld1_167[k]
                   + pb_y[k] * lf_278[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pa_x, pb_y, pb_z, ig0_419, \
                         ig1_419, kf_210, kf_280, kg_419, kg_420, lf_279, \
                         lf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = pb_y[k] * lf_279[k];

        t_419[k] = f_9 * ig0_419[k]
                   - f_10 * ig1_419[k]
                   + pa_x[k] * kg_419[k];

        t_420[k] = f_8 * kf_280[k]
                   + pa_x[k] * kg_420[k];

        t_421[k] = f_7 * kf_210[k]
                   + pb_y[k] * lf_280[k];

        t_422[k] = pb_z[k] * lf_280[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, pa_x, pb_x, pb_z, kf_283, kf_285, \
                         kf_286, kg_423, kg_425, lf_281, lf_283, \
                         lf_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_6 * kf_283[k]
                   + pa_x[k] * kg_423[k];

        t_424[k] = pb_z[k] * lf_281[k];

        t_425[k] = f_6 * kf_285[k]
                   + pa_x[k] * kg_425[k];

        t_426[k] = f_5 * kf_286[k]
                   + pb_x[k] * lf_286[k];

        t_427[k] = pb_z[k] * lf_283[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pa_x, pb_x, pb_z, kf_288, kf_289, \
                         kg_430, kg_432, lf_286, lf_288, lf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_5 * kf_288[k]
                   + pb_x[k] * lf_288[k];

        t_429[k] = f_5 * kf_289[k]
                   + pb_x[k] * lf_289[k];

        t_430[k] = pa_x[k] * kg_430[k];

        t_431[k] = pb_z[k] * lf_286[k];

        t_432[k] = pa_x[k] * kg_432[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, t_438, pa_x, pa_z, pb_z, kf_210, \
                         kg_315, kg_316, kg_318, kg_433, kg_434, \
                         lf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = pa_x[k] * kg_433[k];

        t_434[k] = pa_x[k] * kg_434[k];

        t_435[k] = pa_z[k] * kg_315[k];

        t_436[k] = pa_z[k] * kg_316[k];

        t_437[k] = f_5 * kf_210[k]
                   + pb_z[k] * lf_290[k];

        t_438[k] = pa_z[k] * kg_318[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pa_x, pa_z, pb_x, pb_y, kf_222, kf_295, \
                         kf_297, kg_321, kg_440, lf_292, lf_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_11 * kf_222[k]
                   + pb_y[k] * lf_292[k];

        t_440[k] = f_6 * kf_295[k]
                   + pa_x[k] * kg_440[k];

        t_441[k] = pa_z[k] * kg_321[k];

        t_442[k] = f_5 * kf_297[k]
                   + pb_x[k] * lf_297[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, pa_x, pb_x, kf_298, kf_299, \
                         kg_445, kg_446, kg_447, kg_448, lf_298, \
                         lf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_5 * kf_298[k]
                   + pb_x[k] * lf_298[k];

        t_444[k] = f_5 * kf_299[k]
                   + pb_x[k] * lf_299[k];

        t_445[k] = pa_x[k] * kg_445[k];

        t_446[k] = pa_x[k] * kg_446[k];

        t_447[k] = pa_x[k] * kg_447[k];

        t_448[k] = pa_x[k] * kg_448[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, pa_x, pb_y, pb_z, kf_220, kf_230, \
                         kf_300, kf_303, kg_449, kg_450, kg_453, \
                         lf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_x[k] * kg_449[k];

        t_450[k] = f_8 * kf_300[k]
                   + pa_x[k] * kg_450[k];

        t_451[k] = f_17 * kf_230[k]
                   + pb_y[k] * lf_300[k];

        t_452[k] = f_6 * kf_220[k]
                   + pb_z[k] * lf_300[k];

        t_453[k] = f_6 * kf_303[k]
                   + pa_x[k] * kg_453[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, pa_x, pb_x, pb_y, kf_232, kf_305, kf_306, \
                         kf_307, kg_455, lf_302, lf_306, lf_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_17 * kf_232[k]
                   + pb_y[k] * lf_302[k];

        t_455[k] = f_6 * kf_305[k]
                   + pa_x[k] * kg_455[k];

        t_456[k] = f_5 * kf_306[k]
                   + pb_x[k] * lf_306[k];

        t_457[k] = f_5 * kf_307[k]
                   + pb_x[k] * lf_307[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, t_463, pa_x, pb_x, kf_308, kf_309, \
                         kg_460, kg_461, kg_462, kg_463, lf_308, \
                         lf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_5 * kf_308[k]
                   + pb_x[k] * lf_308[k];

        t_459[k] = f_5 * kf_309[k]
                   + pb_x[k] * lf_309[k];

        t_460[k] = pa_x[k] * kg_460[k];

        t_461[k] = pa_x[k] * kg_461[k];

        t_462[k] = pa_x[k] * kg_462[k];

        t_463[k] = pa_x[k] * kg_463[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, pa_x, pb_y, pb_z, kf_230, kf_240, \
                         kf_310, kf_313, kg_464, kg_465, kg_468, \
                         lf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pa_x[k] * kg_464[k];

        t_465[k] = f_8 * kf_310[k]
                   + pa_x[k] * kg_465[k];

        t_466[k] = f_8 * kf_240[k]
                   + pb_y[k] * lf_310[k];

        t_467[k] = f_16 * kf_230[k]
                   + pb_z[k] * lf_310[k];

        t_468[k] = f_6 * kf_313[k]
                   + pa_x[k] * kg_468[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pa_x, pb_x, pb_y, kf_242, kf_315, kf_316, \
                         kf_317, kg_470, lf_312, lf_316, lf_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_8 * kf_242[k]
                   + pb_y[k] * lf_312[k];

        t_470[k] = f_6 * kf_315[k]
                   + pa_x[k] * kg_470[k];

        t_471[k] = f_5 * kf_316[k]
                   + pb_x[k] * lf_316[k];

        t_472[k] = f_5 * kf_317[k]
                   + pb_x[k] * lf_317[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, t_478, pa_x, pb_x, kf_318, kf_319, \
                         kg_475, kg_476, kg_477, kg_478, lf_318, \
                         lf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_5 * kf_318[k]
                   + pb_x[k] * lf_318[k];

        t_474[k] = f_5 * kf_319[k]
                   + pb_x[k] * lf_319[k];

        t_475[k] = pa_x[k] * kg_475[k];

        t_476[k] = pa_x[k] * kg_476[k];

        t_477[k] = pa_x[k] * kg_477[k];

        t_478[k] = pa_x[k] * kg_478[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, pa_x, pb_y, pb_z, kf_240, kf_250, \
                         kf_320, kf_323, kg_479, kg_480, kg_483, \
                         lf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pa_x[k] * kg_479[k];

        t_480[k] = f_8 * kf_320[k]
                   + pa_x[k] * kg_480[k];

        t_481[k] = f_16 * kf_250[k]
                   + pb_y[k] * lf_320[k];

        t_482[k] = f_8 * kf_240[k]
                   + pb_z[k] * lf_320[k];

        t_483[k] = f_6 * kf_323[k]
                   + pa_x[k] * kg_483[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pa_x, pb_x, pb_y, kf_252, kf_325, kf_326, \
                         kf_327, kg_485, lf_322, lf_326, lf_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_16 * kf_252[k]
                   + pb_y[k] * lf_322[k];

        t_485[k] = f_6 * kf_325[k]
                   + pa_x[k] * kg_485[k];

        t_486[k] = f_5 * kf_326[k]
                   + pb_x[k] * lf_326[k];

        t_487[k] = f_5 * kf_327[k]
                   + pb_x[k] * lf_327[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, t_493, pa_x, pb_x, kf_328, kf_329, \
                         kg_490, kg_491, kg_492, kg_493, lf_328, \
                         lf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * kf_328[k]
                   + pb_x[k] * lf_328[k];

        t_489[k] = f_5 * kf_329[k]
                   + pb_x[k] * lf_329[k];

        t_490[k] = pa_x[k] * kg_490[k];

        t_491[k] = pa_x[k] * kg_491[k];

        t_492[k] = pa_x[k] * kg_492[k];

        t_493[k] = pa_x[k] * kg_493[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, pa_x, pb_y, pb_z, kf_250, kf_260, \
                         kf_330, kf_333, kg_494, kg_495, kg_498, \
                         lf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = pa_x[k] * kg_494[k];

        t_495[k] = f_8 * kf_330[k]
                   + pa_x[k] * kg_495[k];

        t_496[k] = f_6 * kf_260[k]
                   + pb_y[k] * lf_330[k];

        t_497[k] = f_17 * kf_250[k]
                   + pb_z[k] * lf_330[k];

        t_498[k] = f_6 * kf_333[k]
                   + pa_x[k] * kg_498[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, pa_x, pb_x, pb_y, kf_262, kf_335, kf_336, \
                         kf_337, kg_500, lf_332, lf_336, lf_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_6 * kf_262[k]
                   + pb_y[k] * lf_332[k];

        t_500[k] = f_6 * kf_335[k]
                   + pa_x[k] * kg_500[k];

        t_501[k] = f_5 * kf_336[k]
                   + pb_x[k] * lf_336[k];

        t_502[k] = f_5 * kf_337[k]
                   + pb_x[k] * lf_337[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, t_508, pa_x, pb_x, kf_338, kf_339, \
                         kg_505, kg_506, kg_507, kg_508, lf_338, \
                         lf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_5 * kf_338[k]
                   + pb_x[k] * lf_338[k];

        t_504[k] = f_5 * kf_339[k]
                   + pb_x[k] * lf_339[k];

        t_505[k] = pa_x[k] * kg_505[k];

        t_506[k] = pa_x[k] * kg_506[k];

        t_507[k] = pa_x[k] * kg_507[k];

        t_508[k] = pa_x[k] * kg_508[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, pa_x, pa_y, pb_y, kf_270, kf_343, \
                         kg_405, kg_407, kg_509, kg_513, lf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pa_x[k] * kg_509[k];

        t_510[k] = pa_y[k] * kg_405[k];

        t_511[k] = f_5 * kf_270[k]
                   + pb_y[k] * lf_340[k];

        t_512[k] = pa_y[k] * kg_407[k];

        t_513[k] = f_6 * kf_343[k]
                   + pa_x[k] * kg_513[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, pa_y, pb_x, pb_y, kf_272, kf_346, kf_347, \
                         kg_410, lf_342, lf_346, lf_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_5 * kf_272[k]
                   + pb_y[k] * lf_342[k];

        t_515[k] = pa_y[k] * kg_410[k];

        t_516[k] = f_5 * kf_346[k]
                   + pb_x[k] * lf_346[k];

        t_517[k] = f_5 * kf_347[k]
                   + pb_x[k] * lf_347[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, t_522, t_523, pa_x, pa_y, pb_x, kf_348, \
                         kg_414, kg_520, kg_521, kg_522, kg_523, \
                         lf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_5 * kf_348[k]
                   + pb_x[k] * lf_348[k];

        t_519[k] = pa_y[k] * kg_414[k];

        t_520[k] = pa_x[k] * kg_520[k];

        t_521[k] = pa_x[k] * kg_521[k];

        t_522[k] = pa_x[k] * kg_522[k];

        t_523[k] = pa_x[k] * kg_523[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, pa_x, pb_y, pb_z, kf_270, kf_350, \
                         kf_353, kg_524, kg_525, kg_528, lf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = pa_x[k] * kg_524[k];

        t_525[k] = f_8 * kf_350[k]
                   + pa_x[k] * kg_525[k];

        t_526[k] = pb_y[k] * lf_350[k];

        t_527[k] = f_7 * kf_270[k]
                   + pb_z[k] * lf_350[k];

        t_528[k] = f_6 * kf_353[k]
                   + pa_x[k] * kg_528[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, pa_x, pb_x, pb_y, kf_355, kf_356, \
                         kf_357, kg_530, lf_352, lf_355, lf_356, \
                         lf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = pb_y[k] * lf_352[k];

        t_530[k] = f_6 * kf_355[k]
                   + pa_x[k] * kg_530[k];

        t_531[k] = f_5 * kf_356[k]
                   + pb_x[k] * lf_356[k];

        t_532[k] = f_5 * kf_357[k]
                   + pb_x[k] * lf_357[k];

        t_533[k] = pb_y[k] * lf_355[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, t_539, pa_x, pb_x, pb_y, kf_359, \
                         kg_535, kg_536, kg_537, kg_539, lf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_5 * kf_359[k]
                   + pb_x[k] * lf_359[k];

        t_535[k] = pa_x[k] * kg_535[k];

        t_536[k] = pa_x[k] * kg_536[k];

        t_537[k] = pa_x[k] * kg_537[k];

        t_538[k] = pb_y[k] * lf_359[k];

        t_539[k] = pa_x[k] * kg_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, pb_x, pb_y, pb_z, kf_280, ld0_216, \
                         ld0_219, ld1_216, ld1_219, lf_360, lf_361, \
                         lf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * ld0_216[k]
                   - f_2 * ld1_216[k]
                   + pb_x[k] * lf_360[k];

        t_541[k] = f_0 * kf_280[k]
                   + pb_y[k] * lf_360[k];

        t_542[k] = pb_z[k] * lf_360[k];

        t_543[k] = f_3 * ld0_219[k]
                   - f_4 * ld1_219[k]
                   + pb_x[k] * lf_363[k];

        t_544[k] = pb_z[k] * lf_361[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, pb_x, ld0_221, ld1_221, lf_365, \
                         lf_366, lf_367, lf_368, lf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_3 * ld0_221[k]
                   - f_4 * ld1_221[k]
                   + pb_x[k] * lf_365[k];

        t_546[k] = pb_x[k] * lf_366[k];

        t_547[k] = pb_x[k] * lf_367[k];

        t_548[k] = pb_x[k] * lf_368[k];

        t_549[k] = pb_x[k] * lf_369[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, pb_y, pb_z, kf_286, kf_289, \
                         ld0_219, ld0_221, ld1_219, ld1_221, lf_366, lf_367, \
                         lf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_0 * kf_286[k]
                   + f_1 * ld0_219[k]
                   - f_2 * ld1_219[k]
                   + pb_y[k] * lf_366[k];

        t_551[k] = pb_z[k] * lf_366[k];

        t_552[k] = f_3 * ld0_219[k]
                   - f_4 * ld1_219[k]
                   + pb_z[k] * lf_367[k];

        t_553[k] = f_0 * kf_289[k]
                   + pb_y[k] * lf_369[k];

        t_554[k] = f_1 * ld0_221[k]
                   - f_2 * ld1_221[k]
                   + pb_z[k] * lf_369[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, pa_z, pb_y, pb_z, kf_280, kf_292, \
                         kg_420, kg_421, kg_423, lf_370, lf_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_z[k] * kg_420[k];

        t_556[k] = pa_z[k] * kg_421[k];

        t_557[k] = f_5 * kf_280[k]
                   + pb_z[k] * lf_370[k];

        t_558[k] = pa_z[k] * kg_423[k];

        t_559[k] = f_7 * kf_292[k]
                   + pb_y[k] * lf_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, pa_z, pb_x, kf_282, kg_425, \
                         kg_430, lf_376, lf_377, lf_378, lf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_6 * kf_282[k]
                   + pa_z[k] * kg_425[k];

        t_561[k] = pb_x[k] * lf_376[k];

        t_562[k] = pb_x[k] * lf_377[k];

        t_563[k] = pb_x[k] * lf_378[k];

        t_564[k] = pb_x[k] * lf_379[k];

        t_565[k] = pa_z[k] * kg_430[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_z, pb_y, pb_z, kf_286, kf_287, kf_289, \
                         kf_299, kg_432, kg_434, lf_376, lf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_5 * kf_286[k]
                   + pb_z[k] * lf_376[k];

        t_567[k] = f_6 * kf_287[k]
                   + pa_z[k] * kg_432[k];

        t_568[k] = f_7 * kf_299[k]
                   + pb_y[k] * lf_379[k];

        t_569[k] = f_8 * kf_289[k]
                   + pa_z[k] * kg_434[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pb_x, pb_y, pb_z, kf_290, kf_300, \
                         ld0_228, ld0_231, ld1_228, ld1_231, lf_380, \
                         lf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_1 * ld0_228[k]
                   - f_2 * ld1_228[k]
                   + pb_x[k] * lf_380[k];

        t_571[k] = f_11 * kf_300[k]
                   + pb_y[k] * lf_380[k];

        t_572[k] = f_6 * kf_290[k]
                   + pb_z[k] * lf_380[k];

        t_573[k] = f_3 * ld0_231[k]
                   - f_4 * ld1_231[k]
                   + pb_x[k] * lf_383[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, pb_x, pb_y, kf_302, ld0_233, \
                         ld1_233, lf_382, lf_385, lf_386, lf_387, \
                         lf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_11 * kf_302[k]
                   + pb_y[k] * lf_382[k];

        t_575[k] = f_3 * ld0_233[k]
                   - f_4 * ld1_233[k]
                   + pb_x[k] * lf_385[k];

        t_576[k] = pb_x[k] * lf_386[k];

        t_577[k] = pb_x[k] * lf_387[k];

        t_578[k] = pb_x[k] * lf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pa_z, pb_x, pb_z, ig0_325, ig1_325, kf_296, \
                         kg_445, lf_386, lf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = pb_x[k] * lf_389[k];

        t_580[k] = f_9 * ig0_325[k]
                   - f_10 * ig1_325[k]
                   + pa_z[k] * kg_445[k];

        t_581[k] = f_6 * kf_296[k]
                   + pb_z[k] * lf_386[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pa_y, pb_y, ig0_359, ig1_359, kf_308, kf_309, \
                         kg_464, ld0_233, ld1_233, lf_388, lf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_11 * kf_308[k]
                   + f_3 * ld0_233[k]
                   - f_4 * ld1_233[k]
                   + pb_y[k] * lf_388[k];

        t_583[k] = f_11 * kf_309[k]
                   + pb_y[k] * lf_389[k];

        t_584[k] = f_12 * ig0_359[k]
                   - f_13 * ig1_359[k]
                   + pa_y[k] * kg_464[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pb_x, pb_y, pb_z, kf_300, kf_310, \
                         ld0_234, ld0_237, ld1_234, ld1_237, lf_390, \
                         lf_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_1 * ld0_234[k]
                   - f_2 * ld1_234[k]
                   + pb_x[k] * lf_390[k];

        t_586[k] = f_17 * kf_310[k]
                   + pb_y[k] * lf_390[k];

        t_587[k] = f_16 * kf_300[k]
                   + pb_z[k] * lf_390[k];

        t_588[k] = f_3 * ld0_237[k]
                   - f_4 * ld1_237[k]
                   + pb_x[k] * lf_393[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, pb_x, pb_y, kf_312, ld0_239, \
                         ld1_239, lf_392, lf_395, lf_396, lf_397, \
                         lf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_17 * kf_312[k]
                   + pb_y[k] * lf_392[k];

        t_590[k] = f_3 * ld0_239[k]
                   - f_4 * ld1_239[k]
                   + pb_x[k] * lf_395[k];

        t_591[k] = pb_x[k] * lf_396[k];

        t_592[k] = pb_x[k] * lf_397[k];

        t_593[k] = pb_x[k] * lf_398[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pa_z, pb_x, pb_z, ig0_340, ig1_340, kf_306, \
                         kg_460, lf_396, lf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = pb_x[k] * lf_399[k];

        t_595[k] = f_14 * ig0_340[k]
                   - f_15 * ig1_340[k]
                   + pa_z[k] * kg_460[k];

        t_596[k] = f_16 * kf_306[k]
                   + pb_z[k] * lf_396[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pa_y, pb_y, ig0_374, ig1_374, kf_318, kf_319, \
                         kg_479, ld0_239, ld1_239, lf_398, lf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_17 * kf_318[k]
                   + f_3 * ld0_239[k]
                   - f_4 * ld1_239[k]
                   + pb_y[k] * lf_398[k];

        t_598[k] = f_17 * kf_319[k]
                   + pb_y[k] * lf_399[k];

        t_599[k] = f_18 * ig0_374[k]
                   - f_19 * ig1_374[k]
                   + pa_y[k] * kg_479[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pb_x, pb_y, pb_z, kf_310, kf_320, \
                         ld0_240, ld0_243, ld1_240, ld1_243, lf_400, \
                         lf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_1 * ld0_240[k]
                   - f_2 * ld1_240[k]
                   + pb_x[k] * lf_400[k];

        t_601[k] = f_8 * kf_320[k]
                   + pb_y[k] * lf_400[k];

        t_602[k] = f_8 * kf_310[k]
                   + pb_z[k] * lf_400[k];

        t_603[k] = f_3 * ld0_243[k]
                   - f_4 * ld1_243[k]
                   + pb_x[k] * lf_403[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pb_x, pb_y, kf_322, ld0_245, \
                         ld1_245, lf_402, lf_405, lf_406, lf_407, \
                         lf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_8 * kf_322[k]
                   + pb_y[k] * lf_402[k];

        t_605[k] = f_3 * ld0_245[k]
                   - f_4 * ld1_245[k]
                   + pb_x[k] * lf_405[k];

        t_606[k] = pb_x[k] * lf_406[k];

        t_607[k] = pb_x[k] * lf_407[k];

        t_608[k] = pb_x[k] * lf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pa_z, pb_x, pb_z, ig0_355, ig1_355, kf_316, \
                         kg_475, lf_406, lf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pb_x[k] * lf_409[k];

        t_610[k] = f_20 * ig0_355[k]
                   - f_21 * ig1_355[k]
                   + pa_z[k] * kg_475[k];

        t_611[k] = f_8 * kf_316[k]
                   + pb_z[k] * lf_406[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_y, pb_y, ig0_389, ig1_389, kf_328, kf_329, \
                         kg_494, ld0_245, ld1_245, lf_408, lf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_8 * kf_328[k]
                   + f_3 * ld0_245[k]
                   - f_4 * ld1_245[k]
                   + pb_y[k] * lf_408[k];

        t_613[k] = f_8 * kf_329[k]
                   + pb_y[k] * lf_409[k];

        t_614[k] = f_20 * ig0_389[k]
                   - f_21 * ig1_389[k]
                   + pa_y[k] * kg_494[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pb_x, pb_y, pb_z, kf_320, kf_330, \
                         ld0_246, ld0_249, ld1_246, ld1_249, lf_410, \
                         lf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_1 * ld0_246[k]
                   - f_2 * ld1_246[k]
                   + pb_x[k] * lf_410[k];

        t_616[k] = f_16 * kf_330[k]
                   + pb_y[k] * lf_410[k];

        t_617[k] = f_17 * kf_320[k]
                   + pb_z[k] * lf_410[k];

        t_618[k] = f_3 * ld0_249[k]
                   - f_4 * ld1_249[k]
                   + pb_x[k] * lf_413[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, pb_x, pb_y, kf_332, ld0_251, \
                         ld1_251, lf_412, lf_415, lf_416, lf_417, \
                         lf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_16 * kf_332[k]
                   + pb_y[k] * lf_412[k];

        t_620[k] = f_3 * ld0_251[k]
                   - f_4 * ld1_251[k]
                   + pb_x[k] * lf_415[k];

        t_621[k] = pb_x[k] * lf_416[k];

        t_622[k] = pb_x[k] * lf_417[k];

        t_623[k] = pb_x[k] * lf_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pb_x, pb_z, ig0_370, ig1_370, kf_326, \
                         kg_490, lf_416, lf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pb_x[k] * lf_419[k];

        t_625[k] = f_18 * ig0_370[k]
                   - f_19 * ig1_370[k]
                   + pa_z[k] * kg_490[k];

        t_626[k] = f_17 * kf_326[k]
                   + pb_z[k] * lf_416[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_y, pb_y, ig0_404, ig1_404, kf_338, kf_339, \
                         kg_509, ld0_251, ld1_251, lf_418, lf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_16 * kf_338[k]
                   + f_3 * ld0_251[k]
                   - f_4 * ld1_251[k]
                   + pb_y[k] * lf_418[k];

        t_628[k] = f_16 * kf_339[k]
                   + pb_y[k] * lf_419[k];

        t_629[k] = f_14 * ig0_404[k]
                   - f_15 * ig1_404[k]
                   + pa_y[k] * kg_509[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pb_x, pb_y, pb_z, kf_330, kf_340, \
                         ld0_252, ld0_255, ld1_252, ld1_255, lf_420, \
                         lf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_1 * ld0_252[k]
                   - f_2 * ld1_252[k]
                   + pb_x[k] * lf_420[k];

        t_631[k] = f_6 * kf_340[k]
                   + pb_y[k] * lf_420[k];

        t_632[k] = f_11 * kf_330[k]
                   + pb_z[k] * lf_420[k];

        t_633[k] = f_3 * ld0_255[k]
                   - f_4 * ld1_255[k]
                   + pb_x[k] * lf_423[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pb_x, pb_y, kf_342, ld0_257, \
                         ld1_257, lf_422, lf_425, lf_426, lf_427, \
                         lf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_6 * kf_342[k]
                   + pb_y[k] * lf_422[k];

        t_635[k] = f_3 * ld0_257[k]
                   - f_4 * ld1_257[k]
                   + pb_x[k] * lf_425[k];

        t_636[k] = pb_x[k] * lf_426[k];

        t_637[k] = pb_x[k] * lf_427[k];

        t_638[k] = pb_x[k] * lf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pa_z, pb_x, pb_z, ig0_385, ig1_385, kf_336, \
                         kg_505, lf_426, lf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pb_x[k] * lf_429[k];

        t_640[k] = f_12 * ig0_385[k]
                   - f_13 * ig1_385[k]
                   + pa_z[k] * kg_505[k];

        t_641[k] = f_11 * kf_336[k]
                   + pb_z[k] * lf_426[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_y, pb_y, ig0_419, ig1_419, kf_348, \
                         kf_349, kg_524, kg_525, ld0_257, ld1_257, lf_428, \
                         lf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_6 * kf_348[k]
                   + f_3 * ld0_257[k]
                   - f_4 * ld1_257[k]
                   + pb_y[k] * lf_428[k];

        t_643[k] = f_6 * kf_349[k]
                   + pb_y[k] * lf_429[k];

        t_644[k] = f_9 * ig0_419[k]
                   - f_10 * ig1_419[k]
                   + pa_y[k] * kg_524[k];

        t_645[k] = pa_y[k] * kg_525[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, t_650, pa_y, pb_y, kf_350, kf_351, \
                         kf_352, kg_527, kg_528, kg_530, lf_430, \
                         lf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_5 * kf_350[k]
                   + pb_y[k] * lf_430[k];

        t_647[k] = pa_y[k] * kg_527[k];

        t_648[k] = f_6 * kf_351[k]
                   + pa_y[k] * kg_528[k];

        t_649[k] = f_5 * kf_352[k]
                   + pb_y[k] * lf_432[k];

        t_650[k] = pa_y[k] * kg_530[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, t_655, t_656, pa_y, pb_x, pb_z, kf_346, \
                         kf_356, kg_535, lf_436, lf_437, lf_438, \
                         lf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = pb_x[k] * lf_436[k];

        t_652[k] = pb_x[k] * lf_437[k];

        t_653[k] = pb_x[k] * lf_438[k];

        t_654[k] = pb_x[k] * lf_439[k];

        t_655[k] = f_8 * kf_356[k]
                   + pa_y[k] * kg_535[k];

        t_656[k] = f_7 * kf_346[k]
                   + pb_z[k] * lf_436[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, pa_y, pb_x, pb_y, kf_358, kf_359, \
                         kg_537, kg_539, ld0_264, ld1_264, lf_439, \
                         lf_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_6 * kf_358[k]
                   + pa_y[k] * kg_537[k];

        t_658[k] = f_5 * kf_359[k]
                   + pb_y[k] * lf_439[k];

        t_659[k] = pa_y[k] * kg_539[k];

        t_660[k] = f_1 * ld0_264[k]
                   - f_2 * ld1_264[k]
                   + pb_x[k] * lf_440[k];

        t_661[k] = pb_y[k] * lf_440[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pb_x, pb_y, pb_z, kf_350, ld0_267, \
                         ld0_269, ld1_267, ld1_269, lf_440, lf_442, lf_443, \
                         lf_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_0 * kf_350[k]
                   + pb_z[k] * lf_440[k];

        t_663[k] = f_3 * ld0_267[k]
                   - f_4 * ld1_267[k]
                   + pb_x[k] * lf_443[k];

        t_664[k] = pb_y[k] * lf_442[k];

        t_665[k] = f_3 * ld0_269[k]
                   - f_4 * ld1_269[k]
                   + pb_x[k] * lf_445[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, t_670, t_671, pb_x, pb_y, pb_z, kf_356, \
                         ld0_267, ld1_267, lf_446, lf_447, lf_448, \
                         lf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = pb_x[k] * lf_446[k];

        t_667[k] = pb_x[k] * lf_447[k];

        t_668[k] = pb_x[k] * lf_448[k];

        t_669[k] = pb_x[k] * lf_449[k];

        t_670[k] = f_1 * ld0_267[k]
                   - f_2 * ld1_267[k]
                   + pb_y[k] * lf_446[k];

        t_671[k] = f_0 * kf_356[k]
                   + pb_z[k] * lf_446[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, pb_y, pb_z, kf_359, ld0_269, ld1_269, lf_448, \
                         lf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_3 * ld0_269[k]
                   - f_4 * ld1_269[k]
                   + pb_y[k] * lf_448[k];

        t_673[k] = pb_y[k] * lf_449[k];

        t_674[k] = f_0 * kf_359[k]
                   + f_1 * ld0_269[k]
                   - f_2 * ld1_269[k]
                   + pb_z[k] * lf_449[k];
    }
}

}  // namespace simdt2ceri
