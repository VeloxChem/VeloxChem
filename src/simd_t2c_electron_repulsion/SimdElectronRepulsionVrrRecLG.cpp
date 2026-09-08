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
    const auto *ig0_1 = buffer.data(ig0 + 1);
    const auto *ig0_2 = buffer.data(ig0 + 2);
    const auto *ig0_3 = buffer.data(ig0 + 3);
    const auto *ig0_4 = buffer.data(ig0 + 4);
    const auto *ig0_5 = buffer.data(ig0 + 5);
    const auto *ig0_6 = buffer.data(ig0 + 6);
    const auto *ig0_7 = buffer.data(ig0 + 7);
    const auto *ig0_8 = buffer.data(ig0 + 8);
    const auto *ig0_9 = buffer.data(ig0 + 9);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_12 = buffer.data(ig0 + 12);
    const auto *ig0_13 = buffer.data(ig0 + 13);
    const auto *ig0_14 = buffer.data(ig0 + 14);
    const auto *ig0_15 = buffer.data(ig0 + 15);
    const auto *ig0_16 = buffer.data(ig0 + 16);
    const auto *ig0_17 = buffer.data(ig0 + 17);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_19 = buffer.data(ig0 + 19);
    const auto *ig0_20 = buffer.data(ig0 + 20);
    const auto *ig0_21 = buffer.data(ig0 + 21);
    const auto *ig0_22 = buffer.data(ig0 + 22);
    const auto *ig0_23 = buffer.data(ig0 + 23);
    const auto *ig0_24 = buffer.data(ig0 + 24);
    const auto *ig0_25 = buffer.data(ig0 + 25);
    const auto *ig0_26 = buffer.data(ig0 + 26);
    const auto *ig0_27 = buffer.data(ig0 + 27);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_31 = buffer.data(ig0 + 31);
    const auto *ig0_32 = buffer.data(ig0 + 32);
    const auto *ig0_33 = buffer.data(ig0 + 33);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_35 = buffer.data(ig0 + 35);
    const auto *ig0_36 = buffer.data(ig0 + 36);
    const auto *ig0_37 = buffer.data(ig0 + 37);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_39 = buffer.data(ig0 + 39);
    const auto *ig0_40 = buffer.data(ig0 + 40);
    const auto *ig0_41 = buffer.data(ig0 + 41);
    const auto *ig0_42 = buffer.data(ig0 + 42);
    const auto *ig0_43 = buffer.data(ig0 + 43);
    const auto *ig0_44 = buffer.data(ig0 + 44);
    const auto *ig0_45 = buffer.data(ig0 + 45);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);
    const auto *ig0_48 = buffer.data(ig0 + 48);
    const auto *ig0_49 = buffer.data(ig0 + 49);
    const auto *ig0_50 = buffer.data(ig0 + 50);
    const auto *ig0_51 = buffer.data(ig0 + 51);
    const auto *ig0_52 = buffer.data(ig0 + 52);
    const auto *ig0_53 = buffer.data(ig0 + 53);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_1 = buffer.data(ig1 + 1);
    const auto *ig1_2 = buffer.data(ig1 + 2);
    const auto *ig1_3 = buffer.data(ig1 + 3);
    const auto *ig1_4 = buffer.data(ig1 + 4);
    const auto *ig1_5 = buffer.data(ig1 + 5);
    const auto *ig1_6 = buffer.data(ig1 + 6);
    const auto *ig1_7 = buffer.data(ig1 + 7);
    const auto *ig1_8 = buffer.data(ig1 + 8);
    const auto *ig1_9 = buffer.data(ig1 + 9);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_12 = buffer.data(ig1 + 12);
    const auto *ig1_13 = buffer.data(ig1 + 13);
    const auto *ig1_14 = buffer.data(ig1 + 14);
    const auto *ig1_15 = buffer.data(ig1 + 15);
    const auto *ig1_16 = buffer.data(ig1 + 16);
    const auto *ig1_17 = buffer.data(ig1 + 17);
    const auto *ig1_18 = buffer.data(ig1 + 18);
    const auto *ig1_19 = buffer.data(ig1 + 19);
    const auto *ig1_20 = buffer.data(ig1 + 20);
    const auto *ig1_21 = buffer.data(ig1 + 21);
    const auto *ig1_22 = buffer.data(ig1 + 22);
    const auto *ig1_23 = buffer.data(ig1 + 23);
    const auto *ig1_24 = buffer.data(ig1 + 24);
    const auto *ig1_25 = buffer.data(ig1 + 25);
    const auto *ig1_26 = buffer.data(ig1 + 26);
    const auto *ig1_27 = buffer.data(ig1 + 27);
    const auto *ig1_28 = buffer.data(ig1 + 28);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_30 = buffer.data(ig1 + 30);
    const auto *ig1_31 = buffer.data(ig1 + 31);
    const auto *ig1_32 = buffer.data(ig1 + 32);
    const auto *ig1_33 = buffer.data(ig1 + 33);
    const auto *ig1_34 = buffer.data(ig1 + 34);
    const auto *ig1_35 = buffer.data(ig1 + 35);
    const auto *ig1_36 = buffer.data(ig1 + 36);
    const auto *ig1_37 = buffer.data(ig1 + 37);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_39 = buffer.data(ig1 + 39);
    const auto *ig1_40 = buffer.data(ig1 + 40);
    const auto *ig1_41 = buffer.data(ig1 + 41);
    const auto *ig1_42 = buffer.data(ig1 + 42);
    const auto *ig1_43 = buffer.data(ig1 + 43);
    const auto *ig1_44 = buffer.data(ig1 + 44);
    const auto *ig1_45 = buffer.data(ig1 + 45);
    const auto *ig1_46 = buffer.data(ig1 + 46);
    const auto *ig1_47 = buffer.data(ig1 + 47);
    const auto *ig1_48 = buffer.data(ig1 + 48);
    const auto *ig1_49 = buffer.data(ig1 + 49);
    const auto *ig1_50 = buffer.data(ig1 + 50);
    const auto *ig1_51 = buffer.data(ig1 + 51);
    const auto *ig1_52 = buffer.data(ig1 + 52);
    const auto *ig1_53 = buffer.data(ig1 + 53);

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
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
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
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_1 = buffer.data(kg + 1);
    const auto *kg_2 = buffer.data(kg + 2);
    const auto *kg_3 = buffer.data(kg + 3);
    const auto *kg_4 = buffer.data(kg + 4);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_6 = buffer.data(kg + 6);
    const auto *kg_7 = buffer.data(kg + 7);
    const auto *kg_8 = buffer.data(kg + 8);
    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_18 = buffer.data(kg + 18);
    const auto *kg_19 = buffer.data(kg + 19);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_27 = buffer.data(kg + 27);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_36 = buffer.data(kg + 36);
    const auto *kg_37 = buffer.data(kg + 37);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_40 = buffer.data(kg + 40);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_43 = buffer.data(kg + 43);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_45 = buffer.data(kg + 45);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_54 = buffer.data(kg + 54);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_64 = buffer.data(kg + 64);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_79 = buffer.data(kg + 79);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
    const auto *kg_84 = buffer.data(kg + 84);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_92 = buffer.data(kg + 92);
    const auto *kg_93 = buffer.data(kg + 93);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_97 = buffer.data(kg + 97);
    const auto *kg_98 = buffer.data(kg + 98);
    const auto *kg_99 = buffer.data(kg + 99);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_146 = buffer.data(kg + 146);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_148 = buffer.data(kg + 148);
    const auto *kg_149 = buffer.data(kg + 149);
    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_152 = buffer.data(kg + 152);
    const auto *kg_153 = buffer.data(kg + 153);
    const auto *kg_154 = buffer.data(kg + 154);
    const auto *kg_155 = buffer.data(kg + 155);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_157 = buffer.data(kg + 157);
    const auto *kg_158 = buffer.data(kg + 158);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_160 = buffer.data(kg + 160);
    const auto *kg_161 = buffer.data(kg + 161);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_163 = buffer.data(kg + 163);
    const auto *kg_164 = buffer.data(kg + 164);
    const auto *kg_165 = buffer.data(kg + 165);
    const auto *kg_166 = buffer.data(kg + 166);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_169 = buffer.data(kg + 169);
    const auto *kg_170 = buffer.data(kg + 170);
    const auto *kg_171 = buffer.data(kg + 171);
    const auto *kg_172 = buffer.data(kg + 172);
    const auto *kg_173 = buffer.data(kg + 173);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_175 = buffer.data(kg + 175);
    const auto *kg_176 = buffer.data(kg + 176);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_178 = buffer.data(kg + 178);
    const auto *kg_179 = buffer.data(kg + 179);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_181 = buffer.data(kg + 181);
    const auto *kg_182 = buffer.data(kg + 182);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_184 = buffer.data(kg + 184);
    const auto *kg_185 = buffer.data(kg + 185);
    const auto *kg_186 = buffer.data(kg + 186);
    const auto *kg_187 = buffer.data(kg + 187);
    const auto *kg_188 = buffer.data(kg + 188);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_111 = buffer.data(lf + 111);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_113 = buffer.data(lf + 113);
    const auto *lf_114 = buffer.data(lf + 114);
    const auto *lf_115 = buffer.data(lf + 115);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_121 = buffer.data(lf + 121);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_123 = buffer.data(lf + 123);
    const auto *lf_124 = buffer.data(lf + 124);
    const auto *lf_125 = buffer.data(lf + 125);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_131 = buffer.data(lf + 131);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_133 = buffer.data(lf + 133);
    const auto *lf_134 = buffer.data(lf + 134);
    const auto *lf_135 = buffer.data(lf + 135);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_143 = buffer.data(lf + 143);
    const auto *lf_144 = buffer.data(lf + 144);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_154 = buffer.data(lf + 154);
    const auto *lf_155 = buffer.data(lf + 155);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_161 = buffer.data(lf + 161);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_163 = buffer.data(lf + 163);
    const auto *lf_164 = buffer.data(lf + 164);
    const auto *lf_165 = buffer.data(lf + 165);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_171 = buffer.data(lf + 171);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_173 = buffer.data(lf + 173);
    const auto *lf_174 = buffer.data(lf + 174);
    const auto *lf_175 = buffer.data(lf + 175);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_181 = buffer.data(lf + 181);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_183 = buffer.data(lf + 183);
    const auto *lf_184 = buffer.data(lf + 184);
    const auto *lf_185 = buffer.data(lf + 185);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_191 = buffer.data(lf + 191);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_193 = buffer.data(lf + 193);
    const auto *lf_194 = buffer.data(lf + 194);
    const auto *lf_195 = buffer.data(lf + 195);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_201 = buffer.data(lf + 201);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_203 = buffer.data(lf + 203);
    const auto *lf_204 = buffer.data(lf + 204);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_207 = buffer.data(lf + 207);
    const auto *lf_208 = buffer.data(lf + 208);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_211 = buffer.data(lf + 211);
    const auto *lf_212 = buffer.data(lf + 212);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_214 = buffer.data(lf + 214);
    const auto *lf_215 = buffer.data(lf + 215);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_217 = buffer.data(lf + 217);
    const auto *lf_218 = buffer.data(lf + 218);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_221 = buffer.data(lf + 221);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_223 = buffer.data(lf + 223);
    const auto *lf_224 = buffer.data(lf + 224);
    const auto *lf_225 = buffer.data(lf + 225);
    const auto *lf_226 = buffer.data(lf + 226);
    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_231 = buffer.data(lf + 231);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_233 = buffer.data(lf + 233);
    const auto *lf_234 = buffer.data(lf + 234);
    const auto *lf_235 = buffer.data(lf + 235);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_241 = buffer.data(lf + 241);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_243 = buffer.data(lf + 243);
    const auto *lf_244 = buffer.data(lf + 244);
    const auto *lf_245 = buffer.data(lf + 245);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_251 = buffer.data(lf + 251);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_253 = buffer.data(lf + 253);
    const auto *lf_254 = buffer.data(lf + 254);
    const auto *lf_255 = buffer.data(lf + 255);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_261 = buffer.data(lf + 261);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_263 = buffer.data(lf + 263);
    const auto *lf_264 = buffer.data(lf + 264);
    const auto *lf_265 = buffer.data(lf + 265);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_271 = buffer.data(lf + 271);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_273 = buffer.data(lf + 273);
    const auto *lf_274 = buffer.data(lf + 274);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_277 = buffer.data(lf + 277);
    const auto *lf_278 = buffer.data(lf + 278);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_281 = buffer.data(lf + 281);
    const auto *lf_282 = buffer.data(lf + 282);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_284 = buffer.data(lf + 284);
    const auto *lf_285 = buffer.data(lf + 285);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_287 = buffer.data(lf + 287);
    const auto *lf_288 = buffer.data(lf + 288);
    const auto *lf_289 = buffer.data(lf + 289);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_291 = buffer.data(lf + 291);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_293 = buffer.data(lf + 293);
    const auto *lf_294 = buffer.data(lf + 294);
    const auto *lf_295 = buffer.data(lf + 295);
    const auto *lf_296 = buffer.data(lf + 296);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_301 = buffer.data(lf + 301);
    const auto *lf_302 = buffer.data(lf + 302);

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

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, kf_3, kf_6, ld0_1, ld1_1, \
                         lf_3, lf_4, lf_5, lf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * kf_3[k]
                 + pb_x[k] * lf_5[k];

        t_7[k] = pb_z[k] * lf_3[k];

        t_8[k] = pb_y[k] * lf_4[k];

        t_9[k] = f_0 * kf_6[k]
                 + pb_x[k] * lf_7[k];

        t_10[k] = f_1 * ld0_1[k]
                  - f_2 * ld1_1[k]
                  + pb_y[k] * lf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, kg_0, ld0_2, ld1_2, \
                         lf_5, lf_6, lf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * lf_5[k];

        t_12[k] = f_3 * ld0_2[k]
                  - f_4 * ld1_2[k]
                  + pb_y[k] * lf_6[k];

        t_13[k] = pb_y[k] * lf_7[k];

        t_14[k] = f_1 * ld0_2[k]
                  - f_2 * ld1_2[k]
                  + pb_z[k] * lf_7[k];

        t_15[k] = pa_y[k] * kg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, kf_0, kf_1, kg_1, \
                         kg_2, lf_8, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * kf_0[k]
                  + pb_y[k] * lf_8[k];

        t_17[k] = pb_z[k] * lf_8[k];

        t_18[k] = f_6 * kf_1[k]
                  + pa_y[k] * kg_1[k];

        t_19[k] = pb_z[k] * lf_9[k];

        t_20[k] = pa_y[k] * kg_2[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, kf_3, kf_8, kf_9, \
                         kg_4, kg_5, lf_10, lf_11, lf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * kf_8[k]
                  + pb_x[k] * lf_11[k];

        t_22[k] = pb_z[k] * lf_10[k];

        t_23[k] = f_7 * kf_9[k]
                  + pb_x[k] * lf_12[k];

        t_24[k] = pa_y[k] * kg_4[k];

        t_25[k] = f_8 * kf_3[k]
                  + pa_y[k] * kg_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, kf_5, kf_6, \
                         kg_0, kg_6, kg_7, lf_11, lf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * lf_11[k];

        t_27[k] = f_6 * kf_5[k]
                  + pa_y[k] * kg_6[k];

        t_28[k] = f_5 * kf_6[k]
                  + pb_y[k] * lf_13[k];

        t_29[k] = pa_y[k] * kg_7[k];

        t_30[k] = pa_z[k] * kg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, kf_0, kf_2, \
                         kg_1, kg_2, kg_3, lf_14, lf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * lf_14[k];

        t_32[k] = f_5 * kf_0[k]
                  + pb_z[k] * lf_14[k];

        t_33[k] = pa_z[k] * kg_1[k];

        t_34[k] = pb_y[k] * lf_15[k];

        t_35[k] = f_6 * kf_2[k]
                  + pa_z[k] * kg_2[k];

        t_36[k] = pa_z[k] * kg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, kf_14, kf_16, kg_5, lf_16, \
                         lf_18, lf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * kf_14[k]
                  + pb_x[k] * lf_18[k];

        t_38[k] = pb_y[k] * lf_16[k];

        t_39[k] = f_7 * kf_16[k]
                  + pb_x[k] * lf_19[k];

        t_40[k] = pa_z[k] * kg_5[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, kf_3, kf_4, kf_6, kg_6, \
                         kg_7, lf_17, lf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * kf_3[k]
                  + pb_z[k] * lf_17[k];

        t_42[k] = f_6 * kf_4[k]
                  + pa_z[k] * kg_6[k];

        t_43[k] = pb_y[k] * lf_19[k];

        t_44[k] = f_8 * kf_6[k]
                  + pa_z[k] * kg_7[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, ig0_0, ig1_0, kf_7, kg_8, \
                         lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * ig0_0[k]
                  - f_10 * ig1_0[k]
                  + pa_y[k] * kg_8[k];

        t_46[k] = f_6 * kf_7[k]
                  + pb_y[k] * lf_20[k];

        t_47[k] = pb_z[k] * lf_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, kf_19, kf_20, ld0_3, ld0_4, \
                         ld1_3, ld1_4, lf_21, lf_22, lf_23, lf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * kf_19[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_23[k];

        t_49[k] = pb_z[k] * lf_21[k];

        t_50[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_22[k];

        t_51[k] = f_11 * kf_20[k]
                  + pb_x[k] * lf_24[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, ig0_5, ig1_5, kf_22, kf_23, \
                         kg_24, lf_23, lf_26, lf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * lf_23[k];

        t_53[k] = f_11 * kf_22[k]
                  + pb_x[k] * lf_26[k];

        t_54[k] = f_11 * kf_23[k]
                  + pb_x[k] * lf_27[k];

        t_55[k] = f_12 * ig0_5[k]
                  - f_13 * ig1_5[k]
                  + pa_x[k] * kg_24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, kf_10, ld0_4, ld0_5, ld1_4, \
                         ld1_5, lf_24, lf_25, lf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * lf_24[k];

        t_57[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_25[k];

        t_58[k] = f_6 * kf_10[k]
                  + pb_y[k] * lf_27[k];

        t_59[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, kf_12, kg_9, \
                         kg_10, kg_13, kg_14, kg_15, lf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * kg_13[k];

        t_61[k] = pa_z[k] * kg_9[k];

        t_62[k] = pa_y[k] * kg_14[k];

        t_63[k] = pa_z[k] * kg_10[k];

        t_64[k] = f_5 * kf_12[k]
                  + pb_y[k] * lf_28[k];

        t_65[k] = pa_y[k] * kg_15[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, kf_26, kf_27, kg_11, \
                         kg_12, kg_16, lf_30, lf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * kg_11[k];

        t_67[k] = f_11 * kf_26[k]
                  + pb_x[k] * lf_30[k];

        t_68[k] = f_11 * kf_27[k]
                  + pb_x[k] * lf_31[k];

        t_69[k] = pa_y[k] * kg_16[k];

        t_70[k] = pa_z[k] * kg_12[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, kf_8, kf_15, kf_16, kg_17, \
                         kg_18, lf_29, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * kf_8[k]
                  + pb_z[k] * lf_29[k];

        t_72[k] = f_6 * kf_15[k]
                  + pa_y[k] * kg_17[k];

        t_73[k] = f_5 * kf_16[k]
                  + pb_y[k] * lf_32[k];

        t_74[k] = pa_y[k] * kg_18[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, ig0_0, ig1_0, kf_11, kg_13, \
                         ld0_6, ld1_6, lf_33, lf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * ig0_0[k]
                  - f_10 * ig1_0[k]
                  + pa_z[k] * kg_13[k];

        t_76[k] = pb_y[k] * lf_33[k];

        t_77[k] = f_6 * kf_11[k]
                  + pb_z[k] * lf_33[k];

        t_78[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_34[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, kf_32, kf_33, kf_34, ld0_8, \
                         ld1_8, lf_35, lf_36, lf_37, lf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * lf_35[k];

        t_80[k] = f_11 * kf_32[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_36[k];

        t_81[k] = f_11 * kf_33[k]
                  + pb_x[k] * lf_37[k];

        t_82[k] = f_11 * kf_34[k]
                  + pb_x[k] * lf_38[k];

        t_83[k] = pb_y[k] * lf_36[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, kf_13, kf_36, ld0_7, ld0_8, \
                         ld1_7, ld1_8, lf_37, lf_39, lf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * kf_36[k]
                  + pb_x[k] * lf_40[k];

        t_85[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_37[k];

        t_86[k] = f_6 * kf_13[k]
                  + pb_z[k] * lf_37[k];

        t_87[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_39[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, ig0_1, ig0_8, ig1_1, ig1_8, \
                         kf_17, kg_19, kg_34, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * lf_40[k];

        t_89[k] = f_12 * ig0_8[k]
                  - f_13 * ig1_8[k]
                  + pa_x[k] * kg_34[k];

        t_90[k] = f_14 * ig0_1[k]
                  - f_15 * ig1_1[k]
                  + pa_y[k] * kg_19[k];

        t_91[k] = f_16 * kf_17[k]
                  + pb_y[k] * lf_41[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, kf_39, ld0_9, ld0_10, ld1_9, \
                         ld1_10, lf_41, lf_42, lf_43, lf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * lf_41[k];

        t_93[k] = f_17 * kf_39[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_44[k];

        t_94[k] = pb_z[k] * lf_42[k];

        t_95[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_43[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, kf_40, kf_42, kf_43, lf_44, \
                         lf_45, lf_47, lf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_17 * kf_40[k]
                  + pb_x[k] * lf_45[k];

        t_97[k] = pb_z[k] * lf_44[k];

        t_98[k] = f_17 * kf_42[k]
                  + pb_x[k] * lf_47[k];

        t_99[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_48[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, ig0_11, ig1_11, kf_23, \
                         kg_40, ld0_10, ld1_10, lf_45, lf_46, lf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_18 * ig0_11[k]
                   - f_19 * ig1_11[k]
                   + pa_x[k] * kg_40[k];

        t_101[k] = pb_z[k] * lf_45[k];

        t_102[k] = f_3 * ld0_10[k]
                   - f_4 * ld1_10[k]
                   + pb_z[k] * lf_46[k];

        t_103[k] = f_16 * kf_23[k]
                   + pb_y[k] * lf_48[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, kf_17, kg_19, kg_20, \
                         kg_21, ld0_11, ld1_11, lf_48, lf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * ld0_11[k]
                   - f_2 * ld1_11[k]
                   + pb_z[k] * lf_48[k];

        t_105[k] = pa_z[k] * kg_19[k];

        t_106[k] = pa_z[k] * kg_20[k];

        t_107[k] = f_5 * kf_17[k]
                   + pb_z[k] * lf_49[k];

        t_108[k] = pa_z[k] * kg_21[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, kf_18, kf_24, kf_47, \
                         kg_22, kg_23, lf_50, lf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * kf_24[k]
                   + pb_y[k] * lf_50[k];

        t_110[k] = f_6 * kf_18[k]
                   + pa_z[k] * kg_22[k];

        t_111[k] = pa_z[k] * kg_23[k];

        t_112[k] = f_17 * kf_47[k]
                   + pb_x[k] * lf_52[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, kf_20, kf_48, kf_49, \
                         kg_24, lf_51, lf_53, lf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_17 * kf_48[k]
                   + pb_x[k] * lf_53[k];

        t_114[k] = f_17 * kf_49[k]
                   + pb_x[k] * lf_54[k];

        t_115[k] = pa_z[k] * kg_24[k];

        t_116[k] = f_5 * kf_20[k]
                   + pb_z[k] * lf_51[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, kf_21, kf_23, kf_28, \
                         kg_25, kg_26, kg_27, lf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * kf_21[k]
                   + pa_z[k] * kg_25[k];

        t_118[k] = f_6 * kf_28[k]
                   + pb_y[k] * lf_54[k];

        t_119[k] = f_8 * kf_23[k]
                   + pa_z[k] * kg_26[k];

        t_120[k] = pa_y[k] * kg_27[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, kf_29, kf_30, kf_31, \
                         kg_28, kg_29, kg_30, lf_55, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * kf_29[k]
                   + pb_y[k] * lf_55[k];

        t_122[k] = pa_y[k] * kg_28[k];

        t_123[k] = f_6 * kf_30[k]
                   + pa_y[k] * kg_29[k];

        t_124[k] = f_5 * kf_31[k]
                   + pb_y[k] * lf_56[k];

        t_125[k] = pa_y[k] * kg_30[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, kf_33, kf_52, kf_53, \
                         kf_54, kg_31, kg_32, lf_57, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_17 * kf_52[k]
                   + pb_x[k] * lf_57[k];

        t_127[k] = f_17 * kf_53[k]
                   + pb_x[k] * lf_58[k];

        t_128[k] = f_17 * kf_54[k]
                   + pb_x[k] * lf_59[k];

        t_129[k] = pa_y[k] * kg_31[k];

        t_130[k] = f_8 * kf_33[k]
                   + pa_y[k] * kg_32[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, kf_25, kf_35, kf_36, \
                         kg_33, kg_34, lf_57, lf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * kf_25[k]
                   + pb_z[k] * lf_57[k];

        t_132[k] = f_6 * kf_35[k]
                   + pa_y[k] * kg_33[k];

        t_133[k] = f_5 * kf_36[k]
                   + pb_y[k] * lf_60[k];

        t_134[k] = pa_y[k] * kg_34[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, ig0_2, ig1_2, kf_29, \
                         kg_27, ld0_12, ld1_12, lf_61, lf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_14 * ig0_2[k]
                   - f_15 * ig1_2[k]
                   + pa_z[k] * kg_27[k];

        t_136[k] = pb_y[k] * lf_61[k];

        t_137[k] = f_16 * kf_29[k]
                   + pb_z[k] * lf_61[k];

        t_138[k] = f_3 * ld0_12[k]
                   - f_4 * ld1_12[k]
                   + pb_y[k] * lf_62[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, kf_59, kf_60, kf_61, \
                         ld0_14, ld1_14, lf_63, lf_64, lf_65, lf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * lf_63[k];

        t_140[k] = f_17 * kf_59[k]
                   + f_3 * ld0_14[k]
                   - f_4 * ld1_14[k]
                   + pb_x[k] * lf_64[k];

        t_141[k] = f_17 * kf_60[k]
                   + pb_x[k] * lf_65[k];

        t_142[k] = f_17 * kf_61[k]
                   + pb_x[k] * lf_66[k];

        t_143[k] = pb_y[k] * lf_64[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, kf_33, kf_63, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_65, lf_67, lf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_17 * kf_63[k]
                   + pb_x[k] * lf_68[k];

        t_145[k] = f_1 * ld0_13[k]
                   - f_2 * ld1_13[k]
                   + pb_y[k] * lf_65[k];

        t_146[k] = f_16 * kf_33[k]
                   + pb_z[k] * lf_65[k];

        t_147[k] = f_3 * ld0_14[k]
                   - f_4 * ld1_14[k]
                   + pb_y[k] * lf_67[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pa_y, pb_y, ig0_3, ig0_17, ig1_3, \
                         ig1_17, kf_37, kg_35, kg_53, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * lf_68[k];

        t_149[k] = f_18 * ig0_17[k]
                   - f_19 * ig1_17[k]
                   + pa_x[k] * kg_53[k];

        t_150[k] = f_20 * ig0_3[k]
                   - f_21 * ig1_3[k]
                   + pa_y[k] * kg_35[k];

        t_151[k] = f_8 * kf_37[k]
                   + pb_y[k] * lf_69[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_z, kf_66, ld0_15, ld0_16, \
                         ld1_15, ld1_16, lf_69, lf_70, lf_71, lf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_z[k] * lf_69[k];

        t_153[k] = f_8 * kf_66[k]
                   + f_3 * ld0_16[k]
                   - f_4 * ld1_16[k]
                   + pb_x[k] * lf_72[k];

        t_154[k] = pb_z[k] * lf_70[k];

        t_155[k] = f_3 * ld0_15[k]
                   - f_4 * ld1_15[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pb_z, kf_67, kf_69, kf_70, lf_72, \
                         lf_73, lf_75, lf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_8 * kf_67[k]
                   + pb_x[k] * lf_73[k];

        t_157[k] = pb_z[k] * lf_72[k];

        t_158[k] = f_8 * kf_69[k]
                   + pb_x[k] * lf_75[k];

        t_159[k] = f_8 * kf_70[k]
                   + pb_x[k] * lf_76[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_y, pb_z, ig0_20, ig1_20, kf_43, \
                         kg_59, ld0_16, ld1_16, lf_73, lf_74, lf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_20 * ig0_20[k]
                   - f_21 * ig1_20[k]
                   + pa_x[k] * kg_59[k];

        t_161[k] = pb_z[k] * lf_73[k];

        t_162[k] = f_3 * ld0_16[k]
                   - f_4 * ld1_16[k]
                   + pb_z[k] * lf_74[k];

        t_163[k] = f_8 * kf_43[k]
                   + pb_y[k] * lf_76[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_z, pb_z, kf_37, kg_35, kg_36, \
                         kg_37, ld0_17, ld1_17, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * ld0_17[k]
                   - f_2 * ld1_17[k]
                   + pb_z[k] * lf_76[k];

        t_165[k] = pa_z[k] * kg_35[k];

        t_166[k] = pa_z[k] * kg_36[k];

        t_167[k] = f_5 * kf_37[k]
                   + pb_z[k] * lf_77[k];

        t_168[k] = pa_z[k] * kg_37[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_z, pb_x, pb_y, kf_38, kf_45, kf_74, \
                         kg_38, kg_39, lf_78, lf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_16 * kf_45[k]
                   + pb_y[k] * lf_78[k];

        t_170[k] = f_6 * kf_38[k]
                   + pa_z[k] * kg_38[k];

        t_171[k] = pa_z[k] * kg_39[k];

        t_172[k] = f_8 * kf_74[k]
                   + pb_x[k] * lf_80[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pb_x, pb_z, kf_40, kf_75, kf_76, \
                         kg_40, lf_79, lf_81, lf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_8 * kf_75[k]
                   + pb_x[k] * lf_81[k];

        t_174[k] = f_8 * kf_76[k]
                   + pb_x[k] * lf_82[k];

        t_175[k] = pa_z[k] * kg_40[k];

        t_176[k] = f_5 * kf_40[k]
                   + pb_z[k] * lf_79[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pa_z, pb_y, ig0_6, ig1_6, kf_41, \
                         kf_43, kf_49, kg_41, kg_42, kg_44, lf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * kf_41[k]
                   + pa_z[k] * kg_41[k];

        t_178[k] = f_16 * kf_49[k]
                   + pb_y[k] * lf_82[k];

        t_179[k] = f_8 * kf_43[k]
                   + pa_z[k] * kg_42[k];

        t_180[k] = f_9 * ig0_6[k]
                   - f_10 * ig1_6[k]
                   + pa_y[k] * kg_44[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_z, pb_y, pb_z, ig0_4, ig1_4, kf_44, \
                         kf_50, kf_51, kg_43, lf_83, lf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * kf_50[k]
                   + pb_y[k] * lf_83[k];

        t_182[k] = f_6 * kf_44[k]
                   + pb_z[k] * lf_83[k];

        t_183[k] = f_9 * ig0_4[k]
                   - f_10 * ig1_4[k]
                   + pa_z[k] * kg_43[k];

        t_184[k] = f_6 * kf_51[k]
                   + pb_y[k] * lf_84[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pb_x, ig0_7, ig1_7, kf_79, kf_80, \
                         kf_81, kg_45, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_9 * ig0_7[k]
                   - f_10 * ig1_7[k]
                   + pa_y[k] * kg_45[k];

        t_186[k] = f_8 * kf_79[k]
                   + pb_x[k] * lf_85[k];

        t_187[k] = f_8 * kf_80[k]
                   + pb_x[k] * lf_86[k];

        t_188[k] = f_8 * kf_81[k]
                   + pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_x, pb_z, ig0_25, ig1_25, kf_46, kf_82, \
                         kg_66, lf_85, lf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_8 * kf_82[k]
                   + pb_x[k] * lf_88[k];

        t_190[k] = f_20 * ig0_25[k]
                   - f_21 * ig1_25[k]
                   + pa_x[k] * kg_66[k];

        t_191[k] = f_6 * kf_46[k]
                   + pb_z[k] * lf_85[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_x, pa_y, pb_y, ig0_26, ig0_27, ig1_26, \
                         ig1_27, kf_55, kg_46, kg_67, kg_68, lf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_20 * ig0_26[k]
                   - f_21 * ig1_26[k]
                   + pa_x[k] * kg_67[k];

        t_193[k] = f_6 * kf_55[k]
                   + pb_y[k] * lf_88[k];

        t_194[k] = f_20 * ig0_27[k]
                   - f_21 * ig1_27[k]
                   + pa_x[k] * kg_68[k];

        t_195[k] = pa_y[k] * kg_46[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, kf_56, kf_57, kf_58, \
                         kg_47, kg_48, kg_49, lf_89, lf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * kf_56[k]
                   + pb_y[k] * lf_89[k];

        t_197[k] = pa_y[k] * kg_47[k];

        t_198[k] = f_6 * kf_57[k]
                   + pa_y[k] * kg_48[k];

        t_199[k] = f_5 * kf_58[k]
                   + pb_y[k] * lf_90[k];

        t_200[k] = pa_y[k] * kg_49[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pa_y, pb_x, kf_60, kf_85, kf_86, \
                         kf_87, kg_50, kg_51, lf_91, lf_92, lf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_8 * kf_85[k]
                   + pb_x[k] * lf_91[k];

        t_202[k] = f_8 * kf_86[k]
                   + pb_x[k] * lf_92[k];

        t_203[k] = f_8 * kf_87[k]
                   + pb_x[k] * lf_93[k];

        t_204[k] = pa_y[k] * kg_50[k];

        t_205[k] = f_8 * kf_60[k]
                   + pa_y[k] * kg_51[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_y, pb_y, pb_z, kf_52, kf_62, kf_63, \
                         kg_52, kg_53, lf_91, lf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_16 * kf_52[k]
                   + pb_z[k] * lf_91[k];

        t_207[k] = f_6 * kf_62[k]
                   + pa_y[k] * kg_52[k];

        t_208[k] = f_5 * kf_63[k]
                   + pb_y[k] * lf_94[k];

        t_209[k] = pa_y[k] * kg_53[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_z, pb_y, pb_z, ig0_6, ig1_6, kf_56, \
                         kg_46, ld0_18, ld1_18, lf_95, lf_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_20 * ig0_6[k]
                   - f_21 * ig1_6[k]
                   + pa_z[k] * kg_46[k];

        t_211[k] = pb_y[k] * lf_95[k];

        t_212[k] = f_8 * kf_56[k]
                   + pb_z[k] * lf_95[k];

        t_213[k] = f_3 * ld0_18[k]
                   - f_4 * ld1_18[k]
                   + pb_y[k] * lf_96[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, kf_92, kf_93, kf_94, \
                         ld0_20, ld1_20, lf_97, lf_98, lf_99, lf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * lf_97[k];

        t_215[k] = f_8 * kf_92[k]
                   + f_3 * ld0_20[k]
                   - f_4 * ld1_20[k]
                   + pb_x[k] * lf_98[k];

        t_216[k] = f_8 * kf_93[k]
                   + pb_x[k] * lf_99[k];

        t_217[k] = f_8 * kf_94[k]
                   + pb_x[k] * lf_100[k];

        t_218[k] = pb_y[k] * lf_98[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pb_y, pb_z, kf_60, kf_96, ld0_19, \
                         ld0_20, ld1_19, ld1_20, lf_99, lf_101, \
                         lf_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_8 * kf_96[k]
                   + pb_x[k] * lf_102[k];

        t_220[k] = f_1 * ld0_19[k]
                   - f_2 * ld1_19[k]
                   + pb_y[k] * lf_99[k];

        t_221[k] = f_8 * kf_60[k]
                   + pb_z[k] * lf_99[k];

        t_222[k] = f_3 * ld0_20[k]
                   - f_4 * ld1_20[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pa_y, pb_y, ig0_9, ig0_32, ig1_9, \
                         ig1_32, kf_64, kg_54, kg_78, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_y[k] * lf_102[k];

        t_224[k] = f_20 * ig0_32[k]
                   - f_21 * ig1_32[k]
                   + pa_x[k] * kg_78[k];

        t_225[k] = f_18 * ig0_9[k]
                   - f_19 * ig1_9[k]
                   + pa_y[k] * kg_54[k];

        t_226[k] = f_17 * kf_64[k]
                   + pb_y[k] * lf_103[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pb_x, pb_z, kf_99, ld0_21, ld0_22, \
                         ld1_21, ld1_22, lf_103, lf_104, lf_105, \
                         lf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_z[k] * lf_103[k];

        t_228[k] = f_16 * kf_99[k]
                   + f_3 * ld0_22[k]
                   - f_4 * ld1_22[k]
                   + pb_x[k] * lf_106[k];

        t_229[k] = pb_z[k] * lf_104[k];

        t_230[k] = f_3 * ld0_21[k]
                   - f_4 * ld1_21[k]
                   + pb_z[k] * lf_105[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pb_z, kf_100, kf_102, kf_103, \
                         lf_106, lf_107, lf_109, lf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_16 * kf_100[k]
                   + pb_x[k] * lf_107[k];

        t_232[k] = pb_z[k] * lf_106[k];

        t_233[k] = f_16 * kf_102[k]
                   + pb_x[k] * lf_109[k];

        t_234[k] = f_16 * kf_103[k]
                   + pb_x[k] * lf_110[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pa_x, pb_y, pb_z, ig0_33, ig1_33, kf_70, \
                         kg_84, ld0_22, ld1_22, lf_107, lf_108, \
                         lf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_14 * ig0_33[k]
                   - f_15 * ig1_33[k]
                   + pa_x[k] * kg_84[k];

        t_236[k] = pb_z[k] * lf_107[k];

        t_237[k] = f_3 * ld0_22[k]
                   - f_4 * ld1_22[k]
                   + pb_z[k] * lf_108[k];

        t_238[k] = f_17 * kf_70[k]
                   + pb_y[k] * lf_110[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_z, pb_z, kf_64, kg_54, kg_55, \
                         kg_56, ld0_23, ld1_23, lf_110, lf_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * ld0_23[k]
                   - f_2 * ld1_23[k]
                   + pb_z[k] * lf_110[k];

        t_240[k] = pa_z[k] * kg_54[k];

        t_241[k] = pa_z[k] * kg_55[k];

        t_242[k] = f_5 * kf_64[k]
                   + pb_z[k] * lf_111[k];

        t_243[k] = pa_z[k] * kg_56[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pb_x, pb_y, kf_65, kf_72, kf_107, \
                         kg_57, kg_58, lf_112, lf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_8 * kf_72[k]
                   + pb_y[k] * lf_112[k];

        t_245[k] = f_6 * kf_65[k]
                   + pa_z[k] * kg_57[k];

        t_246[k] = pa_z[k] * kg_58[k];

        t_247[k] = f_16 * kf_107[k]
                   + pb_x[k] * lf_114[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_z, pb_x, pb_z, kf_67, kf_108, kf_109, \
                         kg_59, lf_113, lf_115, lf_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_16 * kf_108[k]
                   + pb_x[k] * lf_115[k];

        t_249[k] = f_16 * kf_109[k]
                   + pb_x[k] * lf_116[k];

        t_250[k] = pa_z[k] * kg_59[k];

        t_251[k] = f_5 * kf_67[k]
                   + pb_z[k] * lf_113[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_y, pa_z, pb_y, ig0_13, ig1_13, kf_68, \
                         kf_70, kf_76, kg_60, kg_61, kg_63, lf_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_6 * kf_68[k]
                   + pa_z[k] * kg_60[k];

        t_253[k] = f_8 * kf_76[k]
                   + pb_y[k] * lf_116[k];

        t_254[k] = f_8 * kf_70[k]
                   + pa_z[k] * kg_61[k];

        t_255[k] = f_14 * ig0_13[k]
                   - f_15 * ig1_13[k]
                   + pa_y[k] * kg_63[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_z, pb_y, pb_z, ig0_10, ig1_10, kf_71, \
                         kf_77, kf_78, kg_62, lf_117, lf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_16 * kf_77[k]
                   + pb_y[k] * lf_117[k];

        t_257[k] = f_6 * kf_71[k]
                   + pb_z[k] * lf_117[k];

        t_258[k] = f_9 * ig0_10[k]
                   - f_10 * ig1_10[k]
                   + pa_z[k] * kg_62[k];

        t_259[k] = f_16 * kf_78[k]
                   + pb_y[k] * lf_118[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pb_x, ig0_14, ig1_14, kf_112, \
                         kf_113, kf_114, kg_65, lf_119, lf_120, \
                         lf_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * ig0_14[k]
                   - f_15 * ig1_14[k]
                   + pa_y[k] * kg_65[k];

        t_261[k] = f_16 * kf_112[k]
                   + pb_x[k] * lf_119[k];

        t_262[k] = f_16 * kf_113[k]
                   + pb_x[k] * lf_120[k];

        t_263[k] = f_16 * kf_114[k]
                   + pb_x[k] * lf_121[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_x, pb_x, pb_z, ig0_34, ig1_34, kf_73, kf_115, \
                         kg_91, lf_119, lf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_16 * kf_115[k]
                   + pb_x[k] * lf_122[k];

        t_265[k] = f_14 * ig0_34[k]
                   - f_15 * ig1_34[k]
                   + pa_x[k] * kg_91[k];

        t_266[k] = f_6 * kf_73[k]
                   + pb_z[k] * lf_119[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_x, pb_y, ig0_35, ig0_36, ig1_35, ig1_36, \
                         kf_82, kg_92, kg_93, lf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * ig0_35[k]
                   - f_15 * ig1_35[k]
                   + pa_x[k] * kg_92[k];

        t_268[k] = f_16 * kf_82[k]
                   + pb_y[k] * lf_122[k];

        t_269[k] = f_14 * ig0_36[k]
                   - f_15 * ig1_36[k]
                   + pa_x[k] * kg_93[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_y, pb_y, pb_z, ig0_15, ig1_15, kf_77, kf_83, \
                         kg_69, lf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_9 * ig0_15[k]
                   - f_10 * ig1_15[k]
                   + pa_y[k] * kg_69[k];

        t_271[k] = f_6 * kf_83[k]
                   + pb_y[k] * lf_123[k];

        t_272[k] = f_16 * kf_77[k]
                   + pb_z[k] * lf_123[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_y, pa_z, pb_y, ig0_12, ig0_16, ig1_12, \
                         ig1_16, kf_84, kg_64, kg_70, lf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_14 * ig0_12[k]
                   - f_15 * ig1_12[k]
                   + pa_z[k] * kg_64[k];

        t_274[k] = f_6 * kf_84[k]
                   + pb_y[k] * lf_124[k];

        t_275[k] = f_9 * ig0_16[k]
                   - f_10 * ig1_16[k]
                   + pa_y[k] * kg_70[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, kf_118, kf_119, kf_120, kf_121, \
                         lf_125, lf_126, lf_127, lf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_16 * kf_118[k]
                   + pb_x[k] * lf_125[k];

        t_277[k] = f_16 * kf_119[k]
                   + pb_x[k] * lf_126[k];

        t_278[k] = f_16 * kf_120[k]
                   + pb_x[k] * lf_127[k];

        t_279[k] = f_16 * kf_121[k]
                   + pb_x[k] * lf_128[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_x, pb_z, ig0_37, ig0_38, ig1_37, ig1_38, \
                         kf_79, kg_97, kg_98, lf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_14 * ig0_37[k]
                   - f_15 * ig1_37[k]
                   + pa_x[k] * kg_97[k];

        t_281[k] = f_16 * kf_79[k]
                   + pb_z[k] * lf_125[k];

        t_282[k] = f_14 * ig0_38[k]
                   - f_15 * ig1_38[k]
                   + pa_x[k] * kg_98[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pb_y, ig0_39, ig1_39, kf_88, \
                         kf_89, kg_71, kg_99, lf_128, lf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_6 * kf_88[k]
                   + pb_y[k] * lf_128[k];

        t_284[k] = f_14 * ig0_39[k]
                   - f_15 * ig1_39[k]
                   + pa_x[k] * kg_99[k];

        t_285[k] = pa_y[k] * kg_71[k];

        t_286[k] = f_5 * kf_89[k]
                   + pb_y[k] * lf_129[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_y, pb_x, pb_y, kf_90, kf_91, \
                         kf_124, kg_72, kg_73, kg_74, lf_130, lf_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_y[k] * kg_72[k];

        t_288[k] = f_6 * kf_90[k]
                   + pa_y[k] * kg_73[k];

        t_289[k] = f_5 * kf_91[k]
                   + pb_y[k] * lf_130[k];

        t_290[k] = pa_y[k] * kg_74[k];

        t_291[k] = f_16 * kf_124[k]
                   + pb_x[k] * lf_131[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_x, kf_93, kf_125, kf_126, kg_75, \
                         kg_76, lf_132, lf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_16 * kf_125[k]
                   + pb_x[k] * lf_132[k];

        t_293[k] = f_16 * kf_126[k]
                   + pb_x[k] * lf_133[k];

        t_294[k] = pa_y[k] * kg_75[k];

        t_295[k] = f_8 * kf_93[k]
                   + pa_y[k] * kg_76[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_y, pb_z, kf_85, kf_95, kf_96, \
                         kg_77, kg_78, lf_131, lf_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_8 * kf_85[k]
                   + pb_z[k] * lf_131[k];

        t_297[k] = f_6 * kf_95[k]
                   + pa_y[k] * kg_77[k];

        t_298[k] = f_5 * kf_96[k]
                   + pb_y[k] * lf_134[k];

        t_299[k] = pa_y[k] * kg_78[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_z, pb_y, pb_z, ig0_15, ig1_15, kf_89, \
                         kg_71, ld0_24, ld1_24, lf_135, lf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_18 * ig0_15[k]
                   - f_19 * ig1_15[k]
                   + pa_z[k] * kg_71[k];

        t_301[k] = pb_y[k] * lf_135[k];

        t_302[k] = f_17 * kf_89[k]
                   + pb_z[k] * lf_135[k];

        t_303[k] = f_3 * ld0_24[k]
                   - f_4 * ld1_24[k]
                   + pb_y[k] * lf_136[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, kf_131, kf_132, \
                         kf_133, ld0_26, ld1_26, lf_137, lf_138, lf_139, \
                         lf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pb_y[k] * lf_137[k];

        t_305[k] = f_16 * kf_131[k]
                   + f_3 * ld0_26[k]
                   - f_4 * ld1_26[k]
                   + pb_x[k] * lf_138[k];

        t_306[k] = f_16 * kf_132[k]
                   + pb_x[k] * lf_139[k];

        t_307[k] = f_16 * kf_133[k]
                   + pb_x[k] * lf_140[k];

        t_308[k] = pb_y[k] * lf_138[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pb_y, pb_z, kf_93, kf_135, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_139, lf_141, \
                         lf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_16 * kf_135[k]
                   + pb_x[k] * lf_142[k];

        t_310[k] = f_1 * ld0_25[k]
                   - f_2 * ld1_25[k]
                   + pb_y[k] * lf_139[k];

        t_311[k] = f_17 * kf_93[k]
                   + pb_z[k] * lf_139[k];

        t_312[k] = f_3 * ld0_26[k]
                   - f_4 * ld1_26[k]
                   + pb_y[k] * lf_141[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_x, pa_y, pb_y, ig0_18, ig0_40, ig1_18, \
                         ig1_40, kf_97, kg_79, kg_109, lf_142, lf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * lf_142[k];

        t_314[k] = f_14 * ig0_40[k]
                   - f_15 * ig1_40[k]
                   + pa_x[k] * kg_109[k];

        t_315[k] = f_12 * ig0_18[k]
                   - f_13 * ig1_18[k]
                   + pa_y[k] * kg_79[k];

        t_316[k] = f_11 * kf_97[k]
                   + pb_y[k] * lf_143[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_x, pb_z, kf_137, ld0_27, ld0_28, \
                         ld1_27, ld1_28, lf_143, lf_144, lf_145, \
                         lf_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pb_z[k] * lf_143[k];

        t_318[k] = f_6 * kf_137[k]
                   + f_3 * ld0_28[k]
                   - f_4 * ld1_28[k]
                   + pb_x[k] * lf_146[k];

        t_319[k] = pb_z[k] * lf_144[k];

        t_320[k] = f_3 * ld0_27[k]
                   - f_4 * ld1_27[k]
                   + pb_z[k] * lf_145[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_z, kf_138, kf_139, kf_140, \
                         lf_146, lf_147, lf_149, lf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_6 * kf_138[k]
                   + pb_x[k] * lf_147[k];

        t_322[k] = pb_z[k] * lf_146[k];

        t_323[k] = f_6 * kf_139[k]
                   + pb_x[k] * lf_149[k];

        t_324[k] = f_6 * kf_140[k]
                   + pb_x[k] * lf_150[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_x, pb_y, pb_z, ig0_41, ig1_41, kf_103, \
                         kg_114, ld0_28, ld1_28, lf_147, lf_148, \
                         lf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_9 * ig0_41[k]
                   - f_10 * ig1_41[k]
                   + pa_x[k] * kg_114[k];

        t_326[k] = pb_z[k] * lf_147[k];

        t_327[k] = f_3 * ld0_28[k]
                   - f_4 * ld1_28[k]
                   + pb_z[k] * lf_148[k];

        t_328[k] = f_11 * kf_103[k]
                   + pb_y[k] * lf_150[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, pa_z, pb_z, kf_97, kg_79, kg_80, \
                         kg_81, ld0_29, ld1_29, lf_150, lf_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * ld0_29[k]
                   - f_2 * ld1_29[k]
                   + pb_z[k] * lf_150[k];

        t_330[k] = pa_z[k] * kg_79[k];

        t_331[k] = pa_z[k] * kg_80[k];

        t_332[k] = f_5 * kf_97[k]
                   + pb_z[k] * lf_151[k];

        t_333[k] = pa_z[k] * kg_81[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_z, pb_x, pb_y, kf_98, kf_105, kf_143, \
                         kg_82, kg_83, lf_152, lf_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_17 * kf_105[k]
                   + pb_y[k] * lf_152[k];

        t_335[k] = f_6 * kf_98[k]
                   + pa_z[k] * kg_82[k];

        t_336[k] = pa_z[k] * kg_83[k];

        t_337[k] = f_6 * kf_143[k]
                   + pb_x[k] * lf_154[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_z, pb_x, pb_z, kf_100, kf_144, kf_145, \
                         kg_84, lf_153, lf_155, lf_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_6 * kf_144[k]
                   + pb_x[k] * lf_155[k];

        t_339[k] = f_6 * kf_145[k]
                   + pb_x[k] * lf_156[k];

        t_340[k] = pa_z[k] * kg_84[k];

        t_341[k] = f_5 * kf_100[k]
                   + pb_z[k] * lf_153[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pa_y, pa_z, pb_y, ig0_22, ig1_22, kf_101, \
                         kf_103, kf_109, kg_85, kg_86, kg_88, lf_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_6 * kf_101[k]
                   + pa_z[k] * kg_85[k];

        t_343[k] = f_17 * kf_109[k]
                   + pb_y[k] * lf_156[k];

        t_344[k] = f_8 * kf_103[k]
                   + pa_z[k] * kg_86[k];

        t_345[k] = f_20 * ig0_22[k]
                   - f_21 * ig1_22[k]
                   + pa_y[k] * kg_88[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_z, pb_y, pb_z, ig0_19, ig1_19, kf_104, \
                         kf_110, kf_111, kg_87, lf_157, lf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_8 * kf_110[k]
                   + pb_y[k] * lf_157[k];

        t_347[k] = f_6 * kf_104[k]
                   + pb_z[k] * lf_157[k];

        t_348[k] = f_9 * ig0_19[k]
                   - f_10 * ig1_19[k]
                   + pa_z[k] * kg_87[k];

        t_349[k] = f_8 * kf_111[k]
                   + pb_y[k] * lf_158[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_y, pb_x, ig0_24, ig1_24, kf_148, \
                         kf_149, kf_150, kg_90, lf_159, lf_160, \
                         lf_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_20 * ig0_24[k]
                   - f_21 * ig1_24[k]
                   + pa_y[k] * kg_90[k];

        t_351[k] = f_6 * kf_148[k]
                   + pb_x[k] * lf_159[k];

        t_352[k] = f_6 * kf_149[k]
                   + pb_x[k] * lf_160[k];

        t_353[k] = f_6 * kf_150[k]
                   + pb_x[k] * lf_161[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_x, pb_x, pb_z, ig0_43, ig1_43, kf_106, \
                         kf_151, kg_115, lf_159, lf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_6 * kf_151[k]
                   + pb_x[k] * lf_162[k];

        t_355[k] = f_9 * ig0_43[k]
                   - f_10 * ig1_43[k]
                   + pa_x[k] * kg_115[k];

        t_356[k] = f_6 * kf_106[k]
                   + pb_z[k] * lf_159[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_x, pb_y, ig0_44, ig0_45, ig1_44, ig1_45, \
                         kf_115, kg_116, kg_117, lf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_9 * ig0_44[k]
                   - f_10 * ig1_44[k]
                   + pa_x[k] * kg_116[k];

        t_358[k] = f_8 * kf_115[k]
                   + pb_y[k] * lf_162[k];

        t_359[k] = f_9 * ig0_45[k]
                   - f_10 * ig1_45[k]
                   + pa_x[k] * kg_117[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pb_y, pb_z, ig0_28, ig1_28, kf_110, \
                         kf_116, kg_94, lf_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_14 * ig0_28[k]
                   - f_15 * ig1_28[k]
                   + pa_y[k] * kg_94[k];

        t_361[k] = f_16 * kf_116[k]
                   + pb_y[k] * lf_163[k];

        t_362[k] = f_16 * kf_110[k]
                   + pb_z[k] * lf_163[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_y, pa_z, pb_y, ig0_21, ig0_29, ig1_21, \
                         ig1_29, kf_117, kg_89, kg_96, lf_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * ig0_21[k]
                   - f_15 * ig1_21[k]
                   + pa_z[k] * kg_89[k];

        t_364[k] = f_16 * kf_117[k]
                   + pb_y[k] * lf_164[k];

        t_365[k] = f_14 * ig0_29[k]
                   - f_15 * ig1_29[k]
                   + pa_y[k] * kg_96[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, kf_154, kf_155, kf_156, kf_157, \
                         lf_165, lf_166, lf_167, lf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_6 * kf_154[k]
                   + pb_x[k] * lf_165[k];

        t_367[k] = f_6 * kf_155[k]
                   + pb_x[k] * lf_166[k];

        t_368[k] = f_6 * kf_156[k]
                   + pb_x[k] * lf_167[k];

        t_369[k] = f_6 * kf_157[k]
                   + pb_x[k] * lf_168[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pa_x, pb_z, ig0_46, ig0_47, ig1_46, ig1_47, \
                         kf_112, kg_118, kg_119, lf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_9 * ig0_46[k]
                   - f_10 * ig1_46[k]
                   + pa_x[k] * kg_118[k];

        t_371[k] = f_16 * kf_112[k]
                   + pb_z[k] * lf_165[k];

        t_372[k] = f_9 * ig0_47[k]
                   - f_10 * ig1_47[k]
                   + pa_x[k] * kg_119[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pa_x, pa_y, pb_y, ig0_30, ig0_48, ig1_30, \
                         ig1_48, kf_121, kg_100, kg_120, lf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * kf_121[k]
                   + pb_y[k] * lf_168[k];

        t_374[k] = f_9 * ig0_48[k]
                   - f_10 * ig1_48[k]
                   + pa_x[k] * kg_120[k];

        t_375[k] = f_9 * ig0_30[k]
                   - f_10 * ig1_30[k]
                   + pa_y[k] * kg_100[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pb_y, pb_z, ig0_23, ig1_23, kf_116, \
                         kf_122, kf_123, kg_95, lf_169, lf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_6 * kf_122[k]
                   + pb_y[k] * lf_169[k];

        t_377[k] = f_8 * kf_116[k]
                   + pb_z[k] * lf_169[k];

        t_378[k] = f_20 * ig0_23[k]
                   - f_21 * ig1_23[k]
                   + pa_z[k] * kg_95[k];

        t_379[k] = f_6 * kf_123[k]
                   + pb_y[k] * lf_170[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pa_y, pb_x, ig0_31, ig1_31, kf_160, \
                         kf_161, kf_162, kg_101, lf_171, lf_172, \
                         lf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_9 * ig0_31[k]
                   - f_10 * ig1_31[k]
                   + pa_y[k] * kg_101[k];

        t_381[k] = f_6 * kf_160[k]
                   + pb_x[k] * lf_171[k];

        t_382[k] = f_6 * kf_161[k]
                   + pb_x[k] * lf_172[k];

        t_383[k] = f_6 * kf_162[k]
                   + pb_x[k] * lf_173[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pa_x, pb_x, pb_z, ig0_49, ig1_49, kf_118, \
                         kf_163, kg_121, lf_171, lf_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_6 * kf_163[k]
                   + pb_x[k] * lf_174[k];

        t_385[k] = f_9 * ig0_49[k]
                   - f_10 * ig1_49[k]
                   + pa_x[k] * kg_121[k];

        t_386[k] = f_8 * kf_118[k]
                   + pb_z[k] * lf_171[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_x, pa_y, pb_y, ig0_50, ig0_51, ig1_50, \
                         ig1_51, kf_127, kg_102, kg_122, kg_123, \
                         lf_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_9 * ig0_50[k]
                   - f_10 * ig1_50[k]
                   + pa_x[k] * kg_122[k];

        t_388[k] = f_6 * kf_127[k]
                   + pb_y[k] * lf_174[k];

        t_389[k] = f_9 * ig0_51[k]
                   - f_10 * ig1_51[k]
                   + pa_x[k] * kg_123[k];

        t_390[k] = pa_y[k] * kg_102[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, pa_y, pb_y, kf_128, kf_129, \
                         kf_130, kg_103, kg_104, kg_105, lf_175, \
                         lf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_5 * kf_128[k]
                   + pb_y[k] * lf_175[k];

        t_392[k] = pa_y[k] * kg_103[k];

        t_393[k] = f_6 * kf_129[k]
                   + pa_y[k] * kg_104[k];

        t_394[k] = f_5 * kf_130[k]
                   + pb_y[k] * lf_176[k];

        t_395[k] = pa_y[k] * kg_105[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, pa_y, pb_x, kf_132, kf_166, \
                         kf_167, kf_168, kg_106, kg_107, lf_177, lf_178, \
                         lf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_6 * kf_166[k]
                   + pb_x[k] * lf_177[k];

        t_397[k] = f_6 * kf_167[k]
                   + pb_x[k] * lf_178[k];

        t_398[k] = f_6 * kf_168[k]
                   + pb_x[k] * lf_179[k];

        t_399[k] = pa_y[k] * kg_106[k];

        t_400[k] = f_8 * kf_132[k]
                   + pa_y[k] * kg_107[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pb_y, pb_z, kf_124, kf_134, kf_135, \
                         kg_108, kg_109, lf_177, lf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_17 * kf_124[k]
                   + pb_z[k] * lf_177[k];

        t_402[k] = f_6 * kf_134[k]
                   + pa_y[k] * kg_108[k];

        t_403[k] = f_5 * kf_135[k]
                   + pb_y[k] * lf_180[k];

        t_404[k] = pa_y[k] * kg_109[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_z, pb_y, pb_z, ig0_30, ig1_30, kf_128, \
                         kg_102, ld0_30, ld1_30, lf_181, lf_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_12 * ig0_30[k]
                   - f_13 * ig1_30[k]
                   + pa_z[k] * kg_102[k];

        t_406[k] = pb_y[k] * lf_181[k];

        t_407[k] = f_11 * kf_128[k]
                   + pb_z[k] * lf_181[k];

        t_408[k] = f_3 * ld0_30[k]
                   - f_4 * ld1_30[k]
                   + pb_y[k] * lf_182[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pb_x, pb_y, kf_171, kf_172, \
                         kf_173, ld0_32, ld1_32, lf_183, lf_184, lf_185, \
                         lf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * lf_183[k];

        t_410[k] = f_6 * kf_171[k]
                   + f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_x[k] * lf_184[k];

        t_411[k] = f_6 * kf_172[k]
                   + pb_x[k] * lf_185[k];

        t_412[k] = f_6 * kf_173[k]
                   + pb_x[k] * lf_186[k];

        t_413[k] = pb_y[k] * lf_184[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pb_y, pb_z, kf_132, kf_174, ld0_31, \
                         ld0_32, ld1_31, ld1_32, lf_185, lf_187, \
                         lf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_6 * kf_174[k]
                   + pb_x[k] * lf_188[k];

        t_415[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_185[k];

        t_416[k] = f_11 * kf_132[k]
                   + pb_z[k] * lf_185[k];

        t_417[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_187[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pa_x, pb_y, pb_z, ig0_53, ig1_53, \
                         kf_136, kf_175, kg_128, kg_129, lf_188, \
                         lf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = pb_y[k] * lf_188[k];

        t_419[k] = f_9 * ig0_53[k]
                   - f_10 * ig1_53[k]
                   + pa_x[k] * kg_128[k];

        t_420[k] = f_8 * kf_175[k]
                   + pa_x[k] * kg_129[k];

        t_421[k] = f_7 * kf_136[k]
                   + pb_y[k] * lf_189[k];

        t_422[k] = pb_z[k] * lf_189[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, pa_x, pb_x, pb_z, kf_177, kf_178, \
                         kf_179, kg_131, kg_132, lf_190, lf_191, \
                         lf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_6 * kf_177[k]
                   + pa_x[k] * kg_131[k];

        t_424[k] = pb_z[k] * lf_190[k];

        t_425[k] = f_6 * kf_178[k]
                   + pa_x[k] * kg_132[k];

        t_426[k] = f_5 * kf_179[k]
                   + pb_x[k] * lf_192[k];

        t_427[k] = pb_z[k] * lf_191[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pa_x, pb_x, pb_z, kf_181, kf_182, \
                         kg_133, kg_134, lf_192, lf_193, lf_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_5 * kf_181[k]
                   + pb_x[k] * lf_193[k];

        t_429[k] = f_5 * kf_182[k]
                   + pb_x[k] * lf_194[k];

        t_430[k] = pa_x[k] * kg_133[k];

        t_431[k] = pb_z[k] * lf_192[k];

        t_432[k] = pa_x[k] * kg_134[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, t_438, pa_x, pa_z, pb_z, kf_136, \
                         kg_110, kg_111, kg_112, kg_135, kg_136, \
                         lf_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = pa_x[k] * kg_135[k];

        t_434[k] = pa_x[k] * kg_136[k];

        t_435[k] = pa_z[k] * kg_110[k];

        t_436[k] = pa_z[k] * kg_111[k];

        t_437[k] = f_5 * kf_136[k]
                   + pb_z[k] * lf_195[k];

        t_438[k] = pa_z[k] * kg_112[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pa_x, pa_z, pb_x, pb_y, kf_142, kf_185, \
                         kf_187, kg_113, kg_137, lf_196, lf_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_11 * kf_142[k]
                   + pb_y[k] * lf_196[k];

        t_440[k] = f_6 * kf_185[k]
                   + pa_x[k] * kg_137[k];

        t_441[k] = pa_z[k] * kg_113[k];

        t_442[k] = f_5 * kf_187[k]
                   + pb_x[k] * lf_197[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, pa_x, pb_x, kf_188, kf_189, \
                         kg_138, kg_139, kg_140, kg_141, lf_198, \
                         lf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_5 * kf_188[k]
                   + pb_x[k] * lf_198[k];

        t_444[k] = f_5 * kf_189[k]
                   + pb_x[k] * lf_199[k];

        t_445[k] = pa_x[k] * kg_138[k];

        t_446[k] = pa_x[k] * kg_139[k];

        t_447[k] = pa_x[k] * kg_140[k];

        t_448[k] = pa_x[k] * kg_141[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, pa_x, pb_y, pb_z, kf_141, kf_146, \
                         kf_190, kf_192, kg_142, kg_143, kg_144, \
                         lf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_x[k] * kg_142[k];

        t_450[k] = f_8 * kf_190[k]
                   + pa_x[k] * kg_143[k];

        t_451[k] = f_17 * kf_146[k]
                   + pb_y[k] * lf_200[k];

        t_452[k] = f_6 * kf_141[k]
                   + pb_z[k] * lf_200[k];

        t_453[k] = f_6 * kf_192[k]
                   + pa_x[k] * kg_144[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, pa_x, pb_x, pb_y, kf_147, kf_193, kf_194, \
                         kf_195, kg_145, lf_201, lf_202, lf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_17 * kf_147[k]
                   + pb_y[k] * lf_201[k];

        t_455[k] = f_6 * kf_193[k]
                   + pa_x[k] * kg_145[k];

        t_456[k] = f_5 * kf_194[k]
                   + pb_x[k] * lf_202[k];

        t_457[k] = f_5 * kf_195[k]
                   + pb_x[k] * lf_203[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, t_463, pa_x, pb_x, kf_196, kf_197, \
                         kg_146, kg_147, kg_148, kg_149, lf_204, \
                         lf_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_5 * kf_196[k]
                   + pb_x[k] * lf_204[k];

        t_459[k] = f_5 * kf_197[k]
                   + pb_x[k] * lf_205[k];

        t_460[k] = pa_x[k] * kg_146[k];

        t_461[k] = pa_x[k] * kg_147[k];

        t_462[k] = pa_x[k] * kg_148[k];

        t_463[k] = pa_x[k] * kg_149[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, pa_x, pb_y, pb_z, kf_146, kf_152, \
                         kf_198, kf_200, kg_150, kg_151, kg_152, \
                         lf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pa_x[k] * kg_150[k];

        t_465[k] = f_8 * kf_198[k]
                   + pa_x[k] * kg_151[k];

        t_466[k] = f_8 * kf_152[k]
                   + pb_y[k] * lf_206[k];

        t_467[k] = f_16 * kf_146[k]
                   + pb_z[k] * lf_206[k];

        t_468[k] = f_6 * kf_200[k]
                   + pa_x[k] * kg_152[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pa_x, pb_x, pb_y, kf_153, kf_201, kf_202, \
                         kf_203, kg_153, lf_207, lf_208, lf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_8 * kf_153[k]
                   + pb_y[k] * lf_207[k];

        t_470[k] = f_6 * kf_201[k]
                   + pa_x[k] * kg_153[k];

        t_471[k] = f_5 * kf_202[k]
                   + pb_x[k] * lf_208[k];

        t_472[k] = f_5 * kf_203[k]
                   + pb_x[k] * lf_209[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, t_478, pa_x, pb_x, kf_204, kf_205, \
                         kg_154, kg_155, kg_156, kg_157, lf_210, \
                         lf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_5 * kf_204[k]
                   + pb_x[k] * lf_210[k];

        t_474[k] = f_5 * kf_205[k]
                   + pb_x[k] * lf_211[k];

        t_475[k] = pa_x[k] * kg_154[k];

        t_476[k] = pa_x[k] * kg_155[k];

        t_477[k] = pa_x[k] * kg_156[k];

        t_478[k] = pa_x[k] * kg_157[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, pa_x, pb_y, pb_z, kf_152, kf_158, \
                         kf_206, kf_208, kg_158, kg_159, kg_160, \
                         lf_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pa_x[k] * kg_158[k];

        t_480[k] = f_8 * kf_206[k]
                   + pa_x[k] * kg_159[k];

        t_481[k] = f_16 * kf_158[k]
                   + pb_y[k] * lf_212[k];

        t_482[k] = f_8 * kf_152[k]
                   + pb_z[k] * lf_212[k];

        t_483[k] = f_6 * kf_208[k]
                   + pa_x[k] * kg_160[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pa_x, pb_x, pb_y, kf_159, kf_209, kf_210, \
                         kf_211, kg_161, lf_213, lf_214, lf_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_16 * kf_159[k]
                   + pb_y[k] * lf_213[k];

        t_485[k] = f_6 * kf_209[k]
                   + pa_x[k] * kg_161[k];

        t_486[k] = f_5 * kf_210[k]
                   + pb_x[k] * lf_214[k];

        t_487[k] = f_5 * kf_211[k]
                   + pb_x[k] * lf_215[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, t_493, pa_x, pb_x, kf_212, kf_213, \
                         kg_162, kg_163, kg_164, kg_165, lf_216, \
                         lf_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * kf_212[k]
                   + pb_x[k] * lf_216[k];

        t_489[k] = f_5 * kf_213[k]
                   + pb_x[k] * lf_217[k];

        t_490[k] = pa_x[k] * kg_162[k];

        t_491[k] = pa_x[k] * kg_163[k];

        t_492[k] = pa_x[k] * kg_164[k];

        t_493[k] = pa_x[k] * kg_165[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, pa_x, pb_y, pb_z, kf_158, kf_164, \
                         kf_214, kf_216, kg_166, kg_167, kg_168, \
                         lf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = pa_x[k] * kg_166[k];

        t_495[k] = f_8 * kf_214[k]
                   + pa_x[k] * kg_167[k];

        t_496[k] = f_6 * kf_164[k]
                   + pb_y[k] * lf_218[k];

        t_497[k] = f_17 * kf_158[k]
                   + pb_z[k] * lf_218[k];

        t_498[k] = f_6 * kf_216[k]
                   + pa_x[k] * kg_168[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, pa_x, pb_x, pb_y, kf_165, kf_217, kf_218, \
                         kf_219, kg_169, lf_219, lf_220, lf_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_6 * kf_165[k]
                   + pb_y[k] * lf_219[k];

        t_500[k] = f_6 * kf_217[k]
                   + pa_x[k] * kg_169[k];

        t_501[k] = f_5 * kf_218[k]
                   + pb_x[k] * lf_220[k];

        t_502[k] = f_5 * kf_219[k]
                   + pb_x[k] * lf_221[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, t_508, pa_x, pb_x, kf_220, kf_221, \
                         kg_170, kg_171, kg_172, kg_173, lf_222, \
                         lf_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_5 * kf_220[k]
                   + pb_x[k] * lf_222[k];

        t_504[k] = f_5 * kf_221[k]
                   + pb_x[k] * lf_223[k];

        t_505[k] = pa_x[k] * kg_170[k];

        t_506[k] = pa_x[k] * kg_171[k];

        t_507[k] = pa_x[k] * kg_172[k];

        t_508[k] = pa_x[k] * kg_173[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, pa_x, pa_y, pb_y, kf_169, kf_224, \
                         kg_124, kg_125, kg_174, kg_175, lf_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pa_x[k] * kg_174[k];

        t_510[k] = pa_y[k] * kg_124[k];

        t_511[k] = f_5 * kf_169[k]
                   + pb_y[k] * lf_224[k];

        t_512[k] = pa_y[k] * kg_125[k];

        t_513[k] = f_6 * kf_224[k]
                   + pa_x[k] * kg_175[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, pa_y, pb_x, pb_y, kf_170, kf_225, kf_226, \
                         kg_126, lf_225, lf_226, lf_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_5 * kf_170[k]
                   + pb_y[k] * lf_225[k];

        t_515[k] = pa_y[k] * kg_126[k];

        t_516[k] = f_5 * kf_225[k]
                   + pb_x[k] * lf_226[k];

        t_517[k] = f_5 * kf_226[k]
                   + pb_x[k] * lf_227[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, t_522, t_523, pa_x, pa_y, pb_x, kf_227, \
                         kg_127, kg_176, kg_177, kg_178, kg_179, \
                         lf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_5 * kf_227[k]
                   + pb_x[k] * lf_228[k];

        t_519[k] = pa_y[k] * kg_127[k];

        t_520[k] = pa_x[k] * kg_176[k];

        t_521[k] = pa_x[k] * kg_177[k];

        t_522[k] = pa_x[k] * kg_178[k];

        t_523[k] = pa_x[k] * kg_179[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, pa_x, pb_y, pb_z, kf_169, kf_229, \
                         kf_232, kg_180, kg_181, kg_183, lf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = pa_x[k] * kg_180[k];

        t_525[k] = f_8 * kf_229[k]
                   + pa_x[k] * kg_181[k];

        t_526[k] = pb_y[k] * lf_229[k];

        t_527[k] = f_7 * kf_169[k]
                   + pb_z[k] * lf_229[k];

        t_528[k] = f_6 * kf_232[k]
                   + pa_x[k] * kg_183[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, pa_x, pb_x, pb_y, kf_233, kf_234, \
                         kf_235, kg_184, lf_230, lf_231, lf_232, \
                         lf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = pb_y[k] * lf_230[k];

        t_530[k] = f_6 * kf_233[k]
                   + pa_x[k] * kg_184[k];

        t_531[k] = f_5 * kf_234[k]
                   + pb_x[k] * lf_232[k];

        t_532[k] = f_5 * kf_235[k]
                   + pb_x[k] * lf_233[k];

        t_533[k] = pb_y[k] * lf_231[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, t_539, pa_x, pb_x, pb_y, kf_237, \
                         kg_185, kg_186, kg_187, kg_188, lf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_5 * kf_237[k]
                   + pb_x[k] * lf_234[k];

        t_535[k] = pa_x[k] * kg_185[k];

        t_536[k] = pa_x[k] * kg_186[k];

        t_537[k] = pa_x[k] * kg_187[k];

        t_538[k] = pb_y[k] * lf_234[k];

        t_539[k] = pa_x[k] * kg_188[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, pb_x, pb_y, pb_z, kf_175, ld0_33, \
                         ld0_34, ld1_33, ld1_34, lf_235, lf_236, \
                         lf_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_235[k];

        t_541[k] = f_0 * kf_175[k]
                   + pb_y[k] * lf_235[k];

        t_542[k] = pb_z[k] * lf_235[k];

        t_543[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_237[k];

        t_544[k] = pb_z[k] * lf_236[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, pb_x, ld0_35, ld1_35, lf_238, \
                         lf_239, lf_240, lf_241, lf_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_238[k];

        t_546[k] = pb_x[k] * lf_239[k];

        t_547[k] = pb_x[k] * lf_240[k];

        t_548[k] = pb_x[k] * lf_241[k];

        t_549[k] = pb_x[k] * lf_242[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, pb_y, pb_z, kf_179, kf_182, \
                         ld0_34, ld0_35, ld1_34, ld1_35, lf_239, lf_240, \
                         lf_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_0 * kf_179[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_239[k];

        t_551[k] = pb_z[k] * lf_239[k];

        t_552[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_240[k];

        t_553[k] = f_0 * kf_182[k]
                   + pb_y[k] * lf_242[k];

        t_554[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_242[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, pa_z, pb_y, pb_z, kf_175, kf_184, \
                         kg_129, kg_130, kg_131, lf_243, lf_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_z[k] * kg_129[k];

        t_556[k] = pa_z[k] * kg_130[k];

        t_557[k] = f_5 * kf_175[k]
                   + pb_z[k] * lf_243[k];

        t_558[k] = pa_z[k] * kg_131[k];

        t_559[k] = f_7 * kf_184[k]
                   + pb_y[k] * lf_244[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, pa_z, pb_x, kf_176, kg_132, \
                         kg_133, lf_245, lf_246, lf_247, lf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_6 * kf_176[k]
                   + pa_z[k] * kg_132[k];

        t_561[k] = pb_x[k] * lf_245[k];

        t_562[k] = pb_x[k] * lf_246[k];

        t_563[k] = pb_x[k] * lf_247[k];

        t_564[k] = pb_x[k] * lf_248[k];

        t_565[k] = pa_z[k] * kg_133[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_z, pb_y, pb_z, kf_179, kf_180, kf_182, \
                         kf_189, kg_134, kg_136, lf_245, lf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_5 * kf_179[k]
                   + pb_z[k] * lf_245[k];

        t_567[k] = f_6 * kf_180[k]
                   + pa_z[k] * kg_134[k];

        t_568[k] = f_7 * kf_189[k]
                   + pb_y[k] * lf_248[k];

        t_569[k] = f_8 * kf_182[k]
                   + pa_z[k] * kg_136[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pb_x, pb_y, pb_z, kf_183, kf_190, ld0_36, \
                         ld0_37, ld1_36, ld1_37, lf_249, lf_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_249[k];

        t_571[k] = f_11 * kf_190[k]
                   + pb_y[k] * lf_249[k];

        t_572[k] = f_6 * kf_183[k]
                   + pb_z[k] * lf_249[k];

        t_573[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_251[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, pb_x, pb_y, kf_191, ld0_38, \
                         ld1_38, lf_250, lf_252, lf_253, lf_254, \
                         lf_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_11 * kf_191[k]
                   + pb_y[k] * lf_250[k];

        t_575[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_252[k];

        t_576[k] = pb_x[k] * lf_253[k];

        t_577[k] = pb_x[k] * lf_254[k];

        t_578[k] = pb_x[k] * lf_255[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pa_z, pb_x, pb_z, ig0_41, ig1_41, kf_186, \
                         kg_138, lf_253, lf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = pb_x[k] * lf_256[k];

        t_580[k] = f_9 * ig0_41[k]
                   - f_10 * ig1_41[k]
                   + pa_z[k] * kg_138[k];

        t_581[k] = f_6 * kf_186[k]
                   + pb_z[k] * lf_253[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pa_y, pb_y, ig0_45, ig1_45, kf_196, kf_197, \
                         kg_150, ld0_38, ld1_38, lf_255, lf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_11 * kf_196[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_255[k];

        t_583[k] = f_11 * kf_197[k]
                   + pb_y[k] * lf_256[k];

        t_584[k] = f_12 * ig0_45[k]
                   - f_13 * ig1_45[k]
                   + pa_y[k] * kg_150[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pb_x, pb_y, pb_z, kf_190, kf_198, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_257, lf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_257[k];

        t_586[k] = f_17 * kf_198[k]
                   + pb_y[k] * lf_257[k];

        t_587[k] = f_16 * kf_190[k]
                   + pb_z[k] * lf_257[k];

        t_588[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_259[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, pb_x, pb_y, kf_199, ld0_41, \
                         ld1_41, lf_258, lf_260, lf_261, lf_262, \
                         lf_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_17 * kf_199[k]
                   + pb_y[k] * lf_258[k];

        t_590[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_260[k];

        t_591[k] = pb_x[k] * lf_261[k];

        t_592[k] = pb_x[k] * lf_262[k];

        t_593[k] = pb_x[k] * lf_263[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pa_z, pb_x, pb_z, ig0_42, ig1_42, kf_194, \
                         kg_146, lf_261, lf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = pb_x[k] * lf_264[k];

        t_595[k] = f_14 * ig0_42[k]
                   - f_15 * ig1_42[k]
                   + pa_z[k] * kg_146[k];

        t_596[k] = f_16 * kf_194[k]
                   + pb_z[k] * lf_261[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pa_y, pb_y, ig0_48, ig1_48, kf_204, kf_205, \
                         kg_158, ld0_41, ld1_41, lf_263, lf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_17 * kf_204[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_263[k];

        t_598[k] = f_17 * kf_205[k]
                   + pb_y[k] * lf_264[k];

        t_599[k] = f_18 * ig0_48[k]
                   - f_19 * ig1_48[k]
                   + pa_y[k] * kg_158[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pb_x, pb_y, pb_z, kf_198, kf_206, ld0_42, \
                         ld0_43, ld1_42, ld1_43, lf_265, lf_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_265[k];

        t_601[k] = f_8 * kf_206[k]
                   + pb_y[k] * lf_265[k];

        t_602[k] = f_8 * kf_198[k]
                   + pb_z[k] * lf_265[k];

        t_603[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_267[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pb_x, pb_y, kf_207, ld0_44, \
                         ld1_44, lf_266, lf_268, lf_269, lf_270, \
                         lf_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_8 * kf_207[k]
                   + pb_y[k] * lf_266[k];

        t_605[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_268[k];

        t_606[k] = pb_x[k] * lf_269[k];

        t_607[k] = pb_x[k] * lf_270[k];

        t_608[k] = pb_x[k] * lf_271[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pa_z, pb_x, pb_z, ig0_43, ig1_43, kf_202, \
                         kg_154, lf_269, lf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pb_x[k] * lf_272[k];

        t_610[k] = f_20 * ig0_43[k]
                   - f_21 * ig1_43[k]
                   + pa_z[k] * kg_154[k];

        t_611[k] = f_8 * kf_202[k]
                   + pb_z[k] * lf_269[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_y, pb_y, ig0_51, ig1_51, kf_212, kf_213, \
                         kg_166, ld0_44, ld1_44, lf_271, lf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_8 * kf_212[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_271[k];

        t_613[k] = f_8 * kf_213[k]
                   + pb_y[k] * lf_272[k];

        t_614[k] = f_20 * ig0_51[k]
                   - f_21 * ig1_51[k]
                   + pa_y[k] * kg_166[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pb_x, pb_y, pb_z, kf_206, kf_214, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_273, lf_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_273[k];

        t_616[k] = f_16 * kf_214[k]
                   + pb_y[k] * lf_273[k];

        t_617[k] = f_17 * kf_206[k]
                   + pb_z[k] * lf_273[k];

        t_618[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_275[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, pb_x, pb_y, kf_215, ld0_47, \
                         ld1_47, lf_274, lf_276, lf_277, lf_278, \
                         lf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_16 * kf_215[k]
                   + pb_y[k] * lf_274[k];

        t_620[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_276[k];

        t_621[k] = pb_x[k] * lf_277[k];

        t_622[k] = pb_x[k] * lf_278[k];

        t_623[k] = pb_x[k] * lf_279[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pb_x, pb_z, ig0_46, ig1_46, kf_210, \
                         kg_162, lf_277, lf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pb_x[k] * lf_280[k];

        t_625[k] = f_18 * ig0_46[k]
                   - f_19 * ig1_46[k]
                   + pa_z[k] * kg_162[k];

        t_626[k] = f_17 * kf_210[k]
                   + pb_z[k] * lf_277[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_y, pb_y, ig0_52, ig1_52, kf_220, kf_221, \
                         kg_174, ld0_47, ld1_47, lf_279, lf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_16 * kf_220[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_279[k];

        t_628[k] = f_16 * kf_221[k]
                   + pb_y[k] * lf_280[k];

        t_629[k] = f_14 * ig0_52[k]
                   - f_15 * ig1_52[k]
                   + pa_y[k] * kg_174[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pb_x, pb_y, pb_z, kf_214, kf_222, ld0_48, \
                         ld0_49, ld1_48, ld1_49, lf_281, lf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_281[k];

        t_631[k] = f_6 * kf_222[k]
                   + pb_y[k] * lf_281[k];

        t_632[k] = f_11 * kf_214[k]
                   + pb_z[k] * lf_281[k];

        t_633[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_283[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pb_x, pb_y, kf_223, ld0_50, \
                         ld1_50, lf_282, lf_284, lf_285, lf_286, \
                         lf_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_6 * kf_223[k]
                   + pb_y[k] * lf_282[k];

        t_635[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_284[k];

        t_636[k] = pb_x[k] * lf_285[k];

        t_637[k] = pb_x[k] * lf_286[k];

        t_638[k] = pb_x[k] * lf_287[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pa_z, pb_x, pb_z, ig0_49, ig1_49, kf_218, \
                         kg_170, lf_285, lf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pb_x[k] * lf_288[k];

        t_640[k] = f_12 * ig0_49[k]
                   - f_13 * ig1_49[k]
                   + pa_z[k] * kg_170[k];

        t_641[k] = f_11 * kf_218[k]
                   + pb_z[k] * lf_285[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_y, pb_y, ig0_53, ig1_53, kf_227, \
                         kf_228, kg_180, kg_181, ld0_50, ld1_50, lf_287, \
                         lf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_6 * kf_227[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_287[k];

        t_643[k] = f_6 * kf_228[k]
                   + pb_y[k] * lf_288[k];

        t_644[k] = f_9 * ig0_53[k]
                   - f_10 * ig1_53[k]
                   + pa_y[k] * kg_180[k];

        t_645[k] = pa_y[k] * kg_181[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, t_650, pa_y, pb_y, kf_229, kf_230, \
                         kf_231, kg_182, kg_183, kg_184, lf_289, \
                         lf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_5 * kf_229[k]
                   + pb_y[k] * lf_289[k];

        t_647[k] = pa_y[k] * kg_182[k];

        t_648[k] = f_6 * kf_230[k]
                   + pa_y[k] * kg_183[k];

        t_649[k] = f_5 * kf_231[k]
                   + pb_y[k] * lf_290[k];

        t_650[k] = pa_y[k] * kg_184[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, t_655, t_656, pa_y, pb_x, pb_z, kf_225, \
                         kf_234, kg_185, lf_291, lf_292, lf_293, \
                         lf_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = pb_x[k] * lf_291[k];

        t_652[k] = pb_x[k] * lf_292[k];

        t_653[k] = pb_x[k] * lf_293[k];

        t_654[k] = pb_x[k] * lf_294[k];

        t_655[k] = f_8 * kf_234[k]
                   + pa_y[k] * kg_185[k];

        t_656[k] = f_7 * kf_225[k]
                   + pb_z[k] * lf_291[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, pa_y, pb_x, pb_y, kf_236, kf_237, \
                         kg_187, kg_188, ld0_51, ld1_51, lf_294, \
                         lf_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_6 * kf_236[k]
                   + pa_y[k] * kg_187[k];

        t_658[k] = f_5 * kf_237[k]
                   + pb_y[k] * lf_294[k];

        t_659[k] = pa_y[k] * kg_188[k];

        t_660[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_295[k];

        t_661[k] = pb_y[k] * lf_295[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pb_x, pb_y, pb_z, kf_229, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_295, lf_296, lf_297, \
                         lf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_0 * kf_229[k]
                   + pb_z[k] * lf_295[k];

        t_663[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_297[k];

        t_664[k] = pb_y[k] * lf_296[k];

        t_665[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_298[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, t_670, t_671, pb_x, pb_y, pb_z, kf_234, \
                         ld0_52, ld1_52, lf_299, lf_300, lf_301, \
                         lf_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = pb_x[k] * lf_299[k];

        t_667[k] = pb_x[k] * lf_300[k];

        t_668[k] = pb_x[k] * lf_301[k];

        t_669[k] = pb_x[k] * lf_302[k];

        t_670[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_299[k];

        t_671[k] = f_0 * kf_234[k]
                   + pb_z[k] * lf_299[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, pb_y, pb_z, kf_237, ld0_53, ld1_53, lf_301, \
                         lf_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_301[k];

        t_673[k] = pb_y[k] * lf_302[k];

        t_674[k] = f_0 * kf_237[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_302[k];
    }
}

auto
compute_prim_lg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_1 = buffer.data(ig0 + 1);
    const auto *ig0_2 = buffer.data(ig0 + 2);
    const auto *ig0_3 = buffer.data(ig0 + 3);
    const auto *ig0_4 = buffer.data(ig0 + 4);
    const auto *ig0_6 = buffer.data(ig0 + 6);
    const auto *ig0_7 = buffer.data(ig0 + 7);
    const auto *ig0_8 = buffer.data(ig0 + 8);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_12 = buffer.data(ig0 + 12);
    const auto *ig0_14 = buffer.data(ig0 + 14);
    const auto *ig0_15 = buffer.data(ig0 + 15);
    const auto *ig0_16 = buffer.data(ig0 + 16);
    const auto *ig0_17 = buffer.data(ig0 + 17);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_19 = buffer.data(ig0 + 19);
    const auto *ig0_21 = buffer.data(ig0 + 21);
    const auto *ig0_22 = buffer.data(ig0 + 22);
    const auto *ig0_23 = buffer.data(ig0 + 23);
    const auto *ig0_25 = buffer.data(ig0 + 25);
    const auto *ig0_26 = buffer.data(ig0 + 26);
    const auto *ig0_27 = buffer.data(ig0 + 27);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_31 = buffer.data(ig0 + 31);
    const auto *ig0_32 = buffer.data(ig0 + 32);
    const auto *ig0_33 = buffer.data(ig0 + 33);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_35 = buffer.data(ig0 + 35);
    const auto *ig0_36 = buffer.data(ig0 + 36);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_39 = buffer.data(ig0 + 39);
    const auto *ig0_40 = buffer.data(ig0 + 40);
    const auto *ig0_41 = buffer.data(ig0 + 41);
    const auto *ig0_42 = buffer.data(ig0 + 42);
    const auto *ig0_43 = buffer.data(ig0 + 43);
    const auto *ig0_44 = buffer.data(ig0 + 44);
    const auto *ig0_45 = buffer.data(ig0 + 45);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);
    const auto *ig0_48 = buffer.data(ig0 + 48);
    const auto *ig0_49 = buffer.data(ig0 + 49);
    const auto *ig0_50 = buffer.data(ig0 + 50);
    const auto *ig0_52 = buffer.data(ig0 + 52);
    const auto *ig0_53 = buffer.data(ig0 + 53);
    const auto *ig0_54 = buffer.data(ig0 + 54);
    const auto *ig0_56 = buffer.data(ig0 + 56);
    const auto *ig0_57 = buffer.data(ig0 + 57);
    const auto *ig0_58 = buffer.data(ig0 + 58);
    const auto *ig0_60 = buffer.data(ig0 + 60);
    const auto *ig0_61 = buffer.data(ig0 + 61);
    const auto *ig0_62 = buffer.data(ig0 + 62);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_9 = buffer.data(ig1 + 9);
    const auto *ig1_12 = buffer.data(ig1 + 12);
    const auto *ig1_17 = buffer.data(ig1 + 17);
    const auto *ig1_18 = buffer.data(ig1 + 18);
    const auto *ig1_21 = buffer.data(ig1 + 21);
    const auto *ig1_24 = buffer.data(ig1 + 24);
    const auto *ig1_27 = buffer.data(ig1 + 27);
    const auto *ig1_31 = buffer.data(ig1 + 31);
    const auto *ig1_32 = buffer.data(ig1 + 32);
    const auto *ig1_33 = buffer.data(ig1 + 33);
    const auto *ig1_36 = buffer.data(ig1 + 36);
    const auto *ig1_39 = buffer.data(ig1 + 39);
    const auto *ig1_40 = buffer.data(ig1 + 40);
    const auto *ig1_41 = buffer.data(ig1 + 41);
    const auto *ig1_42 = buffer.data(ig1 + 42);
    const auto *ig1_45 = buffer.data(ig1 + 45);
    const auto *ig1_49 = buffer.data(ig1 + 49);
    const auto *ig1_50 = buffer.data(ig1 + 50);
    const auto *ig1_51 = buffer.data(ig1 + 51);
    const auto *ig1_54 = buffer.data(ig1 + 54);
    const auto *ig1_57 = buffer.data(ig1 + 57);
    const auto *ig1_58 = buffer.data(ig1 + 58);
    const auto *ig1_59 = buffer.data(ig1 + 59);
    const auto *ig1_60 = buffer.data(ig1 + 60);
    const auto *ig1_61 = buffer.data(ig1 + 61);
    const auto *ig1_62 = buffer.data(ig1 + 62);
    const auto *ig1_63 = buffer.data(ig1 + 63);
    const auto *ig1_64 = buffer.data(ig1 + 64);
    const auto *ig1_65 = buffer.data(ig1 + 65);
    const auto *ig1_66 = buffer.data(ig1 + 66);
    const auto *ig1_69 = buffer.data(ig1 + 69);
    const auto *ig1_73 = buffer.data(ig1 + 73);
    const auto *ig1_77 = buffer.data(ig1 + 77);
    const auto *ig1_78 = buffer.data(ig1 + 78);
    const auto *ig1_79 = buffer.data(ig1 + 79);
    const auto *ig1_80 = buffer.data(ig1 + 80);
    const auto *ig1_81 = buffer.data(ig1 + 81);
    const auto *ig1_82 = buffer.data(ig1 + 82);
    const auto *ig1_83 = buffer.data(ig1 + 83);
    const auto *ig1_88 = buffer.data(ig1 + 88);
    const auto *ig1_94 = buffer.data(ig1 + 94);
    const auto *ig1_100 = buffer.data(ig1 + 100);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_110 = buffer.data(ig1 + 110);
    const auto *ig1_112 = buffer.data(ig1 + 112);
    const auto *ig1_116 = buffer.data(ig1 + 116);
    const auto *ig1_118 = buffer.data(ig1 + 118);
    const auto *ig1_120 = buffer.data(ig1 + 120);
    const auto *ig1_124 = buffer.data(ig1 + 124);
    const auto *ig1_126 = buffer.data(ig1 + 126);
    const auto *ig1_128 = buffer.data(ig1 + 128);
    const auto *ig1_134 = buffer.data(ig1 + 134);
    const auto *ig1_145 = buffer.data(ig1 + 145);

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
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
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
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_3 = buffer.data(kg + 3);
    const auto *kg_4 = buffer.data(kg + 4);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_6 = buffer.data(kg + 6);
    const auto *kg_8 = buffer.data(kg + 8);
    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_19 = buffer.data(kg + 19);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_36 = buffer.data(kg + 36);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_43 = buffer.data(kg + 43);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_45 = buffer.data(kg + 45);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_54 = buffer.data(kg + 54);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_79 = buffer.data(kg + 79);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_93 = buffer.data(kg + 93);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_97 = buffer.data(kg + 97);
    const auto *kg_98 = buffer.data(kg + 98);
    const auto *kg_99 = buffer.data(kg + 99);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_146 = buffer.data(kg + 146);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_148 = buffer.data(kg + 148);
    const auto *kg_149 = buffer.data(kg + 149);
    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_152 = buffer.data(kg + 152);
    const auto *kg_153 = buffer.data(kg + 153);
    const auto *kg_154 = buffer.data(kg + 154);
    const auto *kg_157 = buffer.data(kg + 157);
    const auto *kg_158 = buffer.data(kg + 158);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_160 = buffer.data(kg + 160);
    const auto *kg_161 = buffer.data(kg + 161);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_163 = buffer.data(kg + 163);
    const auto *kg_164 = buffer.data(kg + 164);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_169 = buffer.data(kg + 169);
    const auto *kg_170 = buffer.data(kg + 170);
    const auto *kg_171 = buffer.data(kg + 171);
    const auto *kg_172 = buffer.data(kg + 172);
    const auto *kg_173 = buffer.data(kg + 173);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_178 = buffer.data(kg + 178);
    const auto *kg_179 = buffer.data(kg + 179);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_181 = buffer.data(kg + 181);
    const auto *kg_182 = buffer.data(kg + 182);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_184 = buffer.data(kg + 184);
    const auto *kg_187 = buffer.data(kg + 187);
    const auto *kg_188 = buffer.data(kg + 188);
    const auto *kg_189 = buffer.data(kg + 189);
    const auto *kg_190 = buffer.data(kg + 190);
    const auto *kg_191 = buffer.data(kg + 191);
    const auto *kg_192 = buffer.data(kg + 192);
    const auto *kg_193 = buffer.data(kg + 193);
    const auto *kg_194 = buffer.data(kg + 194);
    const auto *kg_195 = buffer.data(kg + 195);
    const auto *kg_196 = buffer.data(kg + 196);
    const auto *kg_197 = buffer.data(kg + 197);
    const auto *kg_198 = buffer.data(kg + 198);
    const auto *kg_199 = buffer.data(kg + 199);
    const auto *kg_200 = buffer.data(kg + 200);
    const auto *kg_201 = buffer.data(kg + 201);
    const auto *kg_204 = buffer.data(kg + 204);
    const auto *kg_205 = buffer.data(kg + 205);
    const auto *kg_206 = buffer.data(kg + 206);
    const auto *kg_208 = buffer.data(kg + 208);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_111 = buffer.data(lf + 111);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_113 = buffer.data(lf + 113);
    const auto *lf_114 = buffer.data(lf + 114);
    const auto *lf_115 = buffer.data(lf + 115);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_121 = buffer.data(lf + 121);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_123 = buffer.data(lf + 123);
    const auto *lf_124 = buffer.data(lf + 124);
    const auto *lf_125 = buffer.data(lf + 125);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_131 = buffer.data(lf + 131);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_133 = buffer.data(lf + 133);
    const auto *lf_134 = buffer.data(lf + 134);
    const auto *lf_135 = buffer.data(lf + 135);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_143 = buffer.data(lf + 143);
    const auto *lf_144 = buffer.data(lf + 144);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_154 = buffer.data(lf + 154);
    const auto *lf_155 = buffer.data(lf + 155);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_161 = buffer.data(lf + 161);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_163 = buffer.data(lf + 163);
    const auto *lf_164 = buffer.data(lf + 164);
    const auto *lf_165 = buffer.data(lf + 165);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, kf_3, kf_6, ld0_1, ld0_2, ld1_1, \
                         ld1_2, lf_3, lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * kf_3[k]
                 + pb_x[k] * lf_3[k];

        t_6[k] = f_0 * kf_6[k]
                 + pb_x[k] * lf_5[k];

        t_7[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_8[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_9[k] = pb_y[k] * lf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, kf_0, kf_1, kg_0, kg_3, \
                         ld0_2, ld1_2, lf_5, lf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * ld0_2[k]
                  - f_2 * ld1_2[k]
                  + pb_z[k] * lf_5[k];

        t_11[k] = pa_y[k] * kg_0[k];

        t_12[k] = f_5 * kf_0[k]
                  + pb_y[k] * lf_6[k];

        t_13[k] = f_6 * kf_1[k]
                  + pa_y[k] * kg_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pb_x, kf_3, kf_5, kf_8, kg_4, kg_5, \
                         kg_6, lf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * kg_4[k];

        t_15[k] = f_7 * kf_8[k]
                  + pb_x[k] * lf_7[k];

        t_16[k] = f_8 * kf_3[k]
                  + pa_y[k] * kg_5[k];

        t_17[k] = f_6 * kf_5[k]
                  + pa_y[k] * kg_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, kf_0, kf_6, \
                         kg_0, kg_3, kg_8, lf_8, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * kf_6[k]
                  + pb_y[k] * lf_8[k];

        t_19[k] = pa_y[k] * kg_8[k];

        t_20[k] = pa_z[k] * kg_0[k];

        t_21[k] = f_5 * kf_0[k]
                  + pb_z[k] * lf_9[k];

        t_22[k] = pa_z[k] * kg_3[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_z, kf_2, kf_3, kf_13, kg_4, \
                         kg_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * kf_2[k]
                  + pa_z[k] * kg_4[k];

        t_24[k] = f_7 * kf_13[k]
                  + pb_x[k] * lf_11[k];

        t_25[k] = pa_z[k] * kg_5[k];

        t_26[k] = f_5 * kf_3[k]
                  + pb_z[k] * lf_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, ig0_0, ig1_0, kf_4, kf_6, \
                         kf_7, kg_6, kg_8, kg_9, lf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * kf_4[k]
                  + pa_z[k] * kg_6[k];

        t_28[k] = f_8 * kf_6[k]
                  + pa_z[k] * kg_8[k];

        t_29[k] = f_9 * ig0_0[k]
                  - f_10 * ig1_0[k]
                  + pa_y[k] * kg_9[k];

        t_30[k] = f_6 * kf_7[k]
                  + pb_y[k] * lf_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, pb_z, kf_16, kf_17, ld0_3, ld0_4, \
                         ld1_3, ld1_4, lf_12, lf_13, lf_14, lf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_z[k] * lf_12[k];

        t_32[k] = f_11 * kf_16[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_14[k];

        t_33[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_13[k];

        t_34[k] = f_11 * kf_17[k]
                  + pb_x[k] * lf_15[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, pb_z, ig0_6, ig1_21, kf_9, kg_22, \
                         ld0_4, ld1_4, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_12 * ig0_6[k]
                  - f_13 * ig1_21[k]
                  + pa_x[k] * kg_22[k];

        t_36[k] = pb_z[k] * lf_15[k];

        t_37[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_16[k];

        t_38[k] = f_6 * kf_9[k]
                  + pb_y[k] * lf_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pa_z, pb_z, kg_10, kg_11, kg_13, \
                         kg_14, ld0_5, ld1_5, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_17[k];

        t_40[k] = pa_y[k] * kg_13[k];

        t_41[k] = pa_z[k] * kg_10[k];

        t_42[k] = pa_y[k] * kg_14[k];

        t_43[k] = pa_z[k] * kg_11[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, kf_8, kf_12, kf_13, kg_15, \
                         kg_16, lf_18, lf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * kf_8[k]
                  + pb_z[k] * lf_18[k];

        t_45[k] = f_6 * kf_12[k]
                  + pa_y[k] * kg_15[k];

        t_46[k] = f_5 * kf_13[k]
                  + pb_y[k] * lf_19[k];

        t_47[k] = pa_y[k] * kg_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, ig0_0, ig1_0, kf_10, kg_12, \
                         ld0_6, ld1_6, lf_20, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_9 * ig0_0[k]
                  - f_10 * ig1_0[k]
                  + pa_z[k] * kg_12[k];

        t_49[k] = pb_y[k] * lf_20[k];

        t_50[k] = f_6 * kf_10[k]
                  + pb_z[k] * lf_20[k];

        t_51[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_21[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, kf_24, kf_27, ld0_7, ld0_8, ld1_7, \
                         ld1_8, lf_22, lf_23, lf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_11 * kf_24[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_22[k];

        t_53[k] = f_11 * kf_27[k]
                  + pb_x[k] * lf_25[k];

        t_54[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_23[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pb_y, pb_z, ig0_10, ig1_31, kf_11, \
                         kg_35, ld0_8, ld1_8, lf_23, lf_24, lf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_6 * kf_11[k]
                  + pb_z[k] * lf_23[k];

        t_56[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_24[k];

        t_57[k] = pb_y[k] * lf_25[k];

        t_58[k] = f_12 * ig0_10[k]
                  - f_13 * ig1_31[k]
                  + pa_x[k] * kg_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_y, pb_y, pb_z, ig0_1, ig1_9, kf_14, kg_17, \
                         lf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_14 * ig0_1[k]
                  - f_15 * ig1_9[k]
                  + pa_y[k] * kg_17[k];

        t_60[k] = f_16 * kf_14[k]
                  + pb_y[k] * lf_26[k];

        t_61[k] = pb_z[k] * lf_26[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pb_x, pb_z, kf_30, kf_31, ld0_9, ld0_10, ld1_9, \
                         ld1_10, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_17 * kf_30[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_28[k];

        t_63[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_27[k];

        t_64[k] = f_17 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_y, pb_z, ig0_14, ig1_36, kf_19, \
                         kg_41, ld0_10, ld1_10, lf_29, lf_30, lf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_18 * ig0_14[k]
                  - f_19 * ig1_36[k]
                  + pa_x[k] * kg_41[k];

        t_66[k] = pb_z[k] * lf_29[k];

        t_67[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_30[k];

        t_68[k] = f_16 * kf_19[k]
                  + pb_y[k] * lf_31[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_z, pb_z, kf_14, kf_15, kg_17, kg_19, \
                         kg_20, ld0_11, ld1_11, lf_31, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_31[k];

        t_70[k] = pa_z[k] * kg_17[k];

        t_71[k] = f_5 * kf_14[k]
                  + pb_z[k] * lf_32[k];

        t_72[k] = pa_z[k] * kg_19[k];

        t_73[k] = f_6 * kf_15[k]
                  + pa_z[k] * kg_20[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pb_y, pb_z, kf_17, kf_18, kf_21, kg_22, \
                         kg_24, lf_33, lf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_z[k] * kg_22[k];

        t_75[k] = f_5 * kf_17[k]
                  + pb_z[k] * lf_33[k];

        t_76[k] = f_6 * kf_18[k]
                  + pa_z[k] * kg_24[k];

        t_77[k] = f_6 * kf_21[k]
                  + pb_y[k] * lf_34[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, kf_19, kf_23, kg_25, kg_26, \
                         kg_28, kg_29, kg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * kf_19[k]
                  + pa_z[k] * kg_25[k];

        t_79[k] = pa_y[k] * kg_26[k];

        t_80[k] = pa_y[k] * kg_28[k];

        t_81[k] = f_6 * kf_23[k]
                  + pa_y[k] * kg_29[k];

        t_82[k] = pa_y[k] * kg_30[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_y, pb_z, kf_20, kf_25, kf_26, kf_27, \
                         kg_32, kg_33, lf_35, lf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_8 * kf_25[k]
                  + pa_y[k] * kg_32[k];

        t_84[k] = f_6 * kf_20[k]
                  + pb_z[k] * lf_35[k];

        t_85[k] = f_6 * kf_26[k]
                  + pa_y[k] * kg_33[k];

        t_86[k] = f_5 * kf_27[k]
                  + pb_y[k] * lf_36[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_y, pb_z, ig0_2, ig1_12, kf_22, \
                         kg_26, kg_35, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_y[k] * kg_35[k];

        t_88[k] = f_14 * ig0_2[k]
                  - f_15 * ig1_12[k]
                  + pa_z[k] * kg_26[k];

        t_89[k] = pb_y[k] * lf_37[k];

        t_90[k] = f_16 * kf_22[k]
                  + pb_z[k] * lf_37[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_y, kf_41, kf_44, ld0_12, ld0_14, ld1_12, \
                         ld1_14, lf_38, lf_39, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_38[k];

        t_92[k] = f_17 * kf_41[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_39[k];

        t_93[k] = f_17 * kf_44[k]
                  + pb_x[k] * lf_42[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_y, pb_z, kf_25, ld0_13, ld0_14, ld1_13, \
                         ld1_14, lf_40, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_40[k];

        t_95[k] = f_16 * kf_25[k]
                  + pb_z[k] * lf_40[k];

        t_96[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_41[k];

        t_97[k] = pb_y[k] * lf_42[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_x, pa_y, pb_y, pb_z, ig0_3, ig0_21, \
                         ig1_17, ig1_49, kf_28, kg_36, kg_57, lf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_18 * ig0_21[k]
                  - f_19 * ig1_49[k]
                  + pa_x[k] * kg_57[k];

        t_99[k] = f_20 * ig0_3[k]
                  - f_21 * ig1_17[k]
                  + pa_y[k] * kg_36[k];

        t_100[k] = f_8 * kf_28[k]
                   + pb_y[k] * lf_43[k];

        t_101[k] = pb_z[k] * lf_43[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_x, pb_z, kf_47, kf_48, ld0_15, ld0_16, \
                         ld1_15, ld1_16, lf_44, lf_45, lf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_8 * kf_47[k]
                   + f_3 * ld0_16[k]
                   - f_4 * ld1_16[k]
                   + pb_x[k] * lf_45[k];

        t_103[k] = f_3 * ld0_15[k]
                   - f_4 * ld1_15[k]
                   + pb_z[k] * lf_44[k];

        t_104[k] = f_8 * kf_48[k]
                   + pb_x[k] * lf_46[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_y, pb_z, ig0_25, ig1_54, kf_33, \
                         kg_63, ld0_16, ld1_16, lf_46, lf_47, lf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_20 * ig0_25[k]
                   - f_21 * ig1_54[k]
                   + pa_x[k] * kg_63[k];

        t_106[k] = pb_z[k] * lf_46[k];

        t_107[k] = f_3 * ld0_16[k]
                   - f_4 * ld1_16[k]
                   + pb_z[k] * lf_47[k];

        t_108[k] = f_8 * kf_33[k]
                   + pb_y[k] * lf_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pa_z, pb_z, kf_28, kf_29, kg_36, \
                         kg_38, kg_39, ld0_17, ld1_17, lf_48, lf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_1 * ld0_17[k]
                   - f_2 * ld1_17[k]
                   + pb_z[k] * lf_48[k];

        t_110[k] = pa_z[k] * kg_36[k];

        t_111[k] = f_5 * kf_28[k]
                   + pb_z[k] * lf_49[k];

        t_112[k] = pa_z[k] * kg_38[k];

        t_113[k] = f_6 * kf_29[k]
                   + pa_z[k] * kg_39[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pb_y, pb_z, kf_31, kf_32, kf_36, \
                         kg_41, kg_43, lf_50, lf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_z[k] * kg_41[k];

        t_115[k] = f_5 * kf_31[k]
                   + pb_z[k] * lf_50[k];

        t_116[k] = f_6 * kf_32[k]
                   + pa_z[k] * kg_43[k];

        t_117[k] = f_16 * kf_36[k]
                   + pb_y[k] * lf_51[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_y, pa_z, pb_z, ig0_7, ig1_24, kf_33, kf_34, \
                         kg_44, kg_46, lf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * kf_33[k]
                   + pa_z[k] * kg_44[k];

        t_119[k] = f_9 * ig0_7[k]
                   - f_10 * ig1_24[k]
                   + pa_y[k] * kg_46[k];

        t_120[k] = f_6 * kf_34[k]
                   + pb_z[k] * lf_52[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_x, pa_y, pa_z, ig0_4, ig0_8, ig0_30, ig1_18, \
                         ig1_27, ig1_61, kg_45, kg_47, kg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_9 * ig0_4[k]
                   - f_10 * ig1_18[k]
                   + pa_z[k] * kg_45[k];

        t_122[k] = f_9 * ig0_8[k]
                   - f_10 * ig1_27[k]
                   + pa_y[k] * kg_47[k];

        t_123[k] = f_20 * ig0_30[k]
                   - f_21 * ig1_61[k]
                   + pa_x[k] * kg_71[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pb_y, pb_z, ig0_31, ig1_62, kf_35, kf_38, \
                         kg_72, lf_53, lf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_6 * kf_35[k]
                   + pb_z[k] * lf_53[k];

        t_125[k] = f_20 * ig0_31[k]
                   - f_21 * ig1_62[k]
                   + pa_x[k] * kg_72[k];

        t_126[k] = f_6 * kf_38[k]
                   + pb_y[k] * lf_54[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_x, pa_y, ig0_32, ig1_63, kf_40, \
                         kg_48, kg_50, kg_51, kg_52, kg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_20 * ig0_32[k]
                   - f_21 * ig1_63[k]
                   + pa_x[k] * kg_73[k];

        t_128[k] = pa_y[k] * kg_48[k];

        t_129[k] = pa_y[k] * kg_50[k];

        t_130[k] = f_6 * kf_40[k]
                   + pa_y[k] * kg_51[k];

        t_131[k] = pa_y[k] * kg_52[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_y, pb_y, pb_z, kf_37, kf_42, kf_43, \
                         kf_44, kg_54, kg_55, lf_55, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_8 * kf_42[k]
                   + pa_y[k] * kg_54[k];

        t_133[k] = f_16 * kf_37[k]
                   + pb_z[k] * lf_55[k];

        t_134[k] = f_6 * kf_43[k]
                   + pa_y[k] * kg_55[k];

        t_135[k] = f_5 * kf_44[k]
                   + pb_y[k] * lf_56[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_y, pa_z, pb_y, pb_z, ig0_7, ig1_24, \
                         kf_39, kg_48, kg_57, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_y[k] * kg_57[k];

        t_137[k] = f_20 * ig0_7[k]
                   - f_21 * ig1_24[k]
                   + pa_z[k] * kg_48[k];

        t_138[k] = pb_y[k] * lf_57[k];

        t_139[k] = f_8 * kf_39[k]
                   + pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, pb_y, kf_61, kf_64, ld0_18, ld0_20, \
                         ld1_18, ld1_20, lf_58, lf_59, lf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_3 * ld0_18[k]
                   - f_4 * ld1_18[k]
                   + pb_y[k] * lf_58[k];

        t_141[k] = f_8 * kf_61[k]
                   + f_3 * ld0_20[k]
                   - f_4 * ld1_20[k]
                   + pb_x[k] * lf_59[k];

        t_142[k] = f_8 * kf_64[k]
                   + pb_x[k] * lf_62[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pb_y, pb_z, kf_42, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_60, lf_61, lf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_1 * ld0_19[k]
                   - f_2 * ld1_19[k]
                   + pb_y[k] * lf_60[k];

        t_144[k] = f_8 * kf_42[k]
                   + pb_z[k] * lf_60[k];

        t_145[k] = f_3 * ld0_20[k]
                   - f_4 * ld1_20[k]
                   + pb_y[k] * lf_61[k];

        t_146[k] = pb_y[k] * lf_62[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_x, pa_y, pb_y, pb_z, ig0_11, ig0_38, \
                         ig1_32, ig1_73, kf_45, kg_58, kg_85, lf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_20 * ig0_38[k]
                   - f_21 * ig1_73[k]
                   + pa_x[k] * kg_85[k];

        t_148[k] = f_18 * ig0_11[k]
                   - f_19 * ig1_32[k]
                   + pa_y[k] * kg_58[k];

        t_149[k] = f_17 * kf_45[k]
                   + pb_y[k] * lf_63[k];

        t_150[k] = pb_z[k] * lf_63[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pb_x, pb_z, kf_67, kf_68, ld0_21, ld0_22, \
                         ld1_21, ld1_22, lf_64, lf_65, lf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_16 * kf_67[k]
                   + f_3 * ld0_22[k]
                   - f_4 * ld1_22[k]
                   + pb_x[k] * lf_65[k];

        t_152[k] = f_3 * ld0_21[k]
                   - f_4 * ld1_21[k]
                   + pb_z[k] * lf_64[k];

        t_153[k] = f_16 * kf_68[k]
                   + pb_x[k] * lf_66[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_y, pb_z, ig0_39, ig1_77, kf_50, \
                         kg_91, ld0_22, ld1_22, lf_66, lf_67, lf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_14 * ig0_39[k]
                   - f_15 * ig1_77[k]
                   + pa_x[k] * kg_91[k];

        t_155[k] = pb_z[k] * lf_66[k];

        t_156[k] = f_3 * ld0_22[k]
                   - f_4 * ld1_22[k]
                   + pb_z[k] * lf_67[k];

        t_157[k] = f_17 * kf_50[k]
                   + pb_y[k] * lf_68[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_z, kf_45, kf_46, kg_58, \
                         kg_60, kg_61, ld0_23, ld1_23, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_1 * ld0_23[k]
                   - f_2 * ld1_23[k]
                   + pb_z[k] * lf_68[k];

        t_159[k] = pa_z[k] * kg_58[k];

        t_160[k] = f_5 * kf_45[k]
                   + pb_z[k] * lf_69[k];

        t_161[k] = pa_z[k] * kg_60[k];

        t_162[k] = f_6 * kf_46[k]
                   + pa_z[k] * kg_61[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_z, pb_y, pb_z, kf_48, kf_49, kf_53, \
                         kg_63, kg_65, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * kg_63[k];

        t_164[k] = f_5 * kf_48[k]
                   + pb_z[k] * lf_70[k];

        t_165[k] = f_6 * kf_49[k]
                   + pa_z[k] * kg_65[k];

        t_166[k] = f_8 * kf_53[k]
                   + pb_y[k] * lf_71[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_y, pa_z, pb_z, ig0_16, ig1_40, kf_50, kf_51, \
                         kg_66, kg_68, lf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_8 * kf_50[k]
                   + pa_z[k] * kg_66[k];

        t_168[k] = f_14 * ig0_16[k]
                   - f_15 * ig1_40[k]
                   + pa_y[k] * kg_68[k];

        t_169[k] = f_6 * kf_51[k]
                   + pb_z[k] * lf_72[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pa_x, pa_y, pa_z, ig0_12, ig0_17, ig0_40, \
                         ig1_33, ig1_41, ig1_78, kg_67, kg_70, kg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_9 * ig0_12[k]
                   - f_10 * ig1_33[k]
                   + pa_z[k] * kg_67[k];

        t_171[k] = f_14 * ig0_17[k]
                   - f_15 * ig1_41[k]
                   + pa_y[k] * kg_70[k];

        t_172[k] = f_14 * ig0_40[k]
                   - f_15 * ig1_78[k]
                   + pa_x[k] * kg_99[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pa_x, pb_y, pb_z, ig0_41, ig1_79, kf_52, kf_56, \
                         kg_100, lf_73, lf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_6 * kf_52[k]
                   + pb_z[k] * lf_73[k];

        t_174[k] = f_14 * ig0_41[k]
                   - f_15 * ig1_79[k]
                   + pa_x[k] * kg_100[k];

        t_175[k] = f_16 * kf_56[k]
                   + pb_y[k] * lf_74[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, pa_y, pb_z, ig0_18, ig0_42, ig1_42, \
                         ig1_80, kf_54, kg_74, kg_101, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_14 * ig0_42[k]
                   - f_15 * ig1_80[k]
                   + pa_x[k] * kg_101[k];

        t_177[k] = f_9 * ig0_18[k]
                   - f_10 * ig1_42[k]
                   + pa_y[k] * kg_74[k];

        t_178[k] = f_16 * kf_54[k]
                   + pb_z[k] * lf_75[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pa_z, ig0_15, ig0_19, ig0_43, \
                         ig1_39, ig1_45, ig1_81, kg_69, kg_75, kg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_14 * ig0_15[k]
                   - f_15 * ig1_39[k]
                   + pa_z[k] * kg_69[k];

        t_180[k] = f_9 * ig0_19[k]
                   - f_10 * ig1_45[k]
                   + pa_y[k] * kg_75[k];

        t_181[k] = f_14 * ig0_43[k]
                   - f_15 * ig1_81[k]
                   + pa_x[k] * kg_105[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pa_x, pb_y, pb_z, ig0_44, ig1_82, kf_55, kf_58, \
                         kg_106, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_16 * kf_55[k]
                   + pb_z[k] * lf_76[k];

        t_183[k] = f_14 * ig0_44[k]
                   - f_15 * ig1_82[k]
                   + pa_x[k] * kg_106[k];

        t_184[k] = f_6 * kf_58[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, pa_x, pa_y, ig0_45, ig1_83, kf_60, \
                         kg_76, kg_78, kg_79, kg_80, kg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * ig0_45[k]
                   - f_15 * ig1_83[k]
                   + pa_x[k] * kg_107[k];

        t_186[k] = pa_y[k] * kg_76[k];

        t_187[k] = pa_y[k] * kg_78[k];

        t_188[k] = f_6 * kf_60[k]
                   + pa_y[k] * kg_79[k];

        t_189[k] = pa_y[k] * kg_80[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_y, pb_y, pb_z, kf_57, kf_62, kf_63, \
                         kf_64, kg_82, kg_83, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * kf_62[k]
                   + pa_y[k] * kg_82[k];

        t_191[k] = f_8 * kf_57[k]
                   + pb_z[k] * lf_78[k];

        t_192[k] = f_6 * kf_63[k]
                   + pa_y[k] * kg_83[k];

        t_193[k] = f_5 * kf_64[k]
                   + pb_y[k] * lf_79[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pa_y, pa_z, pb_y, pb_z, ig0_18, ig1_42, \
                         kf_59, kg_76, kg_85, lf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_y[k] * kg_85[k];

        t_195[k] = f_18 * ig0_18[k]
                   - f_19 * ig1_42[k]
                   + pa_z[k] * kg_76[k];

        t_196[k] = pb_y[k] * lf_80[k];

        t_197[k] = f_17 * kf_59[k]
                   + pb_z[k] * lf_80[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_x, pb_y, kf_84, kf_87, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_81, lf_82, lf_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * ld0_24[k]
                   - f_4 * ld1_24[k]
                   + pb_y[k] * lf_81[k];

        t_199[k] = f_16 * kf_84[k]
                   + f_3 * ld0_26[k]
                   - f_4 * ld1_26[k]
                   + pb_x[k] * lf_82[k];

        t_200[k] = f_16 * kf_87[k]
                   + pb_x[k] * lf_85[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pb_y, pb_z, kf_62, ld0_25, ld0_26, \
                         ld1_25, ld1_26, lf_83, lf_84, lf_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_1 * ld0_25[k]
                   - f_2 * ld1_25[k]
                   + pb_y[k] * lf_83[k];

        t_202[k] = f_17 * kf_62[k]
                   + pb_z[k] * lf_83[k];

        t_203[k] = f_3 * ld0_26[k]
                   - f_4 * ld1_26[k]
                   + pb_y[k] * lf_84[k];

        t_204[k] = pb_y[k] * lf_85[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pa_x, pa_y, pb_y, pb_z, ig0_22, ig0_46, \
                         ig1_50, ig1_88, kf_65, kg_86, kg_119, lf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_14 * ig0_46[k]
                   - f_15 * ig1_88[k]
                   + pa_x[k] * kg_119[k];

        t_206[k] = f_12 * ig0_22[k]
                   - f_13 * ig1_50[k]
                   + pa_y[k] * kg_86[k];

        t_207[k] = f_11 * kf_65[k]
                   + pb_y[k] * lf_86[k];

        t_208[k] = pb_z[k] * lf_86[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_x, pb_z, kf_89, kf_90, ld0_27, ld0_28, \
                         ld1_27, ld1_28, lf_87, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_6 * kf_89[k]
                   + f_3 * ld0_28[k]
                   - f_4 * ld1_28[k]
                   + pb_x[k] * lf_88[k];

        t_210[k] = f_3 * ld0_27[k]
                   - f_4 * ld1_27[k]
                   + pb_z[k] * lf_87[k];

        t_211[k] = f_6 * kf_90[k]
                   + pb_x[k] * lf_89[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pb_y, pb_z, ig0_47, ig1_94, kf_70, \
                         kg_122, ld0_28, ld1_28, lf_89, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_9 * ig0_47[k]
                   - f_10 * ig1_94[k]
                   + pa_x[k] * kg_122[k];

        t_213[k] = pb_z[k] * lf_89[k];

        t_214[k] = f_3 * ld0_28[k]
                   - f_4 * ld1_28[k]
                   + pb_z[k] * lf_90[k];

        t_215[k] = f_11 * kf_70[k]
                   + pb_y[k] * lf_91[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pa_z, pb_z, kf_65, kf_66, kg_86, \
                         kg_88, kg_89, ld0_29, ld1_29, lf_91, lf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_1 * ld0_29[k]
                   - f_2 * ld1_29[k]
                   + pb_z[k] * lf_91[k];

        t_217[k] = pa_z[k] * kg_86[k];

        t_218[k] = f_5 * kf_65[k]
                   + pb_z[k] * lf_92[k];

        t_219[k] = pa_z[k] * kg_88[k];

        t_220[k] = f_6 * kf_66[k]
                   + pa_z[k] * kg_89[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_z, pb_y, pb_z, kf_68, kf_69, kf_73, \
                         kg_91, kg_93, lf_93, lf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pa_z[k] * kg_91[k];

        t_222[k] = f_5 * kf_68[k]
                   + pb_z[k] * lf_93[k];

        t_223[k] = f_6 * kf_69[k]
                   + pa_z[k] * kg_93[k];

        t_224[k] = f_17 * kf_73[k]
                   + pb_y[k] * lf_94[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pa_z, pb_z, ig0_27, ig1_58, kf_70, kf_71, \
                         kg_94, kg_96, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_8 * kf_70[k]
                   + pa_z[k] * kg_94[k];

        t_226[k] = f_20 * ig0_27[k]
                   - f_21 * ig1_58[k]
                   + pa_y[k] * kg_96[k];

        t_227[k] = f_6 * kf_71[k]
                   + pb_z[k] * lf_95[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_x, pa_y, pa_z, ig0_23, ig0_29, ig0_49, \
                         ig1_51, ig1_60, ig1_108, kg_95, kg_98, \
                         kg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_9 * ig0_23[k]
                   - f_10 * ig1_51[k]
                   + pa_z[k] * kg_95[k];

        t_229[k] = f_20 * ig0_29[k]
                   - f_21 * ig1_60[k]
                   + pa_y[k] * kg_98[k];

        t_230[k] = f_9 * ig0_49[k]
                   - f_10 * ig1_108[k]
                   + pa_x[k] * kg_123[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pa_x, pb_y, pb_z, ig0_50, ig1_110, kf_72, kf_76, \
                         kg_124, lf_96, lf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_6 * kf_72[k]
                   + pb_z[k] * lf_96[k];

        t_232[k] = f_9 * ig0_50[k]
                   - f_10 * ig1_110[k]
                   + pa_x[k] * kg_124[k];

        t_233[k] = f_8 * kf_76[k]
                   + pb_y[k] * lf_97[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_x, pa_y, pb_z, ig0_33, ig0_52, ig1_64, \
                         ig1_112, kf_74, kg_102, kg_125, lf_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_9 * ig0_52[k]
                   - f_10 * ig1_112[k]
                   + pa_x[k] * kg_125[k];

        t_235[k] = f_14 * ig0_33[k]
                   - f_15 * ig1_64[k]
                   + pa_y[k] * kg_102[k];

        t_236[k] = f_16 * kf_74[k]
                   + pb_z[k] * lf_98[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_x, pa_y, pa_z, ig0_26, ig0_34, ig0_53, \
                         ig1_57, ig1_65, ig1_116, kg_97, kg_104, \
                         kg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_14 * ig0_26[k]
                   - f_15 * ig1_57[k]
                   + pa_z[k] * kg_97[k];

        t_238[k] = f_14 * ig0_34[k]
                   - f_15 * ig1_65[k]
                   + pa_y[k] * kg_104[k];

        t_239[k] = f_9 * ig0_53[k]
                   - f_10 * ig1_116[k]
                   + pa_x[k] * kg_126[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pa_x, pb_y, pb_z, ig0_54, ig1_118, kf_75, kf_79, \
                         kg_127, lf_99, lf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_16 * kf_75[k]
                   + pb_z[k] * lf_99[k];

        t_241[k] = f_9 * ig0_54[k]
                   - f_10 * ig1_118[k]
                   + pa_x[k] * kg_127[k];

        t_242[k] = f_16 * kf_79[k]
                   + pb_y[k] * lf_100[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_x, pa_y, pb_z, ig0_35, ig0_56, ig1_66, \
                         ig1_120, kf_77, kg_108, kg_128, lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_9 * ig0_56[k]
                   - f_10 * ig1_120[k]
                   + pa_x[k] * kg_128[k];

        t_244[k] = f_9 * ig0_35[k]
                   - f_10 * ig1_66[k]
                   + pa_y[k] * kg_108[k];

        t_245[k] = f_8 * kf_77[k]
                   + pb_z[k] * lf_101[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pa_y, pa_z, ig0_28, ig0_36, ig0_57, \
                         ig1_59, ig1_69, ig1_124, kg_103, kg_109, \
                         kg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_20 * ig0_28[k]
                   - f_21 * ig1_59[k]
                   + pa_z[k] * kg_103[k];

        t_247[k] = f_9 * ig0_36[k]
                   - f_10 * ig1_69[k]
                   + pa_y[k] * kg_109[k];

        t_248[k] = f_9 * ig0_57[k]
                   - f_10 * ig1_124[k]
                   + pa_x[k] * kg_129[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pb_y, pb_z, ig0_58, ig1_126, kf_78, kf_81, \
                         kg_130, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_8 * kf_78[k]
                   + pb_z[k] * lf_102[k];

        t_250[k] = f_9 * ig0_58[k]
                   - f_10 * ig1_126[k]
                   + pa_x[k] * kg_130[k];

        t_251[k] = f_6 * kf_81[k]
                   + pb_y[k] * lf_103[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, pa_x, pa_y, ig0_60, ig1_128, \
                         kf_83, kg_110, kg_112, kg_113, kg_114, \
                         kg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_9 * ig0_60[k]
                   - f_10 * ig1_128[k]
                   + pa_x[k] * kg_131[k];

        t_253[k] = pa_y[k] * kg_110[k];

        t_254[k] = pa_y[k] * kg_112[k];

        t_255[k] = f_6 * kf_83[k]
                   + pa_y[k] * kg_113[k];

        t_256[k] = pa_y[k] * kg_114[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_y, pb_y, pb_z, kf_80, kf_85, kf_86, \
                         kf_87, kg_116, kg_117, lf_104, lf_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_8 * kf_85[k]
                   + pa_y[k] * kg_116[k];

        t_258[k] = f_17 * kf_80[k]
                   + pb_z[k] * lf_104[k];

        t_259[k] = f_6 * kf_86[k]
                   + pa_y[k] * kg_117[k];

        t_260[k] = f_5 * kf_87[k]
                   + pb_y[k] * lf_105[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pa_z, pb_y, pb_z, ig0_35, ig1_66, \
                         kf_82, kg_110, kg_119, lf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_y[k] * kg_119[k];

        t_262[k] = f_12 * ig0_35[k]
                   - f_13 * ig1_66[k]
                   + pa_z[k] * kg_110[k];

        t_263[k] = pb_y[k] * lf_106[k];

        t_264[k] = f_11 * kf_82[k]
                   + pb_z[k] * lf_106[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pb_x, pb_y, kf_96, kf_97, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_107, lf_108, lf_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_3 * ld0_30[k]
                   - f_4 * ld1_30[k]
                   + pb_y[k] * lf_107[k];

        t_266[k] = f_6 * kf_96[k]
                   + f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_x[k] * lf_108[k];

        t_267[k] = f_6 * kf_97[k]
                   + pb_x[k] * lf_111[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pb_y, pb_z, kf_85, ld0_31, ld0_32, \
                         ld1_31, ld1_32, lf_109, lf_110, lf_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_109[k];

        t_269[k] = f_11 * kf_85[k]
                   + pb_z[k] * lf_109[k];

        t_270[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_110[k];

        t_271[k] = pb_y[k] * lf_111[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_x, pb_y, ig0_62, ig1_145, kf_88, \
                         kf_98, kf_100, kg_135, kg_136, kg_137, \
                         lf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_9 * ig0_62[k]
                   - f_10 * ig1_145[k]
                   + pa_x[k] * kg_135[k];

        t_273[k] = f_8 * kf_98[k]
                   + pa_x[k] * kg_136[k];

        t_274[k] = f_7 * kf_88[k]
                   + pb_y[k] * lf_112[k];

        t_275[k] = f_6 * kf_100[k]
                   + pa_x[k] * kg_137[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, pa_x, pb_x, kf_101, kf_102, \
                         kg_138, kg_141, kg_143, kg_144, kg_145, \
                         lf_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * kf_101[k]
                   + pa_x[k] * kg_138[k];

        t_277[k] = f_5 * kf_102[k]
                   + pb_x[k] * lf_113[k];

        t_278[k] = pa_x[k] * kg_141[k];

        t_279[k] = pa_x[k] * kg_143[k];

        t_280[k] = pa_x[k] * kg_144[k];

        t_281[k] = pa_x[k] * kg_145[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, pa_x, pa_z, pb_z, kf_88, kf_106, \
                         kg_120, kg_121, kg_146, kg_148, lf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pa_z[k] * kg_120[k];

        t_283[k] = f_5 * kf_88[k]
                   + pb_z[k] * lf_114[k];

        t_284[k] = pa_z[k] * kg_121[k];

        t_285[k] = f_6 * kf_106[k]
                   + pa_x[k] * kg_146[k];

        t_286[k] = pa_x[k] * kg_148[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_x, pb_z, kf_91, kf_109, kg_149, \
                         kg_150, kg_151, kg_152, lf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_x[k] * kg_149[k];

        t_288[k] = pa_x[k] * kg_150[k];

        t_289[k] = pa_x[k] * kg_151[k];

        t_290[k] = f_8 * kf_109[k]
                   + pa_x[k] * kg_152[k];

        t_291[k] = f_6 * kf_91[k]
                   + pb_z[k] * lf_115[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, t_297, pa_x, kf_110, kf_111, \
                         kg_153, kg_154, kg_157, kg_158, kg_159, \
                         kg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_6 * kf_110[k]
                   + pa_x[k] * kg_153[k];

        t_293[k] = f_6 * kf_111[k]
                   + pa_x[k] * kg_154[k];

        t_294[k] = pa_x[k] * kg_157[k];

        t_295[k] = pa_x[k] * kg_158[k];

        t_296[k] = pa_x[k] * kg_159[k];

        t_297[k] = pa_x[k] * kg_160[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pa_x, pb_z, kf_92, kf_115, kf_116, \
                         kf_117, kg_161, kg_162, kg_163, kg_164, \
                         lf_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pa_x[k] * kg_161[k];

        t_299[k] = f_8 * kf_115[k]
                   + pa_x[k] * kg_162[k];

        t_300[k] = f_16 * kf_92[k]
                   + pb_z[k] * lf_116[k];

        t_301[k] = f_6 * kf_116[k]
                   + pa_x[k] * kg_163[k];

        t_302[k] = f_6 * kf_117[k]
                   + pa_x[k] * kg_164[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, t_308, pa_x, kf_121, kg_167, \
                         kg_168, kg_169, kg_170, kg_171, kg_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pa_x[k] * kg_167[k];

        t_304[k] = pa_x[k] * kg_168[k];

        t_305[k] = pa_x[k] * kg_169[k];

        t_306[k] = pa_x[k] * kg_170[k];

        t_307[k] = pa_x[k] * kg_171[k];

        t_308[k] = f_8 * kf_121[k]
                   + pa_x[k] * kg_172[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pa_x, pb_z, kf_93, kf_122, kf_123, \
                         kg_173, kg_174, kg_177, kg_178, lf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_8 * kf_93[k]
                   + pb_z[k] * lf_117[k];

        t_310[k] = f_6 * kf_122[k]
                   + pa_x[k] * kg_173[k];

        t_311[k] = f_6 * kf_123[k]
                   + pa_x[k] * kg_174[k];

        t_312[k] = pa_x[k] * kg_177[k];

        t_313[k] = pa_x[k] * kg_178[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, pa_x, pb_z, kf_94, kf_127, kg_179, \
                         kg_180, kg_181, kg_182, lf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_x[k] * kg_179[k];

        t_315[k] = pa_x[k] * kg_180[k];

        t_316[k] = pa_x[k] * kg_181[k];

        t_317[k] = f_8 * kf_127[k]
                   + pa_x[k] * kg_182[k];

        t_318[k] = f_17 * kf_94[k]
                   + pb_z[k] * lf_118[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, t_324, pa_x, kf_128, kf_129, \
                         kg_183, kg_184, kg_187, kg_188, kg_189, \
                         kg_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_6 * kf_128[k]
                   + pa_x[k] * kg_183[k];

        t_320[k] = f_6 * kf_129[k]
                   + pa_x[k] * kg_184[k];

        t_321[k] = pa_x[k] * kg_187[k];

        t_322[k] = pa_x[k] * kg_188[k];

        t_323[k] = pa_x[k] * kg_189[k];

        t_324[k] = pa_x[k] * kg_190[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, t_330, pa_x, pa_y, kf_133, kg_132, \
                         kg_133, kg_134, kg_191, kg_192, kg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_x[k] * kg_191[k];

        t_326[k] = pa_y[k] * kg_132[k];

        t_327[k] = pa_y[k] * kg_133[k];

        t_328[k] = f_6 * kf_133[k]
                   + pa_x[k] * kg_192[k];

        t_329[k] = pa_y[k] * kg_134[k];

        t_330[k] = pa_x[k] * kg_193[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, pa_x, pb_z, kf_95, kf_137, kg_194, \
                         kg_195, kg_196, kg_198, lf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_x[k] * kg_194[k];

        t_332[k] = pa_x[k] * kg_195[k];

        t_333[k] = pa_x[k] * kg_196[k];

        t_334[k] = f_8 * kf_137[k]
                   + pa_x[k] * kg_198[k];

        t_335[k] = f_7 * kf_95[k]
                   + pb_z[k] * lf_119[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, pa_x, pb_x, kf_139, kf_140, \
                         kf_143, kg_200, kg_201, kg_204, kg_205, \
                         lf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_6 * kf_139[k]
                   + pa_x[k] * kg_200[k];

        t_337[k] = f_6 * kf_140[k]
                   + pa_x[k] * kg_201[k];

        t_338[k] = f_5 * kf_143[k]
                   + pb_x[k] * lf_120[k];

        t_339[k] = pa_x[k] * kg_204[k];

        t_340[k] = pa_x[k] * kg_205[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_x, pb_x, pb_y, kf_98, kg_206, kg_208, \
                         ld0_33, ld1_33, lf_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pa_x[k] * kg_206[k];

        t_342[k] = pa_x[k] * kg_208[k];

        t_343[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_121[k];

        t_344[k] = f_0 * kf_98[k]
                   + pb_y[k] * lf_121[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pb_x, pb_y, kf_102, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_122, lf_123, lf_124, \
                         lf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_122[k];

        t_346[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_123[k];

        t_347[k] = pb_x[k] * lf_124[k];

        t_348[k] = pb_x[k] * lf_126[k];

        t_349[k] = f_0 * kf_102[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_124[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_y, pb_z, kf_104, ld0_34, ld0_35, \
                         ld1_34, ld1_35, lf_124, lf_125, lf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_z[k] * lf_124[k];

        t_351[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_125[k];

        t_352[k] = f_0 * kf_104[k]
                   + pb_y[k] * lf_126[k];

        t_353[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_126[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pa_z, pb_z, kf_98, kf_99, kg_136, \
                         kg_137, kg_138, kg_141, lf_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pa_z[k] * kg_136[k];

        t_355[k] = f_5 * kf_98[k]
                   + pb_z[k] * lf_127[k];

        t_356[k] = pa_z[k] * kg_137[k];

        t_357[k] = f_6 * kf_99[k]
                   + pa_z[k] * kg_138[k];

        t_358[k] = pa_z[k] * kg_141[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_z, pb_y, pb_z, kf_102, kf_103, kf_104, \
                         kf_108, kg_143, kg_145, lf_128, lf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_5 * kf_102[k]
                   + pb_z[k] * lf_128[k];

        t_360[k] = f_6 * kf_103[k]
                   + pa_z[k] * kg_143[k];

        t_361[k] = f_7 * kf_108[k]
                   + pb_y[k] * lf_129[k];

        t_362[k] = f_8 * kf_104[k]
                   + pa_z[k] * kg_145[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_x, pb_z, kf_105, ld0_36, ld0_37, \
                         ld0_38, ld1_36, ld1_37, ld1_38, lf_130, lf_131, \
                         lf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_130[k];

        t_364[k] = f_6 * kf_105[k]
                   + pb_z[k] * lf_130[k];

        t_365[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_131[k];

        t_366[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_132[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_z, pb_x, pb_z, ig0_47, ig1_94, kf_107, \
                         kg_147, lf_133, lf_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pb_x[k] * lf_133[k];

        t_368[k] = pb_x[k] * lf_135[k];

        t_369[k] = f_9 * ig0_47[k]
                   - f_10 * ig1_94[k]
                   + pa_z[k] * kg_147[k];

        t_370[k] = f_6 * kf_107[k]
                   + pb_z[k] * lf_133[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, pa_y, pb_y, ig0_52, ig1_112, kf_113, kf_114, \
                         kg_161, ld0_38, ld1_38, lf_134, lf_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_11 * kf_113[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_134[k];

        t_372[k] = f_11 * kf_114[k]
                   + pb_y[k] * lf_135[k];

        t_373[k] = f_12 * ig0_52[k]
                   - f_13 * ig1_112[k]
                   + pa_y[k] * kg_161[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pb_x, pb_z, kf_109, ld0_39, ld0_40, \
                         ld0_41, ld1_39, ld1_40, ld1_41, lf_136, lf_137, \
                         lf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_136[k];

        t_375[k] = f_16 * kf_109[k]
                   + pb_z[k] * lf_136[k];

        t_376[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_137[k];

        t_377[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_138[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_z, pb_x, pb_z, ig0_48, ig1_100, \
                         kf_112, kg_157, lf_139, lf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_x[k] * lf_139[k];

        t_379[k] = pb_x[k] * lf_141[k];

        t_380[k] = f_14 * ig0_48[k]
                   - f_15 * ig1_100[k]
                   + pa_z[k] * kg_157[k];

        t_381[k] = f_16 * kf_112[k]
                   + pb_z[k] * lf_139[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pa_y, pb_y, ig0_56, ig1_120, kf_119, kf_120, \
                         kg_171, ld0_41, ld1_41, lf_140, lf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_17 * kf_119[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_140[k];

        t_383[k] = f_17 * kf_120[k]
                   + pb_y[k] * lf_141[k];

        t_384[k] = f_18 * ig0_56[k]
                   - f_19 * ig1_120[k]
                   + pa_y[k] * kg_171[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_x, pb_z, kf_115, ld0_42, ld0_43, \
                         ld0_44, ld1_42, ld1_43, ld1_44, lf_142, lf_143, \
                         lf_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_142[k];

        t_386[k] = f_8 * kf_115[k]
                   + pb_z[k] * lf_142[k];

        t_387[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_143[k];

        t_388[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_144[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_z, pb_x, pb_z, ig0_49, ig1_108, \
                         kf_118, kg_167, lf_145, lf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = pb_x[k] * lf_145[k];

        t_390[k] = pb_x[k] * lf_147[k];

        t_391[k] = f_20 * ig0_49[k]
                   - f_21 * ig1_108[k]
                   + pa_z[k] * kg_167[k];

        t_392[k] = f_8 * kf_118[k]
                   + pb_z[k] * lf_145[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, pa_y, pb_y, ig0_60, ig1_128, kf_125, kf_126, \
                         kg_181, ld0_44, ld1_44, lf_146, lf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_8 * kf_125[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_146[k];

        t_394[k] = f_8 * kf_126[k]
                   + pb_y[k] * lf_147[k];

        t_395[k] = f_20 * ig0_60[k]
                   - f_21 * ig1_128[k]
                   + pa_y[k] * kg_181[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pb_x, pb_z, kf_121, ld0_45, ld0_46, \
                         ld0_47, ld1_45, ld1_46, ld1_47, lf_148, lf_149, \
                         lf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_148[k];

        t_397[k] = f_17 * kf_121[k]
                   + pb_z[k] * lf_148[k];

        t_398[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_149[k];

        t_399[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_150[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, pa_z, pb_x, pb_z, ig0_53, ig1_116, \
                         kf_124, kg_177, lf_151, lf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pb_x[k] * lf_151[k];

        t_401[k] = pb_x[k] * lf_153[k];

        t_402[k] = f_18 * ig0_53[k]
                   - f_19 * ig1_116[k]
                   + pa_z[k] * kg_177[k];

        t_403[k] = f_17 * kf_124[k]
                   + pb_z[k] * lf_151[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_y, pb_y, ig0_61, ig1_134, kf_131, kf_132, \
                         kg_191, ld0_47, ld1_47, lf_152, lf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_16 * kf_131[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_152[k];

        t_405[k] = f_16 * kf_132[k]
                   + pb_y[k] * lf_153[k];

        t_406[k] = f_14 * ig0_61[k]
                   - f_15 * ig1_134[k]
                   + pa_y[k] * kg_191[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pb_x, pb_z, kf_127, ld0_48, ld0_49, \
                         ld0_50, ld1_48, ld1_49, ld1_50, lf_154, lf_155, \
                         lf_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_154[k];

        t_408[k] = f_11 * kf_127[k]
                   + pb_z[k] * lf_154[k];

        t_409[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_155[k];

        t_410[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_156[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, pa_z, pb_x, pb_z, ig0_57, ig1_124, \
                         kf_130, kg_187, lf_157, lf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = pb_x[k] * lf_157[k];

        t_412[k] = pb_x[k] * lf_159[k];

        t_413[k] = f_12 * ig0_57[k]
                   - f_13 * ig1_124[k]
                   + pa_z[k] * kg_187[k];

        t_414[k] = f_11 * kf_130[k]
                   + pb_z[k] * lf_157[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_y, pb_y, ig0_62, ig1_145, kf_135, \
                         kf_136, kg_197, kg_198, ld0_50, ld1_50, lf_158, \
                         lf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_6 * kf_135[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_158[k];

        t_416[k] = f_6 * kf_136[k]
                   + pb_y[k] * lf_159[k];

        t_417[k] = f_9 * ig0_62[k]
                   - f_10 * ig1_145[k]
                   + pa_y[k] * kg_197[k];

        t_418[k] = pa_y[k] * kg_198[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, pa_y, pb_z, kf_134, kf_138, \
                         kf_141, kg_199, kg_200, kg_201, kg_204, \
                         lf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pa_y[k] * kg_199[k];

        t_420[k] = f_6 * kf_138[k]
                   + pa_y[k] * kg_200[k];

        t_421[k] = pa_y[k] * kg_201[k];

        t_422[k] = f_8 * kf_141[k]
                   + pa_y[k] * kg_204[k];

        t_423[k] = f_7 * kf_134[k]
                   + pb_z[k] * lf_160[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pa_y, pb_x, pb_y, kf_142, kf_143, kg_206, \
                         kg_208, ld0_51, ld1_51, lf_161, lf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_6 * kf_142[k]
                   + pa_y[k] * kg_206[k];

        t_425[k] = f_5 * kf_143[k]
                   + pb_y[k] * lf_161[k];

        t_426[k] = pa_y[k] * kg_208[k];

        t_427[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_162[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pb_x, pb_z, kf_137, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_162, lf_163, lf_164, \
                         lf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * kf_137[k]
                   + pb_z[k] * lf_162[k];

        t_429[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_163[k];

        t_430[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_164[k];

        t_431[k] = pb_x[k] * lf_165[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, pb_x, pb_y, pb_z, kf_141, ld0_52, \
                         ld0_53, ld1_52, ld1_53, lf_165, lf_166, \
                         lf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = pb_x[k] * lf_167[k];

        t_433[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_165[k];

        t_434[k] = f_0 * kf_141[k]
                   + pb_z[k] * lf_165[k];

        t_435[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_166[k];

        t_436[k] = pb_y[k] * lf_167[k];
    }

#pragma omp simd aligned(t_437, pb_z, kf_143, ld0_53, ld1_53, lf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_0 * kf_143[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_167[k];
    }
}

auto
compute_prim_lg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_1 = buffer.data(ig0 + 1);
    const auto *ig0_2 = buffer.data(ig0 + 2);
    const auto *ig0_3 = buffer.data(ig0 + 3);
    const auto *ig0_6 = buffer.data(ig0 + 6);
    const auto *ig0_7 = buffer.data(ig0 + 7);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_14 = buffer.data(ig0 + 14);
    const auto *ig0_15 = buffer.data(ig0 + 15);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_19 = buffer.data(ig0 + 19);
    const auto *ig0_22 = buffer.data(ig0 + 22);
    const auto *ig0_23 = buffer.data(ig0 + 23);
    const auto *ig0_24 = buffer.data(ig0 + 24);
    const auto *ig0_27 = buffer.data(ig0 + 27);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_31 = buffer.data(ig0 + 31);
    const auto *ig0_32 = buffer.data(ig0 + 32);
    const auto *ig0_33 = buffer.data(ig0 + 33);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_35 = buffer.data(ig0 + 35);
    const auto *ig0_37 = buffer.data(ig0 + 37);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_39 = buffer.data(ig0 + 39);
    const auto *ig0_41 = buffer.data(ig0 + 41);
    const auto *ig0_42 = buffer.data(ig0 + 42);
    const auto *ig0_43 = buffer.data(ig0 + 43);
    const auto *ig0_45 = buffer.data(ig0 + 45);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_1 = buffer.data(ig1 + 1);
    const auto *ig1_2 = buffer.data(ig1 + 2);
    const auto *ig1_3 = buffer.data(ig1 + 3);
    const auto *ig1_6 = buffer.data(ig1 + 6);
    const auto *ig1_7 = buffer.data(ig1 + 7);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_14 = buffer.data(ig1 + 14);
    const auto *ig1_15 = buffer.data(ig1 + 15);
    const auto *ig1_18 = buffer.data(ig1 + 18);
    const auto *ig1_19 = buffer.data(ig1 + 19);
    const auto *ig1_22 = buffer.data(ig1 + 22);
    const auto *ig1_23 = buffer.data(ig1 + 23);
    const auto *ig1_24 = buffer.data(ig1 + 24);
    const auto *ig1_27 = buffer.data(ig1 + 27);
    const auto *ig1_28 = buffer.data(ig1 + 28);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_30 = buffer.data(ig1 + 30);
    const auto *ig1_31 = buffer.data(ig1 + 31);
    const auto *ig1_32 = buffer.data(ig1 + 32);
    const auto *ig1_33 = buffer.data(ig1 + 33);
    const auto *ig1_34 = buffer.data(ig1 + 34);
    const auto *ig1_35 = buffer.data(ig1 + 35);
    const auto *ig1_37 = buffer.data(ig1 + 37);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_39 = buffer.data(ig1 + 39);
    const auto *ig1_41 = buffer.data(ig1 + 41);
    const auto *ig1_42 = buffer.data(ig1 + 42);
    const auto *ig1_43 = buffer.data(ig1 + 43);
    const auto *ig1_45 = buffer.data(ig1 + 45);
    const auto *ig1_46 = buffer.data(ig1 + 46);
    const auto *ig1_47 = buffer.data(ig1 + 47);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_95 = buffer.data(kf + 95);

    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_8, kg_9, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_9[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_8[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_6, ig1_6, kf_9, kg_16, \
                         ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_9[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_6[k]
                  - f_9 * ig1_6[k]
                  + pa_x[k] * kg_16[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_10, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_14, kf_17, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_14[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_17[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_10, ig1_10, kg_28, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_10[k]
                  - f_9 * ig1_10[k]
                  + pa_x[k] * kg_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_1, ig1_1, kf_20, kg_11, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_1[k]
                  - f_11 * ig1_1[k]
                  + pa_y[k] * kg_11[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_20[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_14, ig1_14, kf_21, \
                         kg_34, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_21[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_14[k]
                  - f_14 * ig1_14[k]
                  + pa_x[k] * kg_34[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_2, ig1_2, kg_20, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_2[k]
                  - f_11 * ig1_2[k]
                  + pa_z[k] * kg_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_26, kf_29, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_26[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_29[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_18, ig1_18, kg_46, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_18[k]
                  - f_14 * ig1_18[k]
                  + pa_x[k] * kg_46[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_3, ig1_3, kf_32, kg_29, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_3[k]
                  - f_16 * ig1_3[k]
                  + pa_y[k] * kg_29[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_32[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_22, ig1_22, kf_33, \
                         kg_52, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_33[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_22[k]
                  - f_16 * ig1_22[k]
                  + pa_x[k] * kg_52[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_23, ig1_23, kg_56, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_23[k]
                  - f_16 * ig1_23[k]
                  + pa_x[k] * kg_56[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_7, ig1_7, kg_38, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_7[k]
                  - f_16 * ig1_7[k]
                  + pa_z[k] * kg_38[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_38, kf_41, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_38[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_41[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_11, ig0_27, \
                         ig1_11, ig1_27, kg_47, kg_65, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_27[k]
                  - f_16 * ig1_27[k]
                  + pa_x[k] * kg_65[k];

        t_64[k] = f_13 * ig0_11[k]
                  - f_14 * ig1_11[k]
                  + pa_y[k] * kg_47[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_44, kf_45, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_44[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_45[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_28, ig1_28, kg_71, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_28[k]
                  - f_11 * ig1_28[k]
                  + pa_x[k] * kg_71[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_15, ig0_29, ig0_30, ig1_15, ig1_29, \
                         ig1_30, kg_57, kg_75, kg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_29[k]
                  - f_11 * ig1_29[k]
                  + pa_x[k] * kg_75[k];

        t_74[k] = f_10 * ig0_30[k]
                  - f_11 * ig1_30[k]
                  + pa_x[k] * kg_76[k];

        t_75[k] = f_13 * ig0_15[k]
                  - f_14 * ig1_15[k]
                  + pa_z[k] * kg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_50, kf_53, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_50[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_53[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_31, ig1_31, kg_85, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_31[k]
                  - f_11 * ig1_31[k]
                  + pa_x[k] * kg_85[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_19, ig1_19, kf_54, kg_66, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_19[k]
                  - f_9 * ig1_19[k]
                  + pa_y[k] * kg_66[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_54[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_32, ig1_32, kf_55, \
                         kg_86, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_55[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_32[k]
                  - f_6 * ig1_32[k]
                  + pa_x[k] * kg_86[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_35, ig1_35, kg_87, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_35[k]
                  - f_6 * ig1_35[k]
                  + pa_x[k] * kg_87[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_24, ig0_39, ig0_43, ig1_24, ig1_39, \
                         ig1_43, kg_77, kg_88, kg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_39[k]
                  - f_6 * ig1_39[k]
                  + pa_x[k] * kg_88[k];

        t_95[k] = f_5 * ig0_43[k]
                  - f_6 * ig1_43[k]
                  + pa_x[k] * kg_89[k];

        t_96[k] = f_8 * ig0_24[k]
                  - f_9 * ig1_24[k]
                  + pa_z[k] * kg_77[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_56, kf_57, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_56[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_57[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_47, ig1_47, kg_90, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_47[k]
                   - f_6 * ig1_47[k]
                   + pa_x[k] * kg_90[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_61, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_61[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_32, ig1_32, kf_68, \
                         kf_69, kg_100, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_32[k]
                   - f_6 * ig1_32[k]
                   + pa_z[k] * kg_100[k];

        t_120[k] = f_7 * kf_68[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_69[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_37, ig1_37, kg_109, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_37[k]
                   - f_9 * ig1_37[k]
                   + pa_y[k] * kg_109[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_33, ig1_33, kg_106, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_33[k]
                   - f_11 * ig1_33[k]
                   + pa_z[k] * kg_106[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_41, ig1_41, kf_74, kf_75, \
                         kg_118, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_74[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_75[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_41[k]
                   - f_14 * ig1_41[k]
                   + pa_y[k] * kg_118[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_34, ig1_34, kf_80, \
                         kf_81, kg_115, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_34[k]
                   - f_16 * ig1_34[k]
                   + pa_z[k] * kg_115[k];

        t_138[k] = f_17 * kf_80[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_81[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_45, ig1_45, kg_127, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_45[k]
                   - f_16 * ig1_45[k]
                   + pa_y[k] * kg_127[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_38, ig1_38, kg_124, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_38[k]
                   - f_14 * ig1_38[k]
                   + pa_z[k] * kg_124[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_46, ig1_46, kf_86, kf_87, \
                         kg_136, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_86[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_87[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_46[k]
                   - f_11 * ig1_46[k]
                   + pa_y[k] * kg_136[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_42, ig1_42, kf_88, \
                         kf_89, kg_133, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_42[k]
                   - f_9 * ig1_42[k]
                   + pa_z[k] * kg_133[k];

        t_156[k] = f_19 * kf_88[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_89[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_47, ig1_47, kg_137, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_47[k]
                   - f_6 * ig1_47[k]
                   + pa_y[k] * kg_137[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_95, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_95[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_1 = buffer.data(ig0 + 1);
    const auto *ig0_2 = buffer.data(ig0 + 2);
    const auto *ig0_3 = buffer.data(ig0 + 3);
    const auto *ig0_6 = buffer.data(ig0 + 6);
    const auto *ig0_7 = buffer.data(ig0 + 7);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_14 = buffer.data(ig0 + 14);
    const auto *ig0_15 = buffer.data(ig0 + 15);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_19 = buffer.data(ig0 + 19);
    const auto *ig0_22 = buffer.data(ig0 + 22);
    const auto *ig0_23 = buffer.data(ig0 + 23);
    const auto *ig0_24 = buffer.data(ig0 + 24);
    const auto *ig0_27 = buffer.data(ig0 + 27);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_31 = buffer.data(ig0 + 31);
    const auto *ig0_32 = buffer.data(ig0 + 32);
    const auto *ig0_33 = buffer.data(ig0 + 33);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_35 = buffer.data(ig0 + 35);
    const auto *ig0_37 = buffer.data(ig0 + 37);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_39 = buffer.data(ig0 + 39);
    const auto *ig0_41 = buffer.data(ig0 + 41);
    const auto *ig0_42 = buffer.data(ig0 + 42);
    const auto *ig0_43 = buffer.data(ig0 + 43);
    const auto *ig0_45 = buffer.data(ig0 + 45);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_13 = buffer.data(ig1 + 13);
    const auto *ig1_18 = buffer.data(ig1 + 18);
    const auto *ig1_23 = buffer.data(ig1 + 23);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_37 = buffer.data(ig1 + 37);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_43 = buffer.data(ig1 + 43);
    const auto *ig1_53 = buffer.data(ig1 + 53);
    const auto *ig1_61 = buffer.data(ig1 + 61);
    const auto *ig1_62 = buffer.data(ig1 + 62);
    const auto *ig1_67 = buffer.data(ig1 + 67);
    const auto *ig1_76 = buffer.data(ig1 + 76);
    const auto *ig1_81 = buffer.data(ig1 + 81);
    const auto *ig1_89 = buffer.data(ig1 + 89);
    const auto *ig1_93 = buffer.data(ig1 + 93);
    const auto *ig1_97 = buffer.data(ig1 + 97);
    const auto *ig1_101 = buffer.data(ig1 + 101);
    const auto *ig1_106 = buffer.data(ig1 + 106);
    const auto *ig1_113 = buffer.data(ig1 + 113);
    const auto *ig1_120 = buffer.data(ig1 + 120);
    const auto *ig1_127 = buffer.data(ig1 + 127);
    const auto *ig1_128 = buffer.data(ig1 + 128);
    const auto *ig1_130 = buffer.data(ig1 + 130);
    const auto *ig1_136 = buffer.data(ig1 + 136);
    const auto *ig1_137 = buffer.data(ig1 + 137);
    const auto *ig1_139 = buffer.data(ig1 + 139);
    const auto *ig1_145 = buffer.data(ig1 + 145);
    const auto *ig1_146 = buffer.data(ig1 + 146);
    const auto *ig1_148 = buffer.data(ig1 + 148);
    const auto *ig1_154 = buffer.data(ig1 + 154);
    const auto *ig1_164 = buffer.data(ig1 + 164);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_98 = buffer.data(kf + 98);

    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_19 = buffer.data(kg + 19);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_148 = buffer.data(kg + 148);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_165 = buffer.data(kg + 165);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_186 = buffer.data(kg + 186);
    const auto *kg_192 = buffer.data(kg + 192);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_10, kg_10, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_10[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_10[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_6, ig1_23, kf_11, \
                         kg_24, ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_11[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_6[k]
                  - f_9 * ig1_23[k]
                  + pa_x[k] * kg_24[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_13, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_16, kf_19, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_16[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_19[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_10, ig1_37, kg_38, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_10[k]
                  - f_9 * ig1_37[k]
                  + pa_x[k] * kg_38[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_1, ig1_10, kf_22, kg_19, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_1[k]
                  - f_11 * ig1_10[k]
                  + pa_y[k] * kg_19[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_22[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_14, ig1_43, kf_23, \
                         kg_44, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_23[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_14[k]
                  - f_14 * ig1_43[k]
                  + pa_x[k] * kg_44[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_2, ig1_13, kg_30, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_2[k]
                  - f_11 * ig1_13[k]
                  + pa_z[k] * kg_30[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_28, kf_31, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_28[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_18, ig1_61, kg_61, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_18[k]
                  - f_14 * ig1_61[k]
                  + pa_x[k] * kg_61[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_3, ig1_18, kf_34, kg_39, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_3[k]
                  - f_16 * ig1_18[k]
                  + pa_y[k] * kg_39[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_34[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_22, ig1_67, kf_35, \
                         kg_67, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_35[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_22[k]
                  - f_16 * ig1_67[k]
                  + pa_x[k] * kg_67[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_23, ig1_76, kg_76, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_23[k]
                  - f_16 * ig1_76[k]
                  + pa_x[k] * kg_76[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_7, ig1_29, kg_53, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_7[k]
                  - f_16 * ig1_29[k]
                  + pa_z[k] * kg_53[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_40, kf_43, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_40[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_11, ig0_27, \
                         ig1_38, ig1_89, kg_62, kg_88, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_27[k]
                  - f_16 * ig1_89[k]
                  + pa_x[k] * kg_88[k];

        t_64[k] = f_13 * ig0_11[k]
                  - f_14 * ig1_38[k]
                  + pa_y[k] * kg_62[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_46, kf_47, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_46[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_47[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_28, ig1_93, kg_94, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_28[k]
                  - f_11 * ig1_93[k]
                  + pa_x[k] * kg_94[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_15, ig0_29, ig0_30, ig1_53, ig1_97, \
                         ig1_101, kg_80, kg_103, kg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_29[k]
                  - f_11 * ig1_97[k]
                  + pa_x[k] * kg_103[k];

        t_74[k] = f_10 * ig0_30[k]
                  - f_11 * ig1_101[k]
                  + pa_x[k] * kg_107[k];

        t_75[k] = f_13 * ig0_15[k]
                  - f_14 * ig1_53[k]
                  + pa_z[k] * kg_80[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_52, kf_55, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_52[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_55[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_31, ig1_106, kg_119, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_31[k]
                  - f_11 * ig1_106[k]
                  + pa_x[k] * kg_119[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_19, ig1_62, kf_56, kg_89, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_19[k]
                  - f_9 * ig1_62[k]
                  + pa_y[k] * kg_89[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_56[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_32, ig1_113, kf_57, \
                         kg_123, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_32[k]
                  - f_6 * ig1_113[k]
                  + pa_x[k] * kg_123[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_35, ig1_128, kg_126, ld0_28, \
                         ld0_29, ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_35[k]
                  - f_6 * ig1_128[k]
                  + pa_x[k] * kg_126[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_24, ig0_39, ig0_43, ig1_81, \
                         ig1_137, ig1_146, kg_111, kg_128, kg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_39[k]
                  - f_6 * ig1_137[k]
                  + pa_x[k] * kg_128[k];

        t_95[k] = f_5 * ig0_43[k]
                  - f_6 * ig1_146[k]
                  + pa_x[k] * kg_130[k];

        t_96[k] = f_8 * ig0_24[k]
                  - f_9 * ig1_81[k]
                  + pa_z[k] * kg_111[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_58, kf_59, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_58[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_59[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_47, ig1_164, kg_134, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_47[k]
                   - f_6 * ig1_164[k]
                   + pa_x[k] * kg_134[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_63, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_63[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_32, ig1_113, kf_71, \
                         kf_72, kg_148, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_32[k]
                   - f_6 * ig1_113[k]
                   + pa_z[k] * kg_148[k];

        t_120[k] = f_7 * kf_71[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_72[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_37, ig1_130, kg_159, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_37[k]
                   - f_9 * ig1_130[k]
                   + pa_y[k] * kg_159[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_33, ig1_120, kg_156, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_33[k]
                   - f_11 * ig1_120[k]
                   + pa_z[k] * kg_156[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_41, ig1_139, kf_77, kf_78, \
                         kg_168, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_77[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_78[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_41[k]
                   - f_14 * ig1_139[k]
                   + pa_y[k] * kg_168[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_34, ig1_127, kf_83, \
                         kf_84, kg_165, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_34[k]
                   - f_16 * ig1_127[k]
                   + pa_z[k] * kg_165[k];

        t_138[k] = f_17 * kf_83[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_84[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_45, ig1_148, kg_177, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_45[k]
                   - f_16 * ig1_148[k]
                   + pa_y[k] * kg_177[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_38, ig1_136, kg_174, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_38[k]
                   - f_14 * ig1_136[k]
                   + pa_z[k] * kg_174[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_46, ig1_154, kf_89, kf_90, \
                         kg_186, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_89[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_90[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_46[k]
                   - f_11 * ig1_154[k]
                   + pa_y[k] * kg_186[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_42, ig1_145, kf_91, \
                         kf_92, kg_183, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_42[k]
                   - f_9 * ig1_145[k]
                   + pa_z[k] * kg_183[k];

        t_156[k] = f_19 * kf_91[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_92[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_47, ig1_164, kg_192, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_47[k]
                   - f_6 * ig1_164[k]
                   + pa_y[k] * kg_192[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_98, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_98[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 2.0 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 3.0 / p;
    const auto f_9 = 2.5 / alpha;
    const auto f_10 = 2.5 * beta / (alpha * p);
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.5 / p;
    const auto f_14 = 2.0 / alpha;
    const auto f_15 = 2.0 * beta / (alpha * p);
    const auto f_16 = 1.5 / alpha;
    const auto f_17 = 1.5 * beta / (alpha * p);
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_13 = buffer.data(ig0 + 13);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_23 = buffer.data(ig0 + 23);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_37 = buffer.data(ig0 + 37);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_43 = buffer.data(ig0 + 43);
    const auto *ig0_50 = buffer.data(ig0 + 50);
    const auto *ig0_53 = buffer.data(ig0 + 53);
    const auto *ig0_61 = buffer.data(ig0 + 61);
    const auto *ig0_62 = buffer.data(ig0 + 62);
    const auto *ig0_67 = buffer.data(ig0 + 67);
    const auto *ig0_74 = buffer.data(ig0 + 74);
    const auto *ig0_75 = buffer.data(ig0 + 75);
    const auto *ig0_76 = buffer.data(ig0 + 76);
    const auto *ig0_77 = buffer.data(ig0 + 77);
    const auto *ig0_78 = buffer.data(ig0 + 78);
    const auto *ig0_81 = buffer.data(ig0 + 81);
    const auto *ig0_89 = buffer.data(ig0 + 89);
    const auto *ig0_93 = buffer.data(ig0 + 93);
    const auto *ig0_96 = buffer.data(ig0 + 96);
    const auto *ig0_97 = buffer.data(ig0 + 97);
    const auto *ig0_98 = buffer.data(ig0 + 98);
    const auto *ig0_100 = buffer.data(ig0 + 100);
    const auto *ig0_101 = buffer.data(ig0 + 101);
    const auto *ig0_102 = buffer.data(ig0 + 102);
    const auto *ig0_106 = buffer.data(ig0 + 106);
    const auto *ig0_113 = buffer.data(ig0 + 113);
    const auto *ig0_120 = buffer.data(ig0 + 120);
    const auto *ig0_127 = buffer.data(ig0 + 127);
    const auto *ig0_128 = buffer.data(ig0 + 128);
    const auto *ig0_130 = buffer.data(ig0 + 130);
    const auto *ig0_136 = buffer.data(ig0 + 136);
    const auto *ig0_137 = buffer.data(ig0 + 137);
    const auto *ig0_139 = buffer.data(ig0 + 139);
    const auto *ig0_145 = buffer.data(ig0 + 145);
    const auto *ig0_146 = buffer.data(ig0 + 146);
    const auto *ig0_148 = buffer.data(ig0 + 148);
    const auto *ig0_154 = buffer.data(ig0 + 154);
    const auto *ig0_164 = buffer.data(ig0 + 164);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_12 = buffer.data(ig1 + 12);
    const auto *ig1_14 = buffer.data(ig1 + 14);
    const auto *ig1_19 = buffer.data(ig1 + 19);
    const auto *ig1_23 = buffer.data(ig1 + 23);
    const auto *ig1_31 = buffer.data(ig1 + 31);
    const auto *ig1_32 = buffer.data(ig1 + 32);
    const auto *ig1_37 = buffer.data(ig1 + 37);
    const auto *ig1_41 = buffer.data(ig1 + 41);
    const auto *ig1_42 = buffer.data(ig1 + 42);
    const auto *ig1_50 = buffer.data(ig1 + 50);
    const auto *ig1_51 = buffer.data(ig1 + 51);
    const auto *ig1_56 = buffer.data(ig1 + 56);
    const auto *ig1_60 = buffer.data(ig1 + 60);
    const auto *ig1_61 = buffer.data(ig1 + 61);
    const auto *ig1_62 = buffer.data(ig1 + 62);
    const auto *ig1_63 = buffer.data(ig1 + 63);
    const auto *ig1_64 = buffer.data(ig1 + 64);
    const auto *ig1_65 = buffer.data(ig1 + 65);
    const auto *ig1_73 = buffer.data(ig1 + 73);
    const auto *ig1_77 = buffer.data(ig1 + 77);
    const auto *ig1_78 = buffer.data(ig1 + 78);
    const auto *ig1_79 = buffer.data(ig1 + 79);
    const auto *ig1_80 = buffer.data(ig1 + 80);
    const auto *ig1_81 = buffer.data(ig1 + 81);
    const auto *ig1_82 = buffer.data(ig1 + 82);
    const auto *ig1_83 = buffer.data(ig1 + 83);
    const auto *ig1_87 = buffer.data(ig1 + 87);
    const auto *ig1_94 = buffer.data(ig1 + 94);
    const auto *ig1_98 = buffer.data(ig1 + 98);
    const auto *ig1_105 = buffer.data(ig1 + 105);
    const auto *ig1_106 = buffer.data(ig1 + 106);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_114 = buffer.data(ig1 + 114);
    const auto *ig1_115 = buffer.data(ig1 + 115);
    const auto *ig1_117 = buffer.data(ig1 + 117);
    const auto *ig1_123 = buffer.data(ig1 + 123);
    const auto *ig1_124 = buffer.data(ig1 + 124);
    const auto *ig1_126 = buffer.data(ig1 + 126);
    const auto *ig1_130 = buffer.data(ig1 + 130);
    const auto *ig1_140 = buffer.data(ig1 + 140);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_101 = buffer.data(kf + 101);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_8 = buffer.data(kg + 8);
    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_18 = buffer.data(kg + 18);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_27 = buffer.data(kg + 27);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_36 = buffer.data(kg + 36);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_40 = buffer.data(kg + 40);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_64 = buffer.data(kg + 64);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_81 = buffer.data(kg + 81);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
    const auto *kg_84 = buffer.data(kg + 84);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_99 = buffer.data(kg + 99);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_161 = buffer.data(kg + 161);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_170 = buffer.data(kg + 170);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, kg_0, ld0_1, ld0_2, ld1_1, \
                         ld1_2, lf_3, lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];

        t_9[k] = pa_y[k] * kg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, ig0_0, ig1_0, kf_3, \
                         kf_5, kg_0, kg_5, kg_8, kg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * kf_3[k]
                  + pa_y[k] * kg_5[k];

        t_11[k] = pa_y[k] * kg_8[k];

        t_12[k] = pa_z[k] * kg_0[k];

        t_13[k] = pa_z[k] * kg_5[k];

        t_14[k] = f_5 * kf_5[k]
                  + pa_z[k] * kg_8[k];

        t_15[k] = f_6 * ig0_0[k]
                  - f_7 * ig1_0[k]
                  + pa_y[k] * kg_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_z, kf_11, kf_12, ld0_3, ld0_4, \
                         ld1_3, ld1_4, lf_6, lf_7, lf_8, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * lf_6[k];

        t_17[k] = f_8 * kf_11[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];

        t_18[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_19[k] = f_8 * kf_12[k]
                  + pb_x[k] * lf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_z, ig0_23, ig1_19, kg_18, ld0_4, \
                         ld0_5, ld1_4, ld1_5, lf_9, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_9 * ig0_23[k]
                  - f_10 * ig1_19[k]
                  + pa_x[k] * kg_18[k];

        t_21[k] = pb_z[k] * lf_9[k];

        t_22[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_23[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pa_z, pb_y, ig0_0, ig1_0, kg_10, kg_11, \
                         kg_12, lf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_z[k] * kg_10[k];

        t_25[k] = pa_y[k] * kg_12[k];

        t_26[k] = f_6 * ig0_0[k]
                  - f_7 * ig1_0[k]
                  + pa_z[k] * kg_11[k];

        t_27[k] = pb_y[k] * lf_12[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, pb_y, kf_17, kf_20, ld0_6, ld0_8, ld1_6, \
                         ld1_8, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_29[k] = f_8 * kf_17[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_30[k] = f_8 * kf_20[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_y, ig0_37, ig1_31, kg_30, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_32[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_33[k] = pb_y[k] * lf_17[k];

        t_34[k] = f_9 * ig0_37[k]
                  - f_10 * ig1_31[k]
                  + pa_x[k] * kg_30[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_x, pb_z, ig0_10, ig1_10, kf_23, kg_13, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_11 * ig0_10[k]
                  - f_12 * ig1_10[k]
                  + pa_y[k] * kg_13[k];

        t_36[k] = pb_z[k] * lf_18[k];

        t_37[k] = f_13 * kf_23[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_x, pb_z, ig0_43, ig1_37, kf_24, \
                         kg_36, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_39[k] = f_13 * kf_24[k]
                  + pb_x[k] * lf_21[k];

        t_40[k] = f_14 * ig0_43[k]
                  - f_15 * ig1_37[k]
                  + pa_x[k] * kg_36[k];

        t_41[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_z, kg_13, kg_18, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_43[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_44[k] = pa_z[k] * kg_13[k];

        t_45[k] = pa_z[k] * kg_18[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, ig0_13, ig1_12, kf_14, kf_18, \
                         kg_21, kg_22, kg_27, kg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_5 * kf_14[k]
                  + pa_z[k] * kg_21[k];

        t_47[k] = f_5 * kf_18[k]
                  + pa_y[k] * kg_27[k];

        t_48[k] = pa_y[k] * kg_30[k];

        t_49[k] = f_11 * ig0_13[k]
                  - f_12 * ig1_12[k]
                  + pa_z[k] * kg_22[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pb_y, kf_29, kf_32, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_y[k] * lf_24[k];

        t_51[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_52[k] = f_13 * kf_29[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_53[k] = f_13 * kf_32[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_y, ig0_61, ig1_50, kg_49, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_55[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_56[k] = pb_y[k] * lf_29[k];

        t_57[k] = f_14 * ig0_61[k]
                  - f_15 * ig1_50[k]
                  + pa_x[k] * kg_49[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pb_x, pb_z, ig0_18, ig1_14, kf_35, kg_31, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_16 * ig0_18[k]
                  - f_17 * ig1_14[k]
                  + pa_y[k] * kg_31[k];

        t_59[k] = pb_z[k] * lf_30[k];

        t_60[k] = f_5 * kf_35[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, ig0_67, ig1_56, kf_36, \
                         kg_55, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_62[k] = f_5 * kf_36[k]
                  + pb_x[k] * lf_33[k];

        t_63[k] = f_16 * ig0_67[k]
                  - f_17 * ig1_56[k]
                  + pa_x[k] * kg_55[k];

        t_64[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_z, kg_31, kg_36, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_66[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_67[k] = pa_z[k] * kg_31[k];

        t_68[k] = pa_z[k] * kg_36[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_x, pa_y, pa_z, ig0_29, ig0_75, ig1_23, ig1_61, \
                         kf_26, kg_39, kg_40, kg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_5 * kf_26[k]
                  + pa_z[k] * kg_39[k];

        t_70[k] = f_6 * ig0_29[k]
                  - f_7 * ig1_23[k]
                  + pa_y[k] * kg_40[k];

        t_71[k] = f_16 * ig0_75[k]
                  - f_17 * ig1_61[k]
                  + pa_x[k] * kg_60[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_x, pa_y, ig0_76, ig0_77, ig1_62, ig1_63, \
                         kf_30, kg_46, kg_49, kg_61, kg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_16 * ig0_76[k]
                  - f_17 * ig1_62[k]
                  + pa_x[k] * kg_61[k];

        t_73[k] = f_16 * ig0_77[k]
                  - f_17 * ig1_63[k]
                  + pa_x[k] * kg_62[k];

        t_74[k] = f_5 * kf_30[k]
                  + pa_y[k] * kg_46[k];

        t_75[k] = pa_y[k] * kg_49[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pb_y, ig0_29, ig1_23, kg_41, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_16 * ig0_29[k]
                  - f_17 * ig1_23[k]
                  + pa_z[k] * kg_41[k];

        t_77[k] = pb_y[k] * lf_36[k];

        t_78[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pb_y, kf_41, kf_44, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * kf_41[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_80[k] = f_5 * kf_44[k]
                  + pb_x[k] * lf_41[k];

        t_81[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_82[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pa_y, pb_y, pb_z, ig0_38, ig0_89, \
                         ig1_32, ig1_73, kg_50, kg_72, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * lf_41[k];

        t_84[k] = f_16 * ig0_89[k]
                  - f_17 * ig1_73[k]
                  + pa_x[k] * kg_72[k];

        t_85[k] = f_14 * ig0_38[k]
                  - f_15 * ig1_32[k]
                  + pa_y[k] * kg_50[k];

        t_86[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_z, kf_47, kf_48, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_18 * kf_47[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_88[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_89[k] = f_18 * kf_48[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pb_z, ig0_93, ig1_77, kg_78, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_11 * ig0_93[k]
                  - f_12 * ig1_77[k]
                  + pa_x[k] * kg_78[k];

        t_91[k] = pb_z[k] * lf_45[k];

        t_92[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_93[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pa_y, pa_z, ig0_50, ig1_41, kf_38, kg_50, \
                         kg_55, kg_58, kg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pa_z[k] * kg_50[k];

        t_95[k] = pa_z[k] * kg_55[k];

        t_96[k] = f_5 * kf_38[k]
                  + pa_z[k] * kg_58[k];

        t_97[k] = f_11 * ig0_50[k]
                  - f_12 * ig1_41[k]
                  + pa_y[k] * kg_59[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_x, ig0_96, ig0_97, ig0_98, ig1_78, ig1_79, \
                         ig1_80, kg_83, kg_84, kg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_11 * ig0_96[k]
                  - f_12 * ig1_78[k]
                  + pa_x[k] * kg_83[k];

        t_99[k] = f_11 * ig0_97[k]
                  - f_12 * ig1_79[k]
                  + pa_x[k] * kg_84[k];

        t_100[k] = f_11 * ig0_98[k]
                   - f_12 * ig1_80[k]
                   + pa_x[k] * kg_85[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_x, pa_y, ig0_53, ig0_100, ig0_101, ig1_42, \
                         ig1_81, ig1_82, kg_63, kg_87, kg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_6 * ig0_53[k]
                   - f_7 * ig1_42[k]
                   + pa_y[k] * kg_63[k];

        t_102[k] = f_11 * ig0_100[k]
                   - f_12 * ig1_81[k]
                   + pa_x[k] * kg_87[k];

        t_103[k] = f_11 * ig0_101[k]
                   - f_12 * ig1_82[k]
                   + pa_x[k] * kg_88[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pa_y, pa_z, ig0_53, ig0_102, \
                         ig1_42, ig1_83, kf_42, kg_64, kg_69, kg_72, \
                         kg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_11 * ig0_102[k]
                   - f_12 * ig1_83[k]
                   + pa_x[k] * kg_89[k];

        t_105[k] = f_5 * kf_42[k]
                   + pa_y[k] * kg_69[k];

        t_106[k] = pa_y[k] * kg_72[k];

        t_107[k] = f_14 * ig0_53[k]
                   - f_15 * ig1_42[k]
                   + pa_z[k] * kg_64[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_x, pb_y, kf_53, kf_56, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_y[k] * lf_48[k];

        t_109[k] = f_3 * ld0_24[k]
                   - f_4 * ld1_24[k]
                   + pb_y[k] * lf_49[k];

        t_110[k] = f_18 * kf_53[k]
                   + f_3 * ld0_26[k]
                   - f_4 * ld1_26[k]
                   + pb_x[k] * lf_50[k];

        t_111[k] = f_18 * kf_56[k]
                   + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_x, pb_y, ig0_106, ig1_87, kg_99, \
                         ld0_25, ld0_26, ld1_25, ld1_26, lf_51, lf_52, \
                         lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * ld0_25[k]
                   - f_2 * ld1_25[k]
                   + pb_y[k] * lf_51[k];

        t_113[k] = f_3 * ld0_26[k]
                   - f_4 * ld1_26[k]
                   + pb_y[k] * lf_52[k];

        t_114[k] = pb_y[k] * lf_53[k];

        t_115[k] = f_11 * ig0_106[k]
                   - f_12 * ig1_87[k]
                   + pa_x[k] * kg_99[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pa_y, pb_x, pb_z, ig0_62, ig1_51, kf_57, kg_73, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_9 * ig0_62[k]
                   - f_10 * ig1_51[k]
                   + pa_y[k] * kg_73[k];

        t_117[k] = pb_z[k] * lf_54[k];

        t_118[k] = f_19 * kf_57[k]
                   + f_3 * ld0_28[k]
                   - f_4 * ld1_28[k]
                   + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_x, pb_x, pb_z, ig0_113, ig1_94, kf_58, \
                         kg_101, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * ld0_27[k]
                   - f_4 * ld1_27[k]
                   + pb_z[k] * lf_55[k];

        t_120[k] = f_19 * kf_58[k]
                   + pb_x[k] * lf_57[k];

        t_121[k] = f_6 * ig0_113[k]
                   - f_7 * ig1_94[k]
                   + pa_x[k] * kg_101[k];

        t_122[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_z, kg_73, kg_78, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_3 * ld0_28[k]
                   - f_4 * ld1_28[k]
                   + pb_z[k] * lf_58[k];

        t_124[k] = f_1 * ld0_29[k]
                   - f_2 * ld1_29[k]
                   + pb_z[k] * lf_59[k];

        t_125[k] = pa_z[k] * kg_73[k];

        t_126[k] = pa_z[k] * kg_78[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_x, pa_y, pa_z, ig0_74, ig0_127, ig1_60, \
                         ig1_105, kf_50, kg_81, kg_82, kg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_5 * kf_50[k]
                   + pa_z[k] * kg_81[k];

        t_128[k] = f_16 * ig0_74[k]
                   - f_17 * ig1_60[k]
                   + pa_y[k] * kg_82[k];

        t_129[k] = f_6 * ig0_127[k]
                   - f_7 * ig1_105[k]
                   + pa_x[k] * kg_102[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_x, pa_y, ig0_78, ig0_128, ig0_130, ig1_64, \
                         ig1_106, ig1_108, kg_86, kg_103, kg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_6 * ig0_128[k]
                   - f_7 * ig1_106[k]
                   + pa_x[k] * kg_103[k];

        t_131[k] = f_6 * ig0_130[k]
                   - f_7 * ig1_108[k]
                   + pa_x[k] * kg_104[k];

        t_132[k] = f_11 * ig0_78[k]
                   - f_12 * ig1_64[k]
                   + pa_y[k] * kg_86[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_x, ig0_136, ig0_137, ig0_139, ig1_114, \
                         ig1_115, ig1_117, kg_105, kg_106, kg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_6 * ig0_136[k]
                   - f_7 * ig1_114[k]
                   + pa_x[k] * kg_105[k];

        t_134[k] = f_6 * ig0_137[k]
                   - f_7 * ig1_115[k]
                   + pa_x[k] * kg_106[k];

        t_135[k] = f_6 * ig0_139[k]
                   - f_7 * ig1_117[k]
                   + pa_x[k] * kg_107[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_x, pa_y, ig0_81, ig0_145, ig0_146, ig1_65, \
                         ig1_123, ig1_124, kg_90, kg_108, kg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_6 * ig0_81[k]
                   - f_7 * ig1_65[k]
                   + pa_y[k] * kg_90[k];

        t_137[k] = f_6 * ig0_145[k]
                   - f_7 * ig1_123[k]
                   + pa_x[k] * kg_108[k];

        t_138[k] = f_6 * ig0_146[k]
                   - f_7 * ig1_124[k]
                   + pa_x[k] * kg_109[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pa_y, pa_z, ig0_81, ig0_148, \
                         ig1_65, ig1_126, kf_54, kg_91, kg_96, kg_99, \
                         kg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_6 * ig0_148[k]
                   - f_7 * ig1_126[k]
                   + pa_x[k] * kg_110[k];

        t_140[k] = f_5 * kf_54[k]
                   + pa_y[k] * kg_96[k];

        t_141[k] = pa_y[k] * kg_99[k];

        t_142[k] = f_9 * ig0_81[k]
                   - f_10 * ig1_65[k]
                   + pa_z[k] * kg_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pb_x, pb_y, kf_59, kf_60, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = pb_y[k] * lf_60[k];

        t_144[k] = f_3 * ld0_30[k]
                   - f_4 * ld1_30[k]
                   + pb_y[k] * lf_61[k];

        t_145[k] = f_19 * kf_59[k]
                   + f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_x[k] * lf_62[k];

        t_146[k] = f_19 * kf_60[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_x, pb_y, ig0_164, ig1_140, kg_112, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_148[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_149[k] = pb_y[k] * lf_65[k];

        t_150[k] = f_6 * ig0_164[k]
                   - f_7 * ig1_140[k]
                   + pa_x[k] * kg_112[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pa_x, pa_z, kf_61, kf_69, kf_75, \
                         kg_100, kg_113, kg_118, kg_124, kg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_5 * kf_61[k]
                   + pa_x[k] * kg_113[k];

        t_152[k] = pa_x[k] * kg_118[k];

        t_153[k] = pa_z[k] * kg_100[k];

        t_154[k] = f_5 * kf_69[k]
                   + pa_x[k] * kg_124[k];

        t_155[k] = f_5 * kf_75[k]
                   + pa_x[k] * kg_133[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pa_x, kf_81, kf_87, kf_96, kg_142, \
                         kg_151, kg_162, kg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_5 * kf_81[k]
                   + pa_x[k] * kg_142[k];

        t_157[k] = f_5 * kf_87[k]
                   + pa_x[k] * kg_151[k];

        t_158[k] = f_5 * kf_96[k]
                   + pa_x[k] * kg_162[k];

        t_159[k] = pa_x[k] * kg_170[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_161[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_162[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_163[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pb_x, pb_y, pb_z, kf_64, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pb_x[k] * lf_71[k];

        t_165[k] = f_0 * kf_64[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_166[k] = pb_z[k] * lf_69[k];

        t_167[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_168[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_z, pb_x, kf_66, kg_113, kg_118, \
                         kg_121, ld0_36, ld1_36, lf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pa_z[k] * kg_113[k];

        t_170[k] = pa_z[k] * kg_118[k];

        t_171[k] = f_5 * kf_66[k]
                   + pa_z[k] * kg_121[k];

        t_172[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pb_x, ld0_37, ld0_38, ld1_37, ld1_38, \
                         lf_73, lf_74, lf_75, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_174[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_175[k] = pb_x[k] * lf_75[k];

        t_176[k] = pb_x[k] * lf_77[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_z, pb_y, ig0_113, ig1_94, kf_73, kf_74, \
                         kg_122, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * ig0_113[k]
                   - f_7 * ig1_94[k]
                   + pa_z[k] * kg_122[k];

        t_178[k] = f_8 * kf_73[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_179[k] = f_8 * kf_74[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_y, pb_x, ig0_130, ig1_108, kg_132, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_9 * ig0_130[k]
                   - f_10 * ig1_108[k]
                   + pa_y[k] * kg_132[k];

        t_181[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_182[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pa_z, pb_x, ig0_120, ig1_98, kg_129, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_184[k] = pb_x[k] * lf_81[k];

        t_185[k] = pb_x[k] * lf_83[k];

        t_186[k] = f_11 * ig0_120[k]
                   - f_12 * ig1_98[k]
                   + pa_z[k] * kg_129[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pb_y, ig0_139, ig1_117, kf_79, kf_80, \
                         kg_141, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_13 * kf_79[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_188[k] = f_13 * kf_80[k]
                   + pb_y[k] * lf_83[k];

        t_189[k] = f_14 * ig0_139[k]
                   - f_15 * ig1_117[k]
                   + pa_y[k] * kg_141[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_191[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_192[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_193[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pa_z, pb_x, pb_y, ig0_127, ig1_105, \
                         kf_85, kf_86, kg_138, ld0_44, ld1_44, lf_88, \
                         lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pb_x[k] * lf_89[k];

        t_195[k] = f_16 * ig0_127[k]
                   - f_17 * ig1_105[k]
                   + pa_z[k] * kg_138[k];

        t_196[k] = f_5 * kf_85[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_197[k] = f_5 * kf_86[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pb_x, ig0_148, ig1_126, kg_150, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * ig0_148[k]
                   - f_17 * ig1_126[k]
                   + pa_y[k] * kg_150[k];

        t_199[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_200[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pa_z, pb_x, ig0_136, ig1_114, kg_147, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_202[k] = pb_x[k] * lf_93[k];

        t_203[k] = pb_x[k] * lf_95[k];

        t_204[k] = f_14 * ig0_136[k]
                   - f_15 * ig1_114[k]
                   + pa_z[k] * kg_147[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pa_y, pb_y, ig0_154, ig1_130, kf_91, kf_92, \
                         kg_159, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_18 * kf_91[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_206[k] = f_18 * kf_92[k]
                   + pb_y[k] * lf_95[k];

        t_207[k] = f_11 * ig0_154[k]
                   - f_12 * ig1_130[k]
                   + pa_y[k] * kg_159[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_209[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_210[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_211[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_z, pb_x, pb_y, ig0_145, ig1_123, \
                         kf_94, kf_95, kg_156, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_x[k] * lf_101[k];

        t_213[k] = f_9 * ig0_145[k]
                   - f_10 * ig1_123[k]
                   + pa_z[k] * kg_156[k];

        t_214[k] = f_19 * kf_94[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_215[k] = f_19 * kf_95[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_y, pb_x, ig0_164, ig1_140, kf_99, \
                         kg_161, kg_167, kg_170, ld0_51, ld1_51, \
                         lf_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_6 * ig0_164[k]
                   - f_7 * ig1_140[k]
                   + pa_y[k] * kg_161[k];

        t_217[k] = f_5 * kf_99[k]
                   + pa_y[k] * kg_167[k];

        t_218[k] = pa_y[k] * kg_170[k];

        t_219[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_103, lf_104, lf_105, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];

        t_221[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_222[k] = pb_x[k] * lf_105[k];

        t_223[k] = pb_x[k] * lf_107[k];

        t_224[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pb_y, pb_z, kf_101, ld0_53, ld1_53, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_226[k] = pb_y[k] * lf_107[k];

        t_227[k] = f_0 * kf_101[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_1 = buffer.data(ig0 + 1);
    const auto *ig0_2 = buffer.data(ig0 + 2);
    const auto *ig0_3 = buffer.data(ig0 + 3);
    const auto *ig0_6 = buffer.data(ig0 + 6);
    const auto *ig0_7 = buffer.data(ig0 + 7);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_14 = buffer.data(ig0 + 14);
    const auto *ig0_15 = buffer.data(ig0 + 15);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_19 = buffer.data(ig0 + 19);
    const auto *ig0_22 = buffer.data(ig0 + 22);
    const auto *ig0_23 = buffer.data(ig0 + 23);
    const auto *ig0_24 = buffer.data(ig0 + 24);
    const auto *ig0_27 = buffer.data(ig0 + 27);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_31 = buffer.data(ig0 + 31);
    const auto *ig0_32 = buffer.data(ig0 + 32);
    const auto *ig0_33 = buffer.data(ig0 + 33);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_35 = buffer.data(ig0 + 35);
    const auto *ig0_37 = buffer.data(ig0 + 37);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_39 = buffer.data(ig0 + 39);
    const auto *ig0_41 = buffer.data(ig0 + 41);
    const auto *ig0_42 = buffer.data(ig0 + 42);
    const auto *ig0_43 = buffer.data(ig0 + 43);
    const auto *ig0_45 = buffer.data(ig0 + 45);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_9 = buffer.data(ig1 + 9);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_16 = buffer.data(ig1 + 16);
    const auto *ig1_20 = buffer.data(ig1 + 20);
    const auto *ig1_28 = buffer.data(ig1 + 28);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_34 = buffer.data(ig1 + 34);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_46 = buffer.data(ig1 + 46);
    const auto *ig1_47 = buffer.data(ig1 + 47);
    const auto *ig1_52 = buffer.data(ig1 + 52);
    const auto *ig1_56 = buffer.data(ig1 + 56);
    const auto *ig1_57 = buffer.data(ig1 + 57);
    const auto *ig1_65 = buffer.data(ig1 + 65);
    const auto *ig1_68 = buffer.data(ig1 + 68);
    const auto *ig1_69 = buffer.data(ig1 + 69);
    const auto *ig1_70 = buffer.data(ig1 + 70);
    const auto *ig1_73 = buffer.data(ig1 + 73);
    const auto *ig1_79 = buffer.data(ig1 + 79);
    const auto *ig1_83 = buffer.data(ig1 + 83);
    const auto *ig1_89 = buffer.data(ig1 + 89);
    const auto *ig1_90 = buffer.data(ig1 + 90);
    const auto *ig1_92 = buffer.data(ig1 + 92);
    const auto *ig1_98 = buffer.data(ig1 + 98);
    const auto *ig1_99 = buffer.data(ig1 + 99);
    const auto *ig1_101 = buffer.data(ig1 + 101);
    const auto *ig1_107 = buffer.data(ig1 + 107);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_110 = buffer.data(ig1 + 110);
    const auto *ig1_113 = buffer.data(ig1 + 113);
    const auto *ig1_122 = buffer.data(ig1 + 122);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_98 = buffer.data(kf + 98);

    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_143 = buffer.data(kg + 143);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_10, kg_9, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_9[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_10[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_6, ig1_16, kf_11, \
                         kg_16, ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_11[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_6[k]
                  - f_9 * ig1_16[k]
                  + pa_x[k] * kg_16[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_10, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_16, kf_19, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_16[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_19[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_10, ig1_28, kg_28, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_10[k]
                  - f_9 * ig1_28[k]
                  + pa_x[k] * kg_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_1, ig1_9, kf_22, kg_11, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_1[k]
                  - f_11 * ig1_9[k]
                  + pa_y[k] * kg_11[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_22[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_14, ig1_34, kf_23, \
                         kg_34, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_23[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_14[k]
                  - f_14 * ig1_34[k]
                  + pa_x[k] * kg_34[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_2, ig1_10, kg_20, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_2[k]
                  - f_11 * ig1_10[k]
                  + pa_z[k] * kg_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_28, kf_31, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_28[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_18, ig1_46, kg_46, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_18[k]
                  - f_14 * ig1_46[k]
                  + pa_x[k] * kg_46[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_3, ig1_11, kf_34, kg_29, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_3[k]
                  - f_16 * ig1_11[k]
                  + pa_y[k] * kg_29[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_34[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_22, ig1_52, kf_35, \
                         kg_52, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_35[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_22[k]
                  - f_16 * ig1_52[k]
                  + pa_x[k] * kg_52[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_23, ig1_56, kg_56, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_23[k]
                  - f_16 * ig1_56[k]
                  + pa_x[k] * kg_56[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_7, ig1_20, kg_38, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_7[k]
                  - f_16 * ig1_20[k]
                  + pa_z[k] * kg_38[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_40, kf_43, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_40[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_11, ig0_27, \
                         ig1_29, ig1_65, kg_47, kg_65, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_27[k]
                  - f_16 * ig1_65[k]
                  + pa_x[k] * kg_65[k];

        t_64[k] = f_13 * ig0_11[k]
                  - f_14 * ig1_29[k]
                  + pa_y[k] * kg_47[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_46, kf_47, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_46[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_47[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_28, ig1_68, kg_71, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_28[k]
                  - f_11 * ig1_68[k]
                  + pa_x[k] * kg_71[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_15, ig0_29, ig0_30, ig1_38, ig1_69, \
                         ig1_70, kg_57, kg_75, kg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_29[k]
                  - f_11 * ig1_69[k]
                  + pa_x[k] * kg_75[k];

        t_74[k] = f_10 * ig0_30[k]
                  - f_11 * ig1_70[k]
                  + pa_x[k] * kg_76[k];

        t_75[k] = f_13 * ig0_15[k]
                  - f_14 * ig1_38[k]
                  + pa_z[k] * kg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_52, kf_55, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_52[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_55[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_31, ig1_73, kg_85, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_31[k]
                  - f_11 * ig1_73[k]
                  + pa_x[k] * kg_85[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_19, ig1_47, kf_56, kg_66, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_19[k]
                  - f_9 * ig1_47[k]
                  + pa_y[k] * kg_66[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_56[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_32, ig1_79, kf_57, \
                         kg_88, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_32[k]
                  - f_6 * ig1_79[k]
                  + pa_x[k] * kg_88[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_35, ig1_90, kg_89, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_35[k]
                  - f_6 * ig1_90[k]
                  + pa_x[k] * kg_89[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_24, ig0_39, ig0_43, ig1_57, ig1_99, \
                         ig1_108, kg_77, kg_90, kg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_39[k]
                  - f_6 * ig1_99[k]
                  + pa_x[k] * kg_90[k];

        t_95[k] = f_5 * ig0_43[k]
                  - f_6 * ig1_108[k]
                  + pa_x[k] * kg_91[k];

        t_96[k] = f_8 * ig0_24[k]
                  - f_9 * ig1_57[k]
                  + pa_z[k] * kg_77[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_58, kf_59, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_58[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_59[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_47, ig1_122, kg_94, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_47[k]
                   - f_6 * ig1_122[k]
                   + pa_x[k] * kg_94[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_63, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_63[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_32, ig1_79, kf_71, \
                         kf_72, kg_104, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_32[k]
                   - f_6 * ig1_79[k]
                   + pa_z[k] * kg_104[k];

        t_120[k] = f_7 * kf_71[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_72[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_37, ig1_92, kg_113, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_37[k]
                   - f_9 * ig1_92[k]
                   + pa_y[k] * kg_113[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_33, ig1_83, kg_110, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_33[k]
                   - f_11 * ig1_83[k]
                   + pa_z[k] * kg_110[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_41, ig1_101, kf_77, kf_78, \
                         kg_122, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_77[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_78[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_41[k]
                   - f_14 * ig1_101[k]
                   + pa_y[k] * kg_122[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_34, ig1_89, kf_83, \
                         kf_84, kg_119, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_34[k]
                   - f_16 * ig1_89[k]
                   + pa_z[k] * kg_119[k];

        t_138[k] = f_17 * kf_83[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_84[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_45, ig1_110, kg_131, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_45[k]
                   - f_16 * ig1_110[k]
                   + pa_y[k] * kg_131[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_38, ig1_98, kg_128, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_38[k]
                   - f_14 * ig1_98[k]
                   + pa_z[k] * kg_128[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_46, ig1_113, kf_89, kf_90, \
                         kg_140, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_89[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_90[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_46[k]
                   - f_11 * ig1_113[k]
                   + pa_y[k] * kg_140[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_42, ig1_107, kf_91, \
                         kf_92, kg_137, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_42[k]
                   - f_9 * ig1_107[k]
                   + pa_z[k] * kg_137[k];

        t_156[k] = f_19 * kf_91[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_92[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_47, ig1_122, kg_143, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_47[k]
                   - f_6 * ig1_122[k]
                   + pa_y[k] * kg_143[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_98, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_98[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_9 = buffer.data(ig0 + 9);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_16 = buffer.data(ig0 + 16);
    const auto *ig0_20 = buffer.data(ig0 + 20);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);
    const auto *ig0_52 = buffer.data(ig0 + 52);
    const auto *ig0_56 = buffer.data(ig0 + 56);
    const auto *ig0_57 = buffer.data(ig0 + 57);
    const auto *ig0_65 = buffer.data(ig0 + 65);
    const auto *ig0_68 = buffer.data(ig0 + 68);
    const auto *ig0_69 = buffer.data(ig0 + 69);
    const auto *ig0_70 = buffer.data(ig0 + 70);
    const auto *ig0_73 = buffer.data(ig0 + 73);
    const auto *ig0_79 = buffer.data(ig0 + 79);
    const auto *ig0_83 = buffer.data(ig0 + 83);
    const auto *ig0_89 = buffer.data(ig0 + 89);
    const auto *ig0_90 = buffer.data(ig0 + 90);
    const auto *ig0_92 = buffer.data(ig0 + 92);
    const auto *ig0_98 = buffer.data(ig0 + 98);
    const auto *ig0_99 = buffer.data(ig0 + 99);
    const auto *ig0_101 = buffer.data(ig0 + 101);
    const auto *ig0_107 = buffer.data(ig0 + 107);
    const auto *ig0_108 = buffer.data(ig0 + 108);
    const auto *ig0_110 = buffer.data(ig0 + 110);
    const auto *ig0_113 = buffer.data(ig0 + 113);
    const auto *ig0_122 = buffer.data(ig0 + 122);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_13 = buffer.data(ig1 + 13);
    const auto *ig1_18 = buffer.data(ig1 + 18);
    const auto *ig1_22 = buffer.data(ig1 + 22);
    const auto *ig1_30 = buffer.data(ig1 + 30);
    const auto *ig1_31 = buffer.data(ig1 + 31);
    const auto *ig1_36 = buffer.data(ig1 + 36);
    const auto *ig1_40 = buffer.data(ig1 + 40);
    const auto *ig1_48 = buffer.data(ig1 + 48);
    const auto *ig1_49 = buffer.data(ig1 + 49);
    const auto *ig1_54 = buffer.data(ig1 + 54);
    const auto *ig1_58 = buffer.data(ig1 + 58);
    const auto *ig1_59 = buffer.data(ig1 + 59);
    const auto *ig1_67 = buffer.data(ig1 + 67);
    const auto *ig1_70 = buffer.data(ig1 + 70);
    const auto *ig1_71 = buffer.data(ig1 + 71);
    const auto *ig1_72 = buffer.data(ig1 + 72);
    const auto *ig1_75 = buffer.data(ig1 + 75);
    const auto *ig1_82 = buffer.data(ig1 + 82);
    const auto *ig1_86 = buffer.data(ig1 + 86);
    const auto *ig1_93 = buffer.data(ig1 + 93);
    const auto *ig1_94 = buffer.data(ig1 + 94);
    const auto *ig1_96 = buffer.data(ig1 + 96);
    const auto *ig1_102 = buffer.data(ig1 + 102);
    const auto *ig1_103 = buffer.data(ig1 + 103);
    const auto *ig1_105 = buffer.data(ig1 + 105);
    const auto *ig1_111 = buffer.data(ig1 + 111);
    const auto *ig1_112 = buffer.data(ig1 + 112);
    const auto *ig1_114 = buffer.data(ig1 + 114);
    const auto *ig1_118 = buffer.data(ig1 + 118);
    const auto *ig1_128 = buffer.data(ig1 + 128);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_98 = buffer.data(kf + 98);

    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_92 = buffer.data(kg + 92);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_146 = buffer.data(kg + 146);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_10, kg_9, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_9[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_10[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_16, ig1_18, kf_11, \
                         kg_17, ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_11[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_16[k]
                  - f_9 * ig1_18[k]
                  + pa_x[k] * kg_17[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_10, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_16, kf_19, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_16[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_19[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_28, ig1_30, kg_29, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_28[k]
                  - f_9 * ig1_30[k]
                  + pa_x[k] * kg_29[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_9, ig1_10, kf_22, kg_12, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_9[k]
                  - f_11 * ig1_10[k]
                  + pa_y[k] * kg_12[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_22[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_34, ig1_36, kf_23, \
                         kg_35, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_23[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_34[k]
                  - f_14 * ig1_36[k]
                  + pa_x[k] * kg_35[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_10, ig1_11, kg_21, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_10[k]
                  - f_11 * ig1_11[k]
                  + pa_z[k] * kg_21[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_28, kf_31, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_28[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_46, ig1_48, kg_47, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_46[k]
                  - f_14 * ig1_48[k]
                  + pa_x[k] * kg_47[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_11, ig1_13, kf_34, kg_30, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_11[k]
                  - f_16 * ig1_13[k]
                  + pa_y[k] * kg_30[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_34[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_52, ig1_54, kf_35, \
                         kg_53, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_35[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_52[k]
                  - f_16 * ig1_54[k]
                  + pa_x[k] * kg_53[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_56, ig1_58, kg_57, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_56[k]
                  - f_16 * ig1_58[k]
                  + pa_x[k] * kg_57[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_20, ig1_22, kg_39, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_20[k]
                  - f_16 * ig1_22[k]
                  + pa_z[k] * kg_39[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_40, kf_43, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_40[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_29, ig0_65, \
                         ig1_31, ig1_67, kg_48, kg_66, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_65[k]
                  - f_16 * ig1_67[k]
                  + pa_x[k] * kg_66[k];

        t_64[k] = f_13 * ig0_29[k]
                  - f_14 * ig1_31[k]
                  + pa_y[k] * kg_48[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_46, kf_47, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_46[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_47[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_68, ig1_70, kg_72, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_68[k]
                  - f_11 * ig1_70[k]
                  + pa_x[k] * kg_72[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_38, ig0_69, ig0_70, ig1_40, ig1_71, \
                         ig1_72, kg_58, kg_76, kg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_69[k]
                  - f_11 * ig1_71[k]
                  + pa_x[k] * kg_76[k];

        t_74[k] = f_10 * ig0_70[k]
                  - f_11 * ig1_72[k]
                  + pa_x[k] * kg_77[k];

        t_75[k] = f_13 * ig0_38[k]
                  - f_14 * ig1_40[k]
                  + pa_z[k] * kg_58[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_52, kf_55, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_52[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_55[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_73, ig1_75, kg_86, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_73[k]
                  - f_11 * ig1_75[k]
                  + pa_x[k] * kg_86[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_47, ig1_49, kf_56, kg_67, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_47[k]
                  - f_9 * ig1_49[k]
                  + pa_y[k] * kg_67[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_56[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_79, ig1_82, kf_57, \
                         kg_89, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_79[k]
                  - f_6 * ig1_82[k]
                  + pa_x[k] * kg_89[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_90, ig1_94, kg_90, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_90[k]
                  - f_6 * ig1_94[k]
                  + pa_x[k] * kg_90[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_57, ig0_99, ig0_108, ig1_59, \
                         ig1_103, ig1_112, kg_78, kg_91, kg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_99[k]
                  - f_6 * ig1_103[k]
                  + pa_x[k] * kg_91[k];

        t_95[k] = f_5 * ig0_108[k]
                  - f_6 * ig1_112[k]
                  + pa_x[k] * kg_92[k];

        t_96[k] = f_8 * ig0_57[k]
                  - f_9 * ig1_59[k]
                  + pa_z[k] * kg_78[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_58, kf_59, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_58[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_59[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_122, ig1_128, kg_95, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_122[k]
                   - f_6 * ig1_128[k]
                   + pa_x[k] * kg_95[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_63, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_63[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_79, ig1_82, kf_71, \
                         kf_72, kg_105, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_79[k]
                   - f_6 * ig1_82[k]
                   + pa_z[k] * kg_105[k];

        t_120[k] = f_7 * kf_71[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_72[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_92, ig1_96, kg_115, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_92[k]
                   - f_9 * ig1_96[k]
                   + pa_y[k] * kg_115[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_83, ig1_86, kg_112, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_83[k]
                   - f_11 * ig1_86[k]
                   + pa_z[k] * kg_112[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_101, ig1_105, kf_77, kf_78, \
                         kg_124, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_77[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_78[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_101[k]
                   - f_14 * ig1_105[k]
                   + pa_y[k] * kg_124[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_89, ig1_93, kf_83, \
                         kf_84, kg_121, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_89[k]
                   - f_16 * ig1_93[k]
                   + pa_z[k] * kg_121[k];

        t_138[k] = f_17 * kf_83[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_84[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_110, ig1_114, kg_133, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_110[k]
                   - f_16 * ig1_114[k]
                   + pa_y[k] * kg_133[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_98, ig1_102, kg_130, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_98[k]
                   - f_14 * ig1_102[k]
                   + pa_z[k] * kg_130[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_113, ig1_118, kf_89, kf_90, \
                         kg_142, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_89[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_90[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_113[k]
                   - f_11 * ig1_118[k]
                   + pa_y[k] * kg_142[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_107, ig1_111, \
                         kf_91, kf_92, kg_139, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_107[k]
                   - f_9 * ig1_111[k]
                   + pa_z[k] * kg_139[k];

        t_156[k] = f_19 * kf_91[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_92[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_122, ig1_128, kg_146, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_122[k]
                   - f_6 * ig1_128[k]
                   + pa_y[k] * kg_146[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_98, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_98[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_13 = buffer.data(ig0 + 13);
    const auto *ig0_18 = buffer.data(ig0 + 18);
    const auto *ig0_22 = buffer.data(ig0 + 22);
    const auto *ig0_30 = buffer.data(ig0 + 30);
    const auto *ig0_31 = buffer.data(ig0 + 31);
    const auto *ig0_36 = buffer.data(ig0 + 36);
    const auto *ig0_40 = buffer.data(ig0 + 40);
    const auto *ig0_48 = buffer.data(ig0 + 48);
    const auto *ig0_49 = buffer.data(ig0 + 49);
    const auto *ig0_54 = buffer.data(ig0 + 54);
    const auto *ig0_58 = buffer.data(ig0 + 58);
    const auto *ig0_59 = buffer.data(ig0 + 59);
    const auto *ig0_67 = buffer.data(ig0 + 67);
    const auto *ig0_70 = buffer.data(ig0 + 70);
    const auto *ig0_71 = buffer.data(ig0 + 71);
    const auto *ig0_72 = buffer.data(ig0 + 72);
    const auto *ig0_75 = buffer.data(ig0 + 75);
    const auto *ig0_82 = buffer.data(ig0 + 82);
    const auto *ig0_86 = buffer.data(ig0 + 86);
    const auto *ig0_93 = buffer.data(ig0 + 93);
    const auto *ig0_94 = buffer.data(ig0 + 94);
    const auto *ig0_96 = buffer.data(ig0 + 96);
    const auto *ig0_102 = buffer.data(ig0 + 102);
    const auto *ig0_103 = buffer.data(ig0 + 103);
    const auto *ig0_105 = buffer.data(ig0 + 105);
    const auto *ig0_111 = buffer.data(ig0 + 111);
    const auto *ig0_112 = buffer.data(ig0 + 112);
    const auto *ig0_114 = buffer.data(ig0 + 114);
    const auto *ig0_118 = buffer.data(ig0 + 118);
    const auto *ig0_128 = buffer.data(ig0 + 128);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_9 = buffer.data(ig1 + 9);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_16 = buffer.data(ig1 + 16);
    const auto *ig1_20 = buffer.data(ig1 + 20);
    const auto *ig1_28 = buffer.data(ig1 + 28);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_34 = buffer.data(ig1 + 34);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_46 = buffer.data(ig1 + 46);
    const auto *ig1_47 = buffer.data(ig1 + 47);
    const auto *ig1_52 = buffer.data(ig1 + 52);
    const auto *ig1_56 = buffer.data(ig1 + 56);
    const auto *ig1_57 = buffer.data(ig1 + 57);
    const auto *ig1_65 = buffer.data(ig1 + 65);
    const auto *ig1_68 = buffer.data(ig1 + 68);
    const auto *ig1_69 = buffer.data(ig1 + 69);
    const auto *ig1_70 = buffer.data(ig1 + 70);
    const auto *ig1_73 = buffer.data(ig1 + 73);
    const auto *ig1_79 = buffer.data(ig1 + 79);
    const auto *ig1_83 = buffer.data(ig1 + 83);
    const auto *ig1_89 = buffer.data(ig1 + 89);
    const auto *ig1_90 = buffer.data(ig1 + 90);
    const auto *ig1_92 = buffer.data(ig1 + 92);
    const auto *ig1_98 = buffer.data(ig1 + 98);
    const auto *ig1_99 = buffer.data(ig1 + 99);
    const auto *ig1_101 = buffer.data(ig1 + 101);
    const auto *ig1_107 = buffer.data(ig1 + 107);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_110 = buffer.data(ig1 + 110);
    const auto *ig1_113 = buffer.data(ig1 + 113);
    const auto *ig1_122 = buffer.data(ig1 + 122);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_98 = buffer.data(kf + 98);

    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_10, kg_9, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_9[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_10[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_18, ig1_16, kf_11, \
                         kg_16, ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_11[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_18[k]
                  - f_9 * ig1_16[k]
                  + pa_x[k] * kg_16[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_10, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_16, kf_19, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_16[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_19[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_30, ig1_28, kg_28, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_30[k]
                  - f_9 * ig1_28[k]
                  + pa_x[k] * kg_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_10, ig1_9, kf_22, kg_11, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_10[k]
                  - f_11 * ig1_9[k]
                  + pa_y[k] * kg_11[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_22[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_36, ig1_34, kf_23, \
                         kg_34, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_23[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_36[k]
                  - f_14 * ig1_34[k]
                  + pa_x[k] * kg_34[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_11, ig1_10, kg_20, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_11[k]
                  - f_11 * ig1_10[k]
                  + pa_z[k] * kg_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_28, kf_31, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_28[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_48, ig1_46, kg_46, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_48[k]
                  - f_14 * ig1_46[k]
                  + pa_x[k] * kg_46[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_13, ig1_11, kf_34, kg_29, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_13[k]
                  - f_16 * ig1_11[k]
                  + pa_y[k] * kg_29[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_34[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_54, ig1_52, kf_35, \
                         kg_52, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_35[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_54[k]
                  - f_16 * ig1_52[k]
                  + pa_x[k] * kg_52[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_58, ig1_56, kg_56, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_58[k]
                  - f_16 * ig1_56[k]
                  + pa_x[k] * kg_56[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_22, ig1_20, kg_38, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_22[k]
                  - f_16 * ig1_20[k]
                  + pa_z[k] * kg_38[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_40, kf_43, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_40[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_31, ig0_67, \
                         ig1_29, ig1_65, kg_47, kg_65, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_67[k]
                  - f_16 * ig1_65[k]
                  + pa_x[k] * kg_65[k];

        t_64[k] = f_13 * ig0_31[k]
                  - f_14 * ig1_29[k]
                  + pa_y[k] * kg_47[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_46, kf_47, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_46[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_47[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_70, ig1_68, kg_71, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_70[k]
                  - f_11 * ig1_68[k]
                  + pa_x[k] * kg_71[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_40, ig0_71, ig0_72, ig1_38, ig1_69, \
                         ig1_70, kg_57, kg_75, kg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_71[k]
                  - f_11 * ig1_69[k]
                  + pa_x[k] * kg_75[k];

        t_74[k] = f_10 * ig0_72[k]
                  - f_11 * ig1_70[k]
                  + pa_x[k] * kg_76[k];

        t_75[k] = f_13 * ig0_40[k]
                  - f_14 * ig1_38[k]
                  + pa_z[k] * kg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_52, kf_55, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_52[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_55[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_75, ig1_73, kg_85, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_75[k]
                  - f_11 * ig1_73[k]
                  + pa_x[k] * kg_85[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_49, ig1_47, kf_56, kg_66, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_49[k]
                  - f_9 * ig1_47[k]
                  + pa_y[k] * kg_66[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_56[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_82, ig1_79, kf_57, \
                         kg_86, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_82[k]
                  - f_6 * ig1_79[k]
                  + pa_x[k] * kg_86[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_94, ig1_90, kg_87, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_94[k]
                  - f_6 * ig1_90[k]
                  + pa_x[k] * kg_87[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_59, ig0_103, ig0_112, ig1_57, \
                         ig1_99, ig1_108, kg_77, kg_88, kg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_103[k]
                  - f_6 * ig1_99[k]
                  + pa_x[k] * kg_88[k];

        t_95[k] = f_5 * ig0_112[k]
                  - f_6 * ig1_108[k]
                  + pa_x[k] * kg_89[k];

        t_96[k] = f_8 * ig0_59[k]
                  - f_9 * ig1_57[k]
                  + pa_z[k] * kg_77[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_58, kf_59, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_58[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_59[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_128, ig1_122, kg_90, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_128[k]
                   - f_6 * ig1_122[k]
                   + pa_x[k] * kg_90[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_63, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_63[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_82, ig1_79, kf_71, \
                         kf_72, kg_100, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_82[k]
                   - f_6 * ig1_79[k]
                   + pa_z[k] * kg_100[k];

        t_120[k] = f_7 * kf_71[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_72[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_96, ig1_92, kg_109, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_96[k]
                   - f_9 * ig1_92[k]
                   + pa_y[k] * kg_109[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_86, ig1_83, kg_106, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_86[k]
                   - f_11 * ig1_83[k]
                   + pa_z[k] * kg_106[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_105, ig1_101, kf_77, kf_78, \
                         kg_118, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_77[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_78[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_105[k]
                   - f_14 * ig1_101[k]
                   + pa_y[k] * kg_118[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_93, ig1_89, kf_83, \
                         kf_84, kg_115, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_93[k]
                   - f_16 * ig1_89[k]
                   + pa_z[k] * kg_115[k];

        t_138[k] = f_17 * kf_83[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_84[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_114, ig1_110, kg_127, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_114[k]
                   - f_16 * ig1_110[k]
                   + pa_y[k] * kg_127[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_102, ig1_98, kg_124, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_102[k]
                   - f_14 * ig1_98[k]
                   + pa_z[k] * kg_124[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_118, ig1_113, kf_89, kf_90, \
                         kg_136, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_89[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_90[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_118[k]
                   - f_11 * ig1_113[k]
                   + pa_y[k] * kg_136[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_111, ig1_107, \
                         kf_91, kf_92, kg_133, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_111[k]
                   - f_9 * ig1_107[k]
                   + pa_z[k] * kg_133[k];

        t_156[k] = f_19 * kf_91[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_92[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_128, ig1_122, kg_137, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_128[k]
                   - f_6 * ig1_122[k]
                   + pa_y[k] * kg_137[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_98, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_98[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_9 = buffer.data(ig0 + 9);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_16 = buffer.data(ig0 + 16);
    const auto *ig0_20 = buffer.data(ig0 + 20);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);
    const auto *ig0_52 = buffer.data(ig0 + 52);
    const auto *ig0_56 = buffer.data(ig0 + 56);
    const auto *ig0_57 = buffer.data(ig0 + 57);
    const auto *ig0_65 = buffer.data(ig0 + 65);
    const auto *ig0_68 = buffer.data(ig0 + 68);
    const auto *ig0_69 = buffer.data(ig0 + 69);
    const auto *ig0_70 = buffer.data(ig0 + 70);
    const auto *ig0_73 = buffer.data(ig0 + 73);
    const auto *ig0_79 = buffer.data(ig0 + 79);
    const auto *ig0_83 = buffer.data(ig0 + 83);
    const auto *ig0_89 = buffer.data(ig0 + 89);
    const auto *ig0_90 = buffer.data(ig0 + 90);
    const auto *ig0_92 = buffer.data(ig0 + 92);
    const auto *ig0_98 = buffer.data(ig0 + 98);
    const auto *ig0_99 = buffer.data(ig0 + 99);
    const auto *ig0_101 = buffer.data(ig0 + 101);
    const auto *ig0_107 = buffer.data(ig0 + 107);
    const auto *ig0_108 = buffer.data(ig0 + 108);
    const auto *ig0_110 = buffer.data(ig0 + 110);
    const auto *ig0_113 = buffer.data(ig0 + 113);
    const auto *ig0_122 = buffer.data(ig0 + 122);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_9 = buffer.data(ig1 + 9);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_16 = buffer.data(ig1 + 16);
    const auto *ig1_20 = buffer.data(ig1 + 20);
    const auto *ig1_28 = buffer.data(ig1 + 28);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_34 = buffer.data(ig1 + 34);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_46 = buffer.data(ig1 + 46);
    const auto *ig1_47 = buffer.data(ig1 + 47);
    const auto *ig1_52 = buffer.data(ig1 + 52);
    const auto *ig1_56 = buffer.data(ig1 + 56);
    const auto *ig1_57 = buffer.data(ig1 + 57);
    const auto *ig1_65 = buffer.data(ig1 + 65);
    const auto *ig1_68 = buffer.data(ig1 + 68);
    const auto *ig1_69 = buffer.data(ig1 + 69);
    const auto *ig1_70 = buffer.data(ig1 + 70);
    const auto *ig1_73 = buffer.data(ig1 + 73);
    const auto *ig1_79 = buffer.data(ig1 + 79);
    const auto *ig1_83 = buffer.data(ig1 + 83);
    const auto *ig1_89 = buffer.data(ig1 + 89);
    const auto *ig1_90 = buffer.data(ig1 + 90);
    const auto *ig1_92 = buffer.data(ig1 + 92);
    const auto *ig1_98 = buffer.data(ig1 + 98);
    const auto *ig1_99 = buffer.data(ig1 + 99);
    const auto *ig1_101 = buffer.data(ig1 + 101);
    const auto *ig1_107 = buffer.data(ig1 + 107);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_110 = buffer.data(ig1 + 110);
    const auto *ig1_113 = buffer.data(ig1 + 113);
    const auto *ig1_122 = buffer.data(ig1 + 122);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_98 = buffer.data(kf + 98);

    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_143 = buffer.data(kg + 143);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_10, kg_9, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_9[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_10[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_16, ig1_16, kf_11, \
                         kg_16, ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_11[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_16[k]
                  - f_9 * ig1_16[k]
                  + pa_x[k] * kg_16[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_10, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_16, kf_19, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_16[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_19[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_28, ig1_28, kg_28, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_28[k]
                  - f_9 * ig1_28[k]
                  + pa_x[k] * kg_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_9, ig1_9, kf_22, kg_11, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_9[k]
                  - f_11 * ig1_9[k]
                  + pa_y[k] * kg_11[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_22[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_34, ig1_34, kf_23, \
                         kg_34, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_23[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_34[k]
                  - f_14 * ig1_34[k]
                  + pa_x[k] * kg_34[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_10, ig1_10, kg_20, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_10[k]
                  - f_11 * ig1_10[k]
                  + pa_z[k] * kg_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_28, kf_31, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_28[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_46, ig1_46, kg_46, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_46[k]
                  - f_14 * ig1_46[k]
                  + pa_x[k] * kg_46[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_11, ig1_11, kf_34, kg_29, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_11[k]
                  - f_16 * ig1_11[k]
                  + pa_y[k] * kg_29[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_34[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_52, ig1_52, kf_35, \
                         kg_52, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_35[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_52[k]
                  - f_16 * ig1_52[k]
                  + pa_x[k] * kg_52[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_56, ig1_56, kg_56, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_56[k]
                  - f_16 * ig1_56[k]
                  + pa_x[k] * kg_56[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_20, ig1_20, kg_38, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_20[k]
                  - f_16 * ig1_20[k]
                  + pa_z[k] * kg_38[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_40, kf_43, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_40[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_29, ig0_65, \
                         ig1_29, ig1_65, kg_47, kg_65, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_65[k]
                  - f_16 * ig1_65[k]
                  + pa_x[k] * kg_65[k];

        t_64[k] = f_13 * ig0_29[k]
                  - f_14 * ig1_29[k]
                  + pa_y[k] * kg_47[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_46, kf_47, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_46[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_47[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_68, ig1_68, kg_71, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_68[k]
                  - f_11 * ig1_68[k]
                  + pa_x[k] * kg_71[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_38, ig0_69, ig0_70, ig1_38, ig1_69, \
                         ig1_70, kg_57, kg_75, kg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_69[k]
                  - f_11 * ig1_69[k]
                  + pa_x[k] * kg_75[k];

        t_74[k] = f_10 * ig0_70[k]
                  - f_11 * ig1_70[k]
                  + pa_x[k] * kg_76[k];

        t_75[k] = f_13 * ig0_38[k]
                  - f_14 * ig1_38[k]
                  + pa_z[k] * kg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_52, kf_55, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_52[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_55[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_73, ig1_73, kg_85, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_73[k]
                  - f_11 * ig1_73[k]
                  + pa_x[k] * kg_85[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_47, ig1_47, kf_56, kg_66, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_47[k]
                  - f_9 * ig1_47[k]
                  + pa_y[k] * kg_66[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_56[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_79, ig1_79, kf_57, \
                         kg_88, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_79[k]
                  - f_6 * ig1_79[k]
                  + pa_x[k] * kg_88[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_90, ig1_90, kg_89, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_90[k]
                  - f_6 * ig1_90[k]
                  + pa_x[k] * kg_89[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_57, ig0_99, ig0_108, ig1_57, \
                         ig1_99, ig1_108, kg_77, kg_90, kg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_99[k]
                  - f_6 * ig1_99[k]
                  + pa_x[k] * kg_90[k];

        t_95[k] = f_5 * ig0_108[k]
                  - f_6 * ig1_108[k]
                  + pa_x[k] * kg_91[k];

        t_96[k] = f_8 * ig0_57[k]
                  - f_9 * ig1_57[k]
                  + pa_z[k] * kg_77[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_58, kf_59, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_58[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_59[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_122, ig1_122, kg_94, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_122[k]
                   - f_6 * ig1_122[k]
                   + pa_x[k] * kg_94[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_63, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_63[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_79, ig1_79, kf_71, \
                         kf_72, kg_104, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_79[k]
                   - f_6 * ig1_79[k]
                   + pa_z[k] * kg_104[k];

        t_120[k] = f_7 * kf_71[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_72[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_92, ig1_92, kg_113, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_92[k]
                   - f_9 * ig1_92[k]
                   + pa_y[k] * kg_113[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_83, ig1_83, kg_110, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_83[k]
                   - f_11 * ig1_83[k]
                   + pa_z[k] * kg_110[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_101, ig1_101, kf_77, kf_78, \
                         kg_122, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_77[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_78[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_101[k]
                   - f_14 * ig1_101[k]
                   + pa_y[k] * kg_122[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_89, ig1_89, kf_83, \
                         kf_84, kg_119, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_89[k]
                   - f_16 * ig1_89[k]
                   + pa_z[k] * kg_119[k];

        t_138[k] = f_17 * kf_83[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_84[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_110, ig1_110, kg_131, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_110[k]
                   - f_16 * ig1_110[k]
                   + pa_y[k] * kg_131[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_98, ig1_98, kg_128, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_98[k]
                   - f_14 * ig1_98[k]
                   + pa_z[k] * kg_128[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_113, ig1_113, kf_89, kf_90, \
                         kg_140, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_89[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_90[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_113[k]
                   - f_11 * ig1_113[k]
                   + pa_y[k] * kg_140[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_107, ig1_107, \
                         kf_91, kf_92, kg_137, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_107[k]
                   - f_9 * ig1_107[k]
                   + pa_z[k] * kg_137[k];

        t_156[k] = f_19 * kf_91[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_92[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_122, ig1_122, kg_143, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_122[k]
                   - f_6 * ig1_122[k]
                   + pa_y[k] * kg_143[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_98, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_98[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

auto
compute_prim_lg_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
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
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / alpha;
    const auto f_14 = 2.0 * beta / (alpha * p);
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);
    const auto f_17 = 2.0 / p;
    const auto f_18 = 1.5 / p;
    const auto f_19 = 1.0 / p;

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

    const auto *ig0_0 = buffer.data(ig0 + 0);
    const auto *ig0_9 = buffer.data(ig0 + 9);
    const auto *ig0_10 = buffer.data(ig0 + 10);
    const auto *ig0_11 = buffer.data(ig0 + 11);
    const auto *ig0_16 = buffer.data(ig0 + 16);
    const auto *ig0_20 = buffer.data(ig0 + 20);
    const auto *ig0_28 = buffer.data(ig0 + 28);
    const auto *ig0_29 = buffer.data(ig0 + 29);
    const auto *ig0_34 = buffer.data(ig0 + 34);
    const auto *ig0_38 = buffer.data(ig0 + 38);
    const auto *ig0_46 = buffer.data(ig0 + 46);
    const auto *ig0_47 = buffer.data(ig0 + 47);
    const auto *ig0_52 = buffer.data(ig0 + 52);
    const auto *ig0_56 = buffer.data(ig0 + 56);
    const auto *ig0_57 = buffer.data(ig0 + 57);
    const auto *ig0_65 = buffer.data(ig0 + 65);
    const auto *ig0_68 = buffer.data(ig0 + 68);
    const auto *ig0_69 = buffer.data(ig0 + 69);
    const auto *ig0_70 = buffer.data(ig0 + 70);
    const auto *ig0_73 = buffer.data(ig0 + 73);
    const auto *ig0_79 = buffer.data(ig0 + 79);
    const auto *ig0_83 = buffer.data(ig0 + 83);
    const auto *ig0_89 = buffer.data(ig0 + 89);
    const auto *ig0_90 = buffer.data(ig0 + 90);
    const auto *ig0_92 = buffer.data(ig0 + 92);
    const auto *ig0_98 = buffer.data(ig0 + 98);
    const auto *ig0_99 = buffer.data(ig0 + 99);
    const auto *ig0_101 = buffer.data(ig0 + 101);
    const auto *ig0_107 = buffer.data(ig0 + 107);
    const auto *ig0_108 = buffer.data(ig0 + 108);
    const auto *ig0_110 = buffer.data(ig0 + 110);
    const auto *ig0_113 = buffer.data(ig0 + 113);
    const auto *ig0_122 = buffer.data(ig0 + 122);

    const auto *ig1_0 = buffer.data(ig1 + 0);
    const auto *ig1_9 = buffer.data(ig1 + 9);
    const auto *ig1_10 = buffer.data(ig1 + 10);
    const auto *ig1_11 = buffer.data(ig1 + 11);
    const auto *ig1_16 = buffer.data(ig1 + 16);
    const auto *ig1_20 = buffer.data(ig1 + 20);
    const auto *ig1_28 = buffer.data(ig1 + 28);
    const auto *ig1_29 = buffer.data(ig1 + 29);
    const auto *ig1_34 = buffer.data(ig1 + 34);
    const auto *ig1_38 = buffer.data(ig1 + 38);
    const auto *ig1_46 = buffer.data(ig1 + 46);
    const auto *ig1_47 = buffer.data(ig1 + 47);
    const auto *ig1_52 = buffer.data(ig1 + 52);
    const auto *ig1_56 = buffer.data(ig1 + 56);
    const auto *ig1_57 = buffer.data(ig1 + 57);
    const auto *ig1_65 = buffer.data(ig1 + 65);
    const auto *ig1_68 = buffer.data(ig1 + 68);
    const auto *ig1_69 = buffer.data(ig1 + 69);
    const auto *ig1_70 = buffer.data(ig1 + 70);
    const auto *ig1_73 = buffer.data(ig1 + 73);
    const auto *ig1_79 = buffer.data(ig1 + 79);
    const auto *ig1_83 = buffer.data(ig1 + 83);
    const auto *ig1_89 = buffer.data(ig1 + 89);
    const auto *ig1_90 = buffer.data(ig1 + 90);
    const auto *ig1_92 = buffer.data(ig1 + 92);
    const auto *ig1_98 = buffer.data(ig1 + 98);
    const auto *ig1_99 = buffer.data(ig1 + 99);
    const auto *ig1_101 = buffer.data(ig1 + 101);
    const auto *ig1_107 = buffer.data(ig1 + 107);
    const auto *ig1_108 = buffer.data(ig1 + 108);
    const auto *ig1_110 = buffer.data(ig1 + 110);
    const auto *ig1_113 = buffer.data(ig1 + 113);
    const auto *ig1_122 = buffer.data(ig1 + 122);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_98 = buffer.data(kf + 98);

    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);

    const auto *ld0_0 = buffer.data(ld0 + 0);
    const auto *ld0_1 = buffer.data(ld0 + 1);
    const auto *ld0_2 = buffer.data(ld0 + 2);
    const auto *ld0_3 = buffer.data(ld0 + 3);
    const auto *ld0_4 = buffer.data(ld0 + 4);
    const auto *ld0_5 = buffer.data(ld0 + 5);
    const auto *ld0_6 = buffer.data(ld0 + 6);
    const auto *ld0_7 = buffer.data(ld0 + 7);
    const auto *ld0_8 = buffer.data(ld0 + 8);
    const auto *ld0_9 = buffer.data(ld0 + 9);
    const auto *ld0_10 = buffer.data(ld0 + 10);
    const auto *ld0_11 = buffer.data(ld0 + 11);
    const auto *ld0_12 = buffer.data(ld0 + 12);
    const auto *ld0_13 = buffer.data(ld0 + 13);
    const auto *ld0_14 = buffer.data(ld0 + 14);
    const auto *ld0_15 = buffer.data(ld0 + 15);
    const auto *ld0_16 = buffer.data(ld0 + 16);
    const auto *ld0_17 = buffer.data(ld0 + 17);
    const auto *ld0_18 = buffer.data(ld0 + 18);
    const auto *ld0_19 = buffer.data(ld0 + 19);
    const auto *ld0_20 = buffer.data(ld0 + 20);
    const auto *ld0_21 = buffer.data(ld0 + 21);
    const auto *ld0_22 = buffer.data(ld0 + 22);
    const auto *ld0_23 = buffer.data(ld0 + 23);
    const auto *ld0_24 = buffer.data(ld0 + 24);
    const auto *ld0_25 = buffer.data(ld0 + 25);
    const auto *ld0_26 = buffer.data(ld0 + 26);
    const auto *ld0_27 = buffer.data(ld0 + 27);
    const auto *ld0_28 = buffer.data(ld0 + 28);
    const auto *ld0_29 = buffer.data(ld0 + 29);
    const auto *ld0_30 = buffer.data(ld0 + 30);
    const auto *ld0_31 = buffer.data(ld0 + 31);
    const auto *ld0_32 = buffer.data(ld0 + 32);
    const auto *ld0_33 = buffer.data(ld0 + 33);
    const auto *ld0_34 = buffer.data(ld0 + 34);
    const auto *ld0_35 = buffer.data(ld0 + 35);
    const auto *ld0_36 = buffer.data(ld0 + 36);
    const auto *ld0_37 = buffer.data(ld0 + 37);
    const auto *ld0_38 = buffer.data(ld0 + 38);
    const auto *ld0_39 = buffer.data(ld0 + 39);
    const auto *ld0_40 = buffer.data(ld0 + 40);
    const auto *ld0_41 = buffer.data(ld0 + 41);
    const auto *ld0_42 = buffer.data(ld0 + 42);
    const auto *ld0_43 = buffer.data(ld0 + 43);
    const auto *ld0_44 = buffer.data(ld0 + 44);
    const auto *ld0_45 = buffer.data(ld0 + 45);
    const auto *ld0_46 = buffer.data(ld0 + 46);
    const auto *ld0_47 = buffer.data(ld0 + 47);
    const auto *ld0_48 = buffer.data(ld0 + 48);
    const auto *ld0_49 = buffer.data(ld0 + 49);
    const auto *ld0_50 = buffer.data(ld0 + 50);
    const auto *ld0_51 = buffer.data(ld0 + 51);
    const auto *ld0_52 = buffer.data(ld0 + 52);
    const auto *ld0_53 = buffer.data(ld0 + 53);

    const auto *ld1_0 = buffer.data(ld1 + 0);
    const auto *ld1_1 = buffer.data(ld1 + 1);
    const auto *ld1_2 = buffer.data(ld1 + 2);
    const auto *ld1_3 = buffer.data(ld1 + 3);
    const auto *ld1_4 = buffer.data(ld1 + 4);
    const auto *ld1_5 = buffer.data(ld1 + 5);
    const auto *ld1_6 = buffer.data(ld1 + 6);
    const auto *ld1_7 = buffer.data(ld1 + 7);
    const auto *ld1_8 = buffer.data(ld1 + 8);
    const auto *ld1_9 = buffer.data(ld1 + 9);
    const auto *ld1_10 = buffer.data(ld1 + 10);
    const auto *ld1_11 = buffer.data(ld1 + 11);
    const auto *ld1_12 = buffer.data(ld1 + 12);
    const auto *ld1_13 = buffer.data(ld1 + 13);
    const auto *ld1_14 = buffer.data(ld1 + 14);
    const auto *ld1_15 = buffer.data(ld1 + 15);
    const auto *ld1_16 = buffer.data(ld1 + 16);
    const auto *ld1_17 = buffer.data(ld1 + 17);
    const auto *ld1_18 = buffer.data(ld1 + 18);
    const auto *ld1_19 = buffer.data(ld1 + 19);
    const auto *ld1_20 = buffer.data(ld1 + 20);
    const auto *ld1_21 = buffer.data(ld1 + 21);
    const auto *ld1_22 = buffer.data(ld1 + 22);
    const auto *ld1_23 = buffer.data(ld1 + 23);
    const auto *ld1_24 = buffer.data(ld1 + 24);
    const auto *ld1_25 = buffer.data(ld1 + 25);
    const auto *ld1_26 = buffer.data(ld1 + 26);
    const auto *ld1_27 = buffer.data(ld1 + 27);
    const auto *ld1_28 = buffer.data(ld1 + 28);
    const auto *ld1_29 = buffer.data(ld1 + 29);
    const auto *ld1_30 = buffer.data(ld1 + 30);
    const auto *ld1_31 = buffer.data(ld1 + 31);
    const auto *ld1_32 = buffer.data(ld1 + 32);
    const auto *ld1_33 = buffer.data(ld1 + 33);
    const auto *ld1_34 = buffer.data(ld1 + 34);
    const auto *ld1_35 = buffer.data(ld1 + 35);
    const auto *ld1_36 = buffer.data(ld1 + 36);
    const auto *ld1_37 = buffer.data(ld1 + 37);
    const auto *ld1_38 = buffer.data(ld1 + 38);
    const auto *ld1_39 = buffer.data(ld1 + 39);
    const auto *ld1_40 = buffer.data(ld1 + 40);
    const auto *ld1_41 = buffer.data(ld1 + 41);
    const auto *ld1_42 = buffer.data(ld1 + 42);
    const auto *ld1_43 = buffer.data(ld1 + 43);
    const auto *ld1_44 = buffer.data(ld1 + 44);
    const auto *ld1_45 = buffer.data(ld1 + 45);
    const auto *ld1_46 = buffer.data(ld1 + 46);
    const auto *ld1_47 = buffer.data(ld1 + 47);
    const auto *ld1_48 = buffer.data(ld1 + 48);
    const auto *ld1_49 = buffer.data(ld1 + 49);
    const auto *ld1_50 = buffer.data(ld1 + 50);
    const auto *ld1_51 = buffer.data(ld1 + 51);
    const auto *ld1_52 = buffer.data(ld1 + 52);
    const auto *ld1_53 = buffer.data(ld1 + 53);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kf_0, ld0_0, ld1_0, lf_0, \
                         lf_1, lf_2 : simd::cache_line_size())
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

        t_4[k] = f_3 * ld0_0[k]
                 - f_4 * ld1_0[k]
                 + pb_z[k] * lf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ld0_1, ld0_2, ld1_1, ld1_2, lf_3, \
                         lf_4, lf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ld0_1[k]
                 - f_2 * ld1_1[k]
                 + pb_y[k] * lf_3[k];

        t_6[k] = f_3 * ld0_2[k]
                 - f_4 * ld1_2[k]
                 + pb_y[k] * lf_4[k];

        t_7[k] = pb_y[k] * lf_5[k];

        t_8[k] = f_1 * ld0_2[k]
                 - f_2 * ld1_2[k]
                 + pb_z[k] * lf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ig0_0, ig1_0, kf_10, kg_9, ld0_4, \
                         ld1_4, lf_6, lf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ig0_0[k]
                 - f_6 * ig1_0[k]
                 + pa_y[k] * kg_9[k];

        t_10[k] = pb_z[k] * lf_6[k];

        t_11[k] = f_7 * kf_10[k]
                  + f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_x[k] * lf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, ig0_16, ig1_16, kf_11, \
                         kg_16, ld0_3, ld1_3, lf_7, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ld0_3[k]
                  - f_4 * ld1_3[k]
                  + pb_z[k] * lf_7[k];

        t_13[k] = f_7 * kf_11[k]
                  + pb_x[k] * lf_9[k];

        t_14[k] = f_8 * ig0_16[k]
                  - f_9 * ig1_16[k]
                  + pa_x[k] * kg_16[k];

        t_15[k] = pb_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, ig0_0, ig1_0, kg_10, ld0_4, ld0_5, \
                         ld1_4, ld1_5, lf_10, lf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ld0_4[k]
                  - f_4 * ld1_4[k]
                  + pb_z[k] * lf_10[k];

        t_17[k] = f_1 * ld0_5[k]
                  - f_2 * ld1_5[k]
                  + pb_z[k] * lf_11[k];

        t_18[k] = f_5 * ig0_0[k]
                  - f_6 * ig1_0[k]
                  + pa_z[k] * kg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, kf_16, kf_19, ld0_6, ld0_8, \
                         ld1_6, ld1_8, lf_12, lf_13, lf_14, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lf_12[k];

        t_20[k] = f_3 * ld0_6[k]
                  - f_4 * ld1_6[k]
                  + pb_y[k] * lf_13[k];

        t_21[k] = f_7 * kf_16[k]
                  + f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_x[k] * lf_14[k];

        t_22[k] = f_7 * kf_19[k]
                  + pb_x[k] * lf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, ig0_28, ig1_28, kg_28, ld0_7, \
                         ld0_8, ld1_7, ld1_8, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ld0_7[k]
                  - f_2 * ld1_7[k]
                  + pb_y[k] * lf_15[k];

        t_24[k] = f_3 * ld0_8[k]
                  - f_4 * ld1_8[k]
                  + pb_y[k] * lf_16[k];

        t_25[k] = pb_y[k] * lf_17[k];

        t_26[k] = f_8 * ig0_28[k]
                  - f_9 * ig1_28[k]
                  + pa_x[k] * kg_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, ig0_9, ig1_9, kf_22, kg_11, \
                         ld0_10, ld1_10, lf_18, lf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * ig0_9[k]
                  - f_11 * ig1_9[k]
                  + pa_y[k] * kg_11[k];

        t_28[k] = pb_z[k] * lf_18[k];

        t_29[k] = f_12 * kf_22[k]
                  + f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_x[k] * lf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ig0_34, ig1_34, kf_23, \
                         kg_34, ld0_9, ld1_9, lf_19, lf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ld0_9[k]
                  - f_4 * ld1_9[k]
                  + pb_z[k] * lf_19[k];

        t_31[k] = f_12 * kf_23[k]
                  + pb_x[k] * lf_21[k];

        t_32[k] = f_13 * ig0_34[k]
                  - f_14 * ig1_34[k]
                  + pa_x[k] * kg_34[k];

        t_33[k] = pb_z[k] * lf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, ig0_10, ig1_10, kg_20, ld0_10, ld0_11, \
                         ld1_10, ld1_11, lf_22, lf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ld0_10[k]
                  - f_4 * ld1_10[k]
                  + pb_z[k] * lf_22[k];

        t_35[k] = f_1 * ld0_11[k]
                  - f_2 * ld1_11[k]
                  + pb_z[k] * lf_23[k];

        t_36[k] = f_10 * ig0_10[k]
                  - f_11 * ig1_10[k]
                  + pa_z[k] * kg_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, kf_28, kf_31, ld0_12, ld0_14, \
                         ld1_12, ld1_14, lf_24, lf_25, lf_26, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * lf_24[k];

        t_38[k] = f_3 * ld0_12[k]
                  - f_4 * ld1_12[k]
                  + pb_y[k] * lf_25[k];

        t_39[k] = f_12 * kf_28[k]
                  + f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_x[k] * lf_26[k];

        t_40[k] = f_12 * kf_31[k]
                  + pb_x[k] * lf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, ig0_46, ig1_46, kg_46, ld0_13, \
                         ld0_14, ld1_13, ld1_14, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ld0_13[k]
                  - f_2 * ld1_13[k]
                  + pb_y[k] * lf_27[k];

        t_42[k] = f_3 * ld0_14[k]
                  - f_4 * ld1_14[k]
                  + pb_y[k] * lf_28[k];

        t_43[k] = pb_y[k] * lf_29[k];

        t_44[k] = f_13 * ig0_46[k]
                  - f_14 * ig1_46[k]
                  + pa_x[k] * kg_46[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_x, pb_z, ig0_11, ig1_11, kf_34, kg_29, \
                         ld0_16, ld1_16, lf_30, lf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * ig0_11[k]
                  - f_16 * ig1_11[k]
                  + pa_y[k] * kg_29[k];

        t_46[k] = pb_z[k] * lf_30[k];

        t_47[k] = f_17 * kf_34[k]
                  + f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_x[k] * lf_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, ig0_52, ig1_52, kf_35, \
                         kg_52, ld0_15, ld1_15, lf_31, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld0_15[k]
                  - f_4 * ld1_15[k]
                  + pb_z[k] * lf_31[k];

        t_49[k] = f_17 * kf_35[k]
                  + pb_x[k] * lf_33[k];

        t_50[k] = f_15 * ig0_52[k]
                  - f_16 * ig1_52[k]
                  + pa_x[k] * kg_52[k];

        t_51[k] = pb_z[k] * lf_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, ig0_56, ig1_56, kg_56, ld0_16, ld0_17, \
                         ld1_16, ld1_17, lf_34, lf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * ld0_16[k]
                  - f_4 * ld1_16[k]
                  + pb_z[k] * lf_34[k];

        t_53[k] = f_1 * ld0_17[k]
                  - f_2 * ld1_17[k]
                  + pb_z[k] * lf_35[k];

        t_54[k] = f_15 * ig0_56[k]
                  - f_16 * ig1_56[k]
                  + pa_x[k] * kg_56[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, ig0_20, ig1_20, kg_38, ld0_18, ld1_18, \
                         lf_36, lf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_15 * ig0_20[k]
                  - f_16 * ig1_20[k]
                  + pa_z[k] * kg_38[k];

        t_56[k] = pb_y[k] * lf_36[k];

        t_57[k] = f_3 * ld0_18[k]
                  - f_4 * ld1_18[k]
                  + pb_y[k] * lf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, kf_40, kf_43, ld0_19, ld0_20, \
                         ld1_19, ld1_20, lf_38, lf_39, lf_40, lf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_17 * kf_40[k]
                  + f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_x[k] * lf_38[k];

        t_59[k] = f_17 * kf_43[k]
                  + pb_x[k] * lf_41[k];

        t_60[k] = f_1 * ld0_19[k]
                  - f_2 * ld1_19[k]
                  + pb_y[k] * lf_39[k];

        t_61[k] = f_3 * ld0_20[k]
                  - f_4 * ld1_20[k]
                  + pb_y[k] * lf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, ig0_29, ig0_65, \
                         ig1_29, ig1_65, kg_47, kg_65, lf_41, lf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * lf_41[k];

        t_63[k] = f_15 * ig0_65[k]
                  - f_16 * ig1_65[k]
                  + pa_x[k] * kg_65[k];

        t_64[k] = f_13 * ig0_29[k]
                  - f_14 * ig1_29[k]
                  + pa_y[k] * kg_47[k];

        t_65[k] = pb_z[k] * lf_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kf_46, kf_47, ld0_21, ld0_22, ld1_21, \
                         ld1_22, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_18 * kf_46[k]
                  + f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_x[k] * lf_44[k];

        t_67[k] = f_3 * ld0_21[k]
                  - f_4 * ld1_21[k]
                  + pb_z[k] * lf_43[k];

        t_68[k] = f_18 * kf_47[k]
                  + pb_x[k] * lf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pb_z, ig0_68, ig1_68, kg_71, ld0_22, \
                         ld0_23, ld1_22, ld1_23, lf_45, lf_46, lf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * ig0_68[k]
                  - f_11 * ig1_68[k]
                  + pa_x[k] * kg_71[k];

        t_70[k] = pb_z[k] * lf_45[k];

        t_71[k] = f_3 * ld0_22[k]
                  - f_4 * ld1_22[k]
                  + pb_z[k] * lf_46[k];

        t_72[k] = f_1 * ld0_23[k]
                  - f_2 * ld1_23[k]
                  + pb_z[k] * lf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_z, ig0_38, ig0_69, ig0_70, ig1_38, ig1_69, \
                         ig1_70, kg_57, kg_75, kg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * ig0_69[k]
                  - f_11 * ig1_69[k]
                  + pa_x[k] * kg_75[k];

        t_74[k] = f_10 * ig0_70[k]
                  - f_11 * ig1_70[k]
                  + pa_x[k] * kg_76[k];

        t_75[k] = f_13 * ig0_38[k]
                  - f_14 * ig1_38[k]
                  + pa_z[k] * kg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, kf_52, kf_55, ld0_24, ld0_26, \
                         ld1_24, ld1_26, lf_48, lf_49, lf_50, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * lf_48[k];

        t_77[k] = f_3 * ld0_24[k]
                  - f_4 * ld1_24[k]
                  + pb_y[k] * lf_49[k];

        t_78[k] = f_18 * kf_52[k]
                  + f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_x[k] * lf_50[k];

        t_79[k] = f_18 * kf_55[k]
                  + pb_x[k] * lf_53[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pb_y, ig0_73, ig1_73, kg_85, ld0_25, \
                         ld0_26, ld1_25, ld1_26, lf_51, lf_52, lf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * ld0_25[k]
                  - f_2 * ld1_25[k]
                  + pb_y[k] * lf_51[k];

        t_81[k] = f_3 * ld0_26[k]
                  - f_4 * ld1_26[k]
                  + pb_y[k] * lf_52[k];

        t_82[k] = pb_y[k] * lf_53[k];

        t_83[k] = f_10 * ig0_73[k]
                  - f_11 * ig1_73[k]
                  + pa_x[k] * kg_85[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, pb_z, ig0_47, ig1_47, kf_56, kg_66, \
                         ld0_28, ld1_28, lf_54, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * ig0_47[k]
                  - f_9 * ig1_47[k]
                  + pa_y[k] * kg_66[k];

        t_85[k] = pb_z[k] * lf_54[k];

        t_86[k] = f_19 * kf_56[k]
                  + f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_x[k] * lf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pb_x, pb_z, ig0_79, ig1_79, kf_57, \
                         kg_86, ld0_27, ld1_27, lf_55, lf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * ld0_27[k]
                  - f_4 * ld1_27[k]
                  + pb_z[k] * lf_55[k];

        t_88[k] = f_19 * kf_57[k]
                  + pb_x[k] * lf_57[k];

        t_89[k] = f_5 * ig0_79[k]
                  - f_6 * ig1_79[k]
                  + pa_x[k] * kg_86[k];

        t_90[k] = pb_z[k] * lf_57[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pb_z, ig0_90, ig1_90, kg_87, ld0_28, ld0_29, \
                         ld1_28, ld1_29, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * ld0_28[k]
                  - f_4 * ld1_28[k]
                  + pb_z[k] * lf_58[k];

        t_92[k] = f_1 * ld0_29[k]
                  - f_2 * ld1_29[k]
                  + pb_z[k] * lf_59[k];

        t_93[k] = f_5 * ig0_90[k]
                  - f_6 * ig1_90[k]
                  + pa_x[k] * kg_87[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, ig0_57, ig0_99, ig0_108, ig1_57, \
                         ig1_99, ig1_108, kg_77, kg_88, kg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * ig0_99[k]
                  - f_6 * ig1_99[k]
                  + pa_x[k] * kg_88[k];

        t_95[k] = f_5 * ig0_108[k]
                  - f_6 * ig1_108[k]
                  + pa_x[k] * kg_89[k];

        t_96[k] = f_8 * ig0_57[k]
                  - f_9 * ig1_57[k]
                  + pa_z[k] * kg_77[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, kf_58, kf_59, ld0_30, ld0_32, \
                         ld1_30, ld1_32, lf_60, lf_61, lf_62, lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pb_y[k] * lf_60[k];

        t_98[k] = f_3 * ld0_30[k]
                  - f_4 * ld1_30[k]
                  + pb_y[k] * lf_61[k];

        t_99[k] = f_19 * kf_58[k]
                  + f_3 * ld0_32[k]
                  - f_4 * ld1_32[k]
                  + pb_x[k] * lf_62[k];

        t_100[k] = f_19 * kf_59[k]
                   + pb_x[k] * lf_65[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, ig0_122, ig1_122, kg_90, \
                         ld0_31, ld0_32, ld1_31, ld1_32, lf_63, lf_64, \
                         lf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_1 * ld0_31[k]
                   - f_2 * ld1_31[k]
                   + pb_y[k] * lf_63[k];

        t_102[k] = f_3 * ld0_32[k]
                   - f_4 * ld1_32[k]
                   + pb_y[k] * lf_64[k];

        t_103[k] = pb_y[k] * lf_65[k];

        t_104[k] = f_5 * ig0_122[k]
                   - f_6 * ig1_122[k]
                   + pa_x[k] * kg_90[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, ld0_33, ld0_34, ld0_35, ld1_33, \
                         ld1_34, ld1_35, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * ld0_33[k]
                   - f_2 * ld1_33[k]
                   + pb_x[k] * lf_66[k];

        t_106[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_x[k] * lf_67[k];

        t_107[k] = f_3 * ld0_35[k]
                   - f_4 * ld1_35[k]
                   + pb_x[k] * lf_68[k];

        t_108[k] = pb_x[k] * lf_69[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, kf_63, ld0_34, \
                         ld0_35, ld1_34, ld1_35, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * lf_71[k];

        t_110[k] = f_0 * kf_63[k]
                   + f_1 * ld0_34[k]
                   - f_2 * ld1_34[k]
                   + pb_y[k] * lf_69[k];

        t_111[k] = pb_z[k] * lf_69[k];

        t_112[k] = f_3 * ld0_34[k]
                   - f_4 * ld1_34[k]
                   + pb_z[k] * lf_70[k];

        t_113[k] = f_1 * ld0_35[k]
                   - f_2 * ld1_35[k]
                   + pb_z[k] * lf_71[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, ld0_36, ld0_37, ld0_38, ld1_36, \
                         ld1_37, ld1_38, lf_72, lf_73, lf_74, lf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ld0_36[k]
                   - f_2 * ld1_36[k]
                   + pb_x[k] * lf_72[k];

        t_115[k] = f_3 * ld0_37[k]
                   - f_4 * ld1_37[k]
                   + pb_x[k] * lf_73[k];

        t_116[k] = f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_x[k] * lf_74[k];

        t_117[k] = pb_x[k] * lf_75[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_z, pb_x, pb_y, ig0_79, ig1_79, kf_71, \
                         kf_72, kg_100, ld0_38, ld1_38, lf_76, lf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * lf_77[k];

        t_119[k] = f_5 * ig0_79[k]
                   - f_6 * ig1_79[k]
                   + pa_z[k] * kg_100[k];

        t_120[k] = f_7 * kf_71[k]
                   + f_3 * ld0_38[k]
                   - f_4 * ld1_38[k]
                   + pb_y[k] * lf_76[k];

        t_121[k] = f_7 * kf_72[k]
                   + pb_y[k] * lf_77[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pb_x, ig0_92, ig1_92, kg_109, ld0_39, \
                         ld0_40, ld1_39, ld1_40, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ig0_92[k]
                   - f_9 * ig1_92[k]
                   + pa_y[k] * kg_109[k];

        t_123[k] = f_1 * ld0_39[k]
                   - f_2 * ld1_39[k]
                   + pb_x[k] * lf_78[k];

        t_124[k] = f_3 * ld0_40[k]
                   - f_4 * ld1_40[k]
                   + pb_x[k] * lf_79[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, ig0_83, ig1_83, kg_106, \
                         ld0_41, ld1_41, lf_80, lf_81, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_x[k] * lf_80[k];

        t_126[k] = pb_x[k] * lf_81[k];

        t_127[k] = pb_x[k] * lf_83[k];

        t_128[k] = f_10 * ig0_83[k]
                   - f_11 * ig1_83[k]
                   + pa_z[k] * kg_106[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_y, ig0_101, ig1_101, kf_77, kf_78, \
                         kg_118, ld0_41, ld1_41, lf_82, lf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * kf_77[k]
                   + f_3 * ld0_41[k]
                   - f_4 * ld1_41[k]
                   + pb_y[k] * lf_82[k];

        t_130[k] = f_12 * kf_78[k]
                   + pb_y[k] * lf_83[k];

        t_131[k] = f_13 * ig0_101[k]
                   - f_14 * ig1_101[k]
                   + pa_y[k] * kg_118[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, ld0_42, ld0_43, ld0_44, ld1_42, \
                         ld1_43, ld1_44, lf_84, lf_85, lf_86, lf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ld0_42[k]
                   - f_2 * ld1_42[k]
                   + pb_x[k] * lf_84[k];

        t_133[k] = f_3 * ld0_43[k]
                   - f_4 * ld1_43[k]
                   + pb_x[k] * lf_85[k];

        t_134[k] = f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_x[k] * lf_86[k];

        t_135[k] = pb_x[k] * lf_87[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_z, pb_x, pb_y, ig0_89, ig1_89, kf_83, \
                         kf_84, kg_115, ld0_44, ld1_44, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * lf_89[k];

        t_137[k] = f_15 * ig0_89[k]
                   - f_16 * ig1_89[k]
                   + pa_z[k] * kg_115[k];

        t_138[k] = f_17 * kf_83[k]
                   + f_3 * ld0_44[k]
                   - f_4 * ld1_44[k]
                   + pb_y[k] * lf_88[k];

        t_139[k] = f_17 * kf_84[k]
                   + pb_y[k] * lf_89[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_y, pb_x, ig0_110, ig1_110, kg_127, ld0_45, \
                         ld0_46, ld1_45, ld1_46, lf_90, lf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_15 * ig0_110[k]
                   - f_16 * ig1_110[k]
                   + pa_y[k] * kg_127[k];

        t_141[k] = f_1 * ld0_45[k]
                   - f_2 * ld1_45[k]
                   + pb_x[k] * lf_90[k];

        t_142[k] = f_3 * ld0_46[k]
                   - f_4 * ld1_46[k]
                   + pb_x[k] * lf_91[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, ig0_98, ig1_98, kg_124, \
                         ld0_47, ld1_47, lf_92, lf_93, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_x[k] * lf_92[k];

        t_144[k] = pb_x[k] * lf_93[k];

        t_145[k] = pb_x[k] * lf_95[k];

        t_146[k] = f_13 * ig0_98[k]
                   - f_14 * ig1_98[k]
                   + pa_z[k] * kg_124[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_y, ig0_113, ig1_113, kf_89, kf_90, \
                         kg_136, ld0_47, ld1_47, lf_94, lf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * kf_89[k]
                   + f_3 * ld0_47[k]
                   - f_4 * ld1_47[k]
                   + pb_y[k] * lf_94[k];

        t_148[k] = f_18 * kf_90[k]
                   + pb_y[k] * lf_95[k];

        t_149[k] = f_10 * ig0_113[k]
                   - f_11 * ig1_113[k]
                   + pa_y[k] * kg_136[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_x, ld0_48, ld0_49, ld0_50, ld1_48, \
                         ld1_49, ld1_50, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * ld0_48[k]
                   - f_2 * ld1_48[k]
                   + pb_x[k] * lf_96[k];

        t_151[k] = f_3 * ld0_49[k]
                   - f_4 * ld1_49[k]
                   + pb_x[k] * lf_97[k];

        t_152[k] = f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_x[k] * lf_98[k];

        t_153[k] = pb_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_z, pb_x, pb_y, ig0_107, ig1_107, \
                         kf_91, kf_92, kg_133, ld0_50, ld1_50, lf_100, \
                         lf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * lf_101[k];

        t_155[k] = f_8 * ig0_107[k]
                   - f_9 * ig1_107[k]
                   + pa_z[k] * kg_133[k];

        t_156[k] = f_19 * kf_91[k]
                   + f_3 * ld0_50[k]
                   - f_4 * ld1_50[k]
                   + pb_y[k] * lf_100[k];

        t_157[k] = f_19 * kf_92[k]
                   + pb_y[k] * lf_101[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, ig0_122, ig1_122, kg_137, ld0_51, \
                         ld0_52, ld1_51, ld1_52, lf_102, lf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * ig0_122[k]
                   - f_6 * ig1_122[k]
                   + pa_y[k] * kg_137[k];

        t_159[k] = f_1 * ld0_51[k]
                   - f_2 * ld1_51[k]
                   + pb_x[k] * lf_102[k];

        t_160[k] = f_3 * ld0_52[k]
                   - f_4 * ld1_52[k]
                   + pb_x[k] * lf_103[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pb_x, pb_y, ld0_52, ld0_53, \
                         ld1_52, ld1_53, lf_104, lf_105, lf_106, \
                         lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_x[k] * lf_104[k];

        t_162[k] = pb_x[k] * lf_105[k];

        t_163[k] = pb_x[k] * lf_107[k];

        t_164[k] = f_1 * ld0_52[k]
                   - f_2 * ld1_52[k]
                   + pb_y[k] * lf_105[k];

        t_165[k] = f_3 * ld0_53[k]
                   - f_4 * ld1_53[k]
                   + pb_y[k] * lf_106[k];

        t_166[k] = pb_y[k] * lf_107[k];
    }

#pragma omp simd aligned(t_167, pb_z, kf_98, ld0_53, ld1_53, lf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * kf_98[k]
                   + f_1 * ld0_53[k]
                   - f_2 * ld1_53[k]
                   + pb_z[k] * lf_107[k];
    }
}

}  // namespace simdt2ceri
