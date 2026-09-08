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


#include "SimdElectronRepulsionVrrRecFL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_fl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 2.0 / p;
    const auto f_16 = 2.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 4.0 / p;

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

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);
    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_1 = buffer.data(fi0 + 1);
    const auto *fi0_2 = buffer.data(fi0 + 2);
    const auto *fi0_3 = buffer.data(fi0 + 3);
    const auto *fi0_4 = buffer.data(fi0 + 4);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_6 = buffer.data(fi0 + 6);
    const auto *fi0_7 = buffer.data(fi0 + 7);
    const auto *fi0_8 = buffer.data(fi0 + 8);
    const auto *fi0_9 = buffer.data(fi0 + 9);
    const auto *fi0_10 = buffer.data(fi0 + 10);
    const auto *fi0_11 = buffer.data(fi0 + 11);
    const auto *fi0_12 = buffer.data(fi0 + 12);
    const auto *fi0_13 = buffer.data(fi0 + 13);
    const auto *fi0_14 = buffer.data(fi0 + 14);
    const auto *fi0_15 = buffer.data(fi0 + 15);
    const auto *fi0_16 = buffer.data(fi0 + 16);
    const auto *fi0_17 = buffer.data(fi0 + 17);
    const auto *fi0_18 = buffer.data(fi0 + 18);
    const auto *fi0_19 = buffer.data(fi0 + 19);
    const auto *fi0_20 = buffer.data(fi0 + 20);
    const auto *fi0_21 = buffer.data(fi0 + 21);
    const auto *fi0_22 = buffer.data(fi0 + 22);
    const auto *fi0_23 = buffer.data(fi0 + 23);
    const auto *fi0_24 = buffer.data(fi0 + 24);
    const auto *fi0_25 = buffer.data(fi0 + 25);
    const auto *fi0_26 = buffer.data(fi0 + 26);
    const auto *fi0_27 = buffer.data(fi0 + 27);
    const auto *fi0_28 = buffer.data(fi0 + 28);
    const auto *fi0_29 = buffer.data(fi0 + 29);
    const auto *fi0_30 = buffer.data(fi0 + 30);
    const auto *fi0_31 = buffer.data(fi0 + 31);
    const auto *fi0_32 = buffer.data(fi0 + 32);
    const auto *fi0_33 = buffer.data(fi0 + 33);
    const auto *fi0_34 = buffer.data(fi0 + 34);
    const auto *fi0_35 = buffer.data(fi0 + 35);
    const auto *fi0_36 = buffer.data(fi0 + 36);
    const auto *fi0_37 = buffer.data(fi0 + 37);
    const auto *fi0_38 = buffer.data(fi0 + 38);
    const auto *fi0_39 = buffer.data(fi0 + 39);
    const auto *fi0_40 = buffer.data(fi0 + 40);
    const auto *fi0_41 = buffer.data(fi0 + 41);
    const auto *fi0_42 = buffer.data(fi0 + 42);
    const auto *fi0_43 = buffer.data(fi0 + 43);
    const auto *fi0_44 = buffer.data(fi0 + 44);
    const auto *fi0_45 = buffer.data(fi0 + 45);
    const auto *fi0_46 = buffer.data(fi0 + 46);
    const auto *fi0_47 = buffer.data(fi0 + 47);
    const auto *fi0_48 = buffer.data(fi0 + 48);
    const auto *fi0_49 = buffer.data(fi0 + 49);
    const auto *fi0_50 = buffer.data(fi0 + 50);
    const auto *fi0_51 = buffer.data(fi0 + 51);
    const auto *fi0_52 = buffer.data(fi0 + 52);
    const auto *fi0_53 = buffer.data(fi0 + 53);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_1 = buffer.data(fi1 + 1);
    const auto *fi1_2 = buffer.data(fi1 + 2);
    const auto *fi1_3 = buffer.data(fi1 + 3);
    const auto *fi1_4 = buffer.data(fi1 + 4);
    const auto *fi1_5 = buffer.data(fi1 + 5);
    const auto *fi1_6 = buffer.data(fi1 + 6);
    const auto *fi1_7 = buffer.data(fi1 + 7);
    const auto *fi1_8 = buffer.data(fi1 + 8);
    const auto *fi1_9 = buffer.data(fi1 + 9);
    const auto *fi1_10 = buffer.data(fi1 + 10);
    const auto *fi1_11 = buffer.data(fi1 + 11);
    const auto *fi1_12 = buffer.data(fi1 + 12);
    const auto *fi1_13 = buffer.data(fi1 + 13);
    const auto *fi1_14 = buffer.data(fi1 + 14);
    const auto *fi1_15 = buffer.data(fi1 + 15);
    const auto *fi1_16 = buffer.data(fi1 + 16);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_18 = buffer.data(fi1 + 18);
    const auto *fi1_19 = buffer.data(fi1 + 19);
    const auto *fi1_20 = buffer.data(fi1 + 20);
    const auto *fi1_21 = buffer.data(fi1 + 21);
    const auto *fi1_22 = buffer.data(fi1 + 22);
    const auto *fi1_23 = buffer.data(fi1 + 23);
    const auto *fi1_24 = buffer.data(fi1 + 24);
    const auto *fi1_25 = buffer.data(fi1 + 25);
    const auto *fi1_26 = buffer.data(fi1 + 26);
    const auto *fi1_27 = buffer.data(fi1 + 27);
    const auto *fi1_28 = buffer.data(fi1 + 28);
    const auto *fi1_29 = buffer.data(fi1 + 29);
    const auto *fi1_30 = buffer.data(fi1 + 30);
    const auto *fi1_31 = buffer.data(fi1 + 31);
    const auto *fi1_32 = buffer.data(fi1 + 32);
    const auto *fi1_33 = buffer.data(fi1 + 33);
    const auto *fi1_34 = buffer.data(fi1 + 34);
    const auto *fi1_35 = buffer.data(fi1 + 35);
    const auto *fi1_36 = buffer.data(fi1 + 36);
    const auto *fi1_37 = buffer.data(fi1 + 37);
    const auto *fi1_38 = buffer.data(fi1 + 38);
    const auto *fi1_39 = buffer.data(fi1 + 39);
    const auto *fi1_40 = buffer.data(fi1 + 40);
    const auto *fi1_41 = buffer.data(fi1 + 41);
    const auto *fi1_42 = buffer.data(fi1 + 42);
    const auto *fi1_43 = buffer.data(fi1 + 43);
    const auto *fi1_44 = buffer.data(fi1 + 44);
    const auto *fi1_45 = buffer.data(fi1 + 45);
    const auto *fi1_46 = buffer.data(fi1 + 46);
    const auto *fi1_47 = buffer.data(fi1 + 47);
    const auto *fi1_48 = buffer.data(fi1 + 48);
    const auto *fi1_49 = buffer.data(fi1 + 49);
    const auto *fi1_50 = buffer.data(fi1 + 50);
    const auto *fi1_51 = buffer.data(fi1 + 51);
    const auto *fi1_52 = buffer.data(fi1 + 52);
    const auto *fi1_53 = buffer.data(fi1 + 53);

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
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dk_0, fi0_0, fi1_0, \
                         fk_0, fk_1, fk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pb_y[k] * fk_0[k];

        t_2[k] = pb_z[k] * fk_0[k];

        t_3[k] = f_3 * fi0_0[k]
                 - f_4 * fi1_0[k]
                 + pb_y[k] * fk_1[k];

        t_4[k] = pb_y[k] * fk_2[k];

        t_5[k] = f_3 * fi0_0[k]
                 - f_4 * fi1_0[k]
                 + pb_z[k] * fk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, fi0_1, fi0_2, fi0_3, fi1_1, \
                         fi1_2, fi1_3, fk_3, fk_4, fk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * fi0_1[k]
                 - f_6 * fi1_1[k]
                 + pb_y[k] * fk_3[k];

        t_7[k] = pb_z[k] * fk_3[k];

        t_8[k] = pb_y[k] * fk_4[k];

        t_9[k] = f_5 * fi0_2[k]
                 - f_6 * fi1_2[k]
                 + pb_z[k] * fk_4[k];

        t_10[k] = f_7 * fi0_3[k]
                  - f_8 * fi1_3[k]
                  + pb_y[k] * fk_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, fi0_4, fi0_5, fi1_4, \
                         fi1_5, fk_5, fk_6, fk_7, fk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * fk_5[k];

        t_12[k] = f_3 * fi0_4[k]
                  - f_4 * fi1_4[k]
                  + pb_y[k] * fk_6[k];

        t_13[k] = pb_y[k] * fk_7[k];

        t_14[k] = f_7 * fi0_4[k]
                  - f_8 * fi1_4[k]
                  + pb_z[k] * fk_7[k];

        t_15[k] = f_9 * fi0_5[k]
                  - f_10 * fi1_5[k]
                  + pb_y[k] * fk_8[k];

        t_16[k] = pb_z[k] * fk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, fi0_6, fi0_7, fi1_6, fi1_7, fk_9, \
                         fk_10, fk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * fi0_6[k]
                  - f_6 * fi1_6[k]
                  + pb_y[k] * fk_9[k];

        t_18[k] = f_3 * fi0_7[k]
                  - f_4 * fi1_7[k]
                  + pb_y[k] * fk_10[k];

        t_19[k] = pb_y[k] * fk_11[k];

        t_20[k] = f_9 * fi0_7[k]
                  - f_10 * fi1_7[k]
                  + pb_z[k] * fk_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, fi0_8, fi0_9, fi0_10, fi1_8, \
                         fi1_9, fi1_10, fk_12, fk_13, fk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * fi0_8[k]
                  - f_12 * fi1_8[k]
                  + pb_y[k] * fk_12[k];

        t_22[k] = pb_z[k] * fk_12[k];

        t_23[k] = f_7 * fi0_9[k]
                  - f_8 * fi1_9[k]
                  + pb_y[k] * fk_13[k];

        t_24[k] = f_5 * fi0_10[k]
                  - f_6 * fi1_10[k]
                  + pb_y[k] * fk_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, dk_20, fi0_11, \
                         fi1_11, fk_15, fk_16, fk_17, fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * fi0_11[k]
                  - f_4 * fi1_11[k]
                  + pb_y[k] * fk_15[k];

        t_26[k] = pb_y[k] * fk_16[k];

        t_27[k] = f_11 * fi0_11[k]
                  - f_12 * fi1_11[k]
                  + pb_z[k] * fk_16[k];

        t_28[k] = f_0 * dk_20[k]
                  + pb_x[k] * fk_19[k];

        t_29[k] = pb_z[k] * fk_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, dk_22, dk_23, dk_24, dk_25, \
                         fk_18, fk_20, fk_21, fk_22, fk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dk_22[k]
                  + pb_x[k] * fk_20[k];

        t_31[k] = f_0 * dk_23[k]
                  + pb_x[k] * fk_21[k];

        t_32[k] = f_0 * dk_24[k]
                  + pb_x[k] * fk_22[k];

        t_33[k] = f_0 * dk_25[k]
                  + pb_x[k] * fk_23[k];

        t_34[k] = pb_y[k] * fk_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, dk_27, fi0_12, fi0_13, \
                         fi1_12, fi1_13, fk_19, fk_20, fk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * dk_27[k]
                  + pb_x[k] * fk_25[k];

        t_36[k] = f_1 * fi0_12[k]
                  - f_2 * fi1_12[k]
                  + pb_y[k] * fk_19[k];

        t_37[k] = pb_z[k] * fk_19[k];

        t_38[k] = f_11 * fi0_13[k]
                  - f_12 * fi1_13[k]
                  + pb_y[k] * fk_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, fi0_14, fi0_15, fi0_16, fi1_14, fi1_15, \
                         fi1_16, fk_21, fk_22, fk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * fi0_14[k]
                  - f_10 * fi1_14[k]
                  + pb_y[k] * fk_21[k];

        t_40[k] = f_7 * fi0_15[k]
                  - f_8 * fi1_15[k]
                  + pb_y[k] * fk_22[k];

        t_41[k] = f_5 * fi0_16[k]
                  - f_6 * fi1_16[k]
                  + pb_y[k] * fk_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, dk_0, dl_0, \
                         fi0_17, fi1_17, fk_24, fk_25, fk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * fi0_17[k]
                  - f_4 * fi1_17[k]
                  + pb_y[k] * fk_24[k];

        t_43[k] = pb_y[k] * fk_25[k];

        t_44[k] = f_1 * fi0_17[k]
                  - f_2 * fi1_17[k]
                  + pb_z[k] * fk_25[k];

        t_45[k] = pa_y[k] * dl_0[k];

        t_46[k] = f_13 * dk_0[k]
                  + pb_y[k] * fk_26[k];

        t_47[k] = pb_z[k] * fk_26[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_z, dk_1, dk_3, dl_1, dl_2, \
                         dl_3, fk_27, fk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * dk_1[k]
                  + pa_y[k] * dl_1[k];

        t_49[k] = pb_z[k] * fk_27[k];

        t_50[k] = pa_y[k] * dl_2[k];

        t_51[k] = f_0 * dk_3[k]
                  + pa_y[k] * dl_3[k];

        t_52[k] = pb_z[k] * fk_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_y, pb_z, dk_4, dk_5, dk_7, \
                         dl_4, dl_5, dl_6, fk_29, fk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * dk_4[k]
                  + pb_y[k] * fk_29[k];

        t_54[k] = pa_y[k] * dl_4[k];

        t_55[k] = f_15 * dk_5[k]
                  + pa_y[k] * dl_5[k];

        t_56[k] = pb_z[k] * fk_30[k];

        t_57[k] = f_14 * dk_7[k]
                  + pa_y[k] * dl_6[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, dk_8, dk_9, dk_11, \
                         dl_7, dl_8, dl_9, fk_31, fk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_13 * dk_8[k]
                  + pb_y[k] * fk_31[k];

        t_59[k] = pa_y[k] * dl_7[k];

        t_60[k] = f_16 * dk_9[k]
                  + pa_y[k] * dl_8[k];

        t_61[k] = pb_z[k] * fk_32[k];

        t_62[k] = f_0 * dk_11[k]
                  + pa_y[k] * dl_9[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, dk_12, dk_13, dk_14, \
                         dl_10, dl_11, dl_12, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * dk_12[k]
                  + pa_y[k] * dl_10[k];

        t_64[k] = f_13 * dk_13[k]
                  + pb_y[k] * fk_33[k];

        t_65[k] = pa_y[k] * dl_11[k];

        t_66[k] = f_17 * dk_14[k]
                  + pa_y[k] * dl_12[k];

        t_67[k] = pb_z[k] * fk_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pb_y, dk_16, dk_17, dk_18, dk_19, \
                         dl_13, dl_14, dl_15, dl_16, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_15 * dk_16[k]
                  + pa_y[k] * dl_13[k];

        t_69[k] = f_0 * dk_17[k]
                  + pa_y[k] * dl_14[k];

        t_70[k] = f_14 * dk_18[k]
                  + pa_y[k] * dl_15[k];

        t_71[k] = f_13 * dk_19[k]
                  + pb_y[k] * fk_35[k];

        t_72[k] = pa_y[k] * dl_16[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, dk_37, dk_38, dk_39, dk_40, \
                         fk_36, fk_37, fk_38, fk_39, fk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_14 * dk_37[k]
                  + pb_x[k] * fk_37[k];

        t_74[k] = pb_z[k] * fk_36[k];

        t_75[k] = f_14 * dk_38[k]
                  + pb_x[k] * fk_38[k];

        t_76[k] = f_14 * dk_39[k]
                  + pb_x[k] * fk_39[k];

        t_77[k] = f_14 * dk_40[k]
                  + pb_x[k] * fk_40[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, dk_20, dk_41, dk_42, \
                         dl_18, dl_19, fk_37, fk_41, fk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_14 * dk_41[k]
                  + pb_x[k] * fk_41[k];

        t_79[k] = f_14 * dk_42[k]
                  + pb_x[k] * fk_42[k];

        t_80[k] = pa_y[k] * dl_18[k];

        t_81[k] = f_18 * dk_20[k]
                  + pa_y[k] * dl_19[k];

        t_82[k] = pb_z[k] * fk_37[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, dk_22, dk_23, dk_24, dk_25, \
                         dk_26, dl_20, dl_21, dl_22, dl_23, dl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_17 * dk_22[k]
                  + pa_y[k] * dl_20[k];

        t_84[k] = f_16 * dk_23[k]
                  + pa_y[k] * dl_21[k];

        t_85[k] = f_15 * dk_24[k]
                  + pa_y[k] * dl_22[k];

        t_86[k] = f_0 * dk_25[k]
                  + pa_y[k] * dl_23[k];

        t_87[k] = f_14 * dk_26[k]
                  + pa_y[k] * dl_24[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, dk_0, dk_27, \
                         dl_0, dl_25, fk_43, fk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_13 * dk_27[k]
                  + pb_y[k] * fk_43[k];

        t_89[k] = pa_y[k] * dl_25[k];

        t_90[k] = pa_z[k] * dl_0[k];

        t_91[k] = pb_y[k] * fk_44[k];

        t_92[k] = f_13 * dk_0[k]
                  + pb_z[k] * fk_44[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, dk_2, dk_3, dl_1, \
                         dl_2, dl_3, fk_45, fk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * dl_1[k];

        t_94[k] = pb_y[k] * fk_45[k];

        t_95[k] = f_14 * dk_2[k]
                  + pa_z[k] * dl_2[k];

        t_96[k] = pa_z[k] * dl_3[k];

        t_97[k] = f_13 * dk_3[k]
                  + pb_z[k] * fk_46[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, dk_4, dk_5, dk_6, \
                         dl_4, dl_5, dl_6, fk_47, fk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * fk_47[k];

        t_99[k] = f_0 * dk_4[k]
                  + pa_z[k] * dl_4[k];

        t_100[k] = pa_z[k] * dl_5[k];

        t_101[k] = f_13 * dk_5[k]
                   + pb_z[k] * fk_48[k];

        t_102[k] = f_14 * dk_6[k]
                   + pa_z[k] * dl_6[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, dk_8, dk_9, \
                         dk_10, dl_7, dl_8, dl_9, fk_49, fk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * fk_49[k];

        t_104[k] = f_15 * dk_8[k]
                   + pa_z[k] * dl_7[k];

        t_105[k] = pa_z[k] * dl_8[k];

        t_106[k] = f_13 * dk_9[k]
                   + pb_z[k] * fk_50[k];

        t_107[k] = f_14 * dk_10[k]
                   + pa_z[k] * dl_9[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, dk_11, dk_13, \
                         dk_14, dl_10, dl_11, dl_12, fk_51, fk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * dk_11[k]
                   + pa_z[k] * dl_10[k];

        t_109[k] = pb_y[k] * fk_51[k];

        t_110[k] = f_16 * dk_13[k]
                   + pa_z[k] * dl_11[k];

        t_111[k] = pa_z[k] * dl_12[k];

        t_112[k] = f_13 * dk_14[k]
                   + pb_z[k] * fk_52[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, dk_15, dk_16, dk_17, \
                         dk_19, dl_13, dl_14, dl_15, dl_16, fk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * dk_15[k]
                   + pa_z[k] * dl_13[k];

        t_114[k] = f_0 * dk_16[k]
                   + pa_z[k] * dl_14[k];

        t_115[k] = f_15 * dk_17[k]
                   + pa_z[k] * dl_15[k];

        t_116[k] = pb_y[k] * fk_53[k];

        t_117[k] = f_17 * dk_19[k]
                   + pa_z[k] * dl_16[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, dk_53, dk_54, dk_55, \
                         dk_56, dl_17, fk_56, fk_57, fk_58, fk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * dl_17[k];

        t_119[k] = f_14 * dk_53[k]
                   + pb_x[k] * fk_56[k];

        t_120[k] = f_14 * dk_54[k]
                   + pb_x[k] * fk_57[k];

        t_121[k] = f_14 * dk_55[k]
                   + pb_x[k] * fk_58[k];

        t_122[k] = f_14 * dk_56[k]
                   + pb_x[k] * fk_59[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, pb_y, dk_57, dk_58, dl_19, \
                         fk_54, fk_60, fk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_14 * dk_57[k]
                   + pb_x[k] * fk_60[k];

        t_124[k] = pb_y[k] * fk_54[k];

        t_125[k] = f_14 * dk_58[k]
                   + pb_x[k] * fk_61[k];

        t_126[k] = pa_z[k] * dl_19[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_z, pb_z, dk_20, dk_21, dk_22, dk_23, \
                         dl_20, dl_21, dl_22, fk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * dk_20[k]
                   + pb_z[k] * fk_55[k];

        t_128[k] = f_14 * dk_21[k]
                   + pa_z[k] * dl_20[k];

        t_129[k] = f_0 * dk_22[k]
                   + pa_z[k] * dl_21[k];

        t_130[k] = f_15 * dk_23[k]
                   + pa_z[k] * dl_22[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_z, pb_y, dk_24, dk_25, dk_27, dl_23, \
                         dl_24, dl_25, fk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_16 * dk_24[k]
                   + pa_z[k] * dl_23[k];

        t_132[k] = f_17 * dk_25[k]
                   + pa_z[k] * dl_24[k];

        t_133[k] = pb_y[k] * fk_61[k];

        t_134[k] = f_18 * dk_27[k]
                   + pa_z[k] * dl_25[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pa_x, pb_y, pb_z, dk_28, dk_59, \
                         dk_61, dl_41, dl_43, fk_62, fk_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_18 * dk_59[k]
                   + pa_x[k] * dl_41[k];

        t_136[k] = f_14 * dk_28[k]
                   + pb_y[k] * fk_62[k];

        t_137[k] = pb_z[k] * fk_62[k];

        t_138[k] = f_17 * dk_61[k]
                   + pa_x[k] * dl_43[k];

        t_139[k] = pb_z[k] * fk_63[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pb_y, pb_z, dk_30, dk_62, dk_63, \
                         dl_44, dl_45, fk_64, fk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_17 * dk_62[k]
                   + pa_x[k] * dl_44[k];

        t_141[k] = f_16 * dk_63[k]
                   + pa_x[k] * dl_45[k];

        t_142[k] = pb_z[k] * fk_64[k];

        t_143[k] = f_14 * dk_30[k]
                   + pb_y[k] * fk_65[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_x, pb_z, dk_65, dk_66, dk_68, dl_46, \
                         dl_47, dl_48, fk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_16 * dk_65[k]
                   + pa_x[k] * dl_46[k];

        t_145[k] = f_15 * dk_66[k]
                   + pa_x[k] * dl_47[k];

        t_146[k] = pb_z[k] * fk_66[k];

        t_147[k] = f_15 * dk_68[k]
                   + pa_x[k] * dl_48[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pb_y, pb_z, dk_32, dk_69, dk_70, \
                         dl_49, dl_50, fk_67, fk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * dk_32[k]
                   + pb_y[k] * fk_67[k];

        t_149[k] = f_15 * dk_69[k]
                   + pa_x[k] * dl_49[k];

        t_150[k] = f_0 * dk_70[k]
                   + pa_x[k] * dl_50[k];

        t_151[k] = pb_z[k] * fk_68[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_x, pb_y, dk_34, dk_72, dk_73, dk_74, \
                         dl_51, dl_52, dl_53, fk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_0 * dk_72[k]
                   + pa_x[k] * dl_51[k];

        t_153[k] = f_0 * dk_73[k]
                   + pa_x[k] * dl_52[k];

        t_154[k] = f_14 * dk_34[k]
                   + pb_y[k] * fk_69[k];

        t_155[k] = f_0 * dk_74[k]
                   + pa_x[k] * dl_53[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pa_x, pb_z, dk_75, dk_76, dk_77, \
                         dk_78, dl_54, dl_55, dl_56, dl_57, fk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_14 * dk_75[k]
                   + pa_x[k] * dl_54[k];

        t_157[k] = pb_z[k] * fk_70[k];

        t_158[k] = f_14 * dk_76[k]
                   + pa_x[k] * dl_55[k];

        t_159[k] = f_14 * dk_77[k]
                   + pa_x[k] * dl_56[k];

        t_160[k] = f_14 * dk_78[k]
                   + pa_x[k] * dl_57[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pa_x, pb_x, pb_y, pb_z, dk_36, dk_79, \
                         dk_80, dl_58, fk_71, fk_72, fk_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_14 * dk_36[k]
                   + pb_y[k] * fk_71[k];

        t_162[k] = f_14 * dk_79[k]
                   + pa_x[k] * dl_58[k];

        t_163[k] = f_13 * dk_80[k]
                   + pb_x[k] * fk_73[k];

        t_164[k] = pb_z[k] * fk_72[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pb_x, dk_82, dk_83, dk_84, dk_85, \
                         dk_86, fk_74, fk_75, fk_76, fk_77, fk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_13 * dk_82[k]
                   + pb_x[k] * fk_74[k];

        t_166[k] = f_13 * dk_83[k]
                   + pb_x[k] * fk_75[k];

        t_167[k] = f_13 * dk_84[k]
                   + pb_x[k] * fk_76[k];

        t_168[k] = f_13 * dk_85[k]
                   + pb_x[k] * fk_77[k];

        t_169[k] = f_13 * dk_86[k]
                   + pb_x[k] * fk_78[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, pa_x, pb_x, pb_z, dk_87, \
                         dl_59, dl_60, dl_61, dl_62, fk_73, fk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * dk_87[k]
                   + pb_x[k] * fk_79[k];

        t_171[k] = pa_x[k] * dl_59[k];

        t_172[k] = pb_z[k] * fk_73[k];

        t_173[k] = pa_x[k] * dl_60[k];

        t_174[k] = pa_x[k] * dl_61[k];

        t_175[k] = pa_x[k] * dl_62[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, t_181, pa_x, pa_y, pa_z, dl_26, \
                         dl_33, dl_63, dl_64, dl_65, dl_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_x[k] * dl_63[k];

        t_177[k] = pa_x[k] * dl_64[k];

        t_178[k] = pa_x[k] * dl_65[k];

        t_179[k] = pa_x[k] * dl_66[k];

        t_180[k] = pa_y[k] * dl_33[k];

        t_181[k] = pa_z[k] * dl_26[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, pa_y, pa_z, pb_y, dk_44, dl_27, \
                         dl_28, dl_34, dl_35, fk_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = pa_y[k] * dl_34[k];

        t_183[k] = pa_z[k] * dl_27[k];

        t_184[k] = f_13 * dk_44[k]
                   + pb_y[k] * fk_80[k];

        t_185[k] = pa_y[k] * dl_35[k];

        t_186[k] = pa_z[k] * dl_28[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, dk_29, dk_46, \
                         dl_29, dl_36, fk_81, fk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_13 * dk_29[k]
                   + pb_z[k] * fk_81[k];

        t_188[k] = f_13 * dk_46[k]
                   + pb_y[k] * fk_82[k];

        t_189[k] = pa_y[k] * dl_36[k];

        t_190[k] = pa_z[k] * dl_29[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_x, pa_y, pb_y, pb_z, dk_31, dk_48, \
                         dk_94, dl_37, dl_67, fk_83, fk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * dk_31[k]
                   + pb_z[k] * fk_83[k];

        t_192[k] = f_15 * dk_94[k]
                   + pa_x[k] * dl_67[k];

        t_193[k] = f_13 * dk_48[k]
                   + pb_y[k] * fk_84[k];

        t_194[k] = pa_y[k] * dl_37[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_x, pa_z, pb_z, dk_33, dk_97, dk_98, \
                         dl_30, dl_68, dl_69, fk_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * dl_30[k];

        t_196[k] = f_13 * dk_33[k]
                   + pb_z[k] * fk_85[k];

        t_197[k] = f_0 * dk_97[k]
                   + pa_x[k] * dl_68[k];

        t_198[k] = f_0 * dk_98[k]
                   + pa_x[k] * dl_69[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, dk_35, dk_50, \
                         dl_31, dl_38, fk_86, fk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * dk_50[k]
                   + pb_y[k] * fk_86[k];

        t_200[k] = pa_y[k] * dl_38[k];

        t_201[k] = pa_z[k] * dl_31[k];

        t_202[k] = f_13 * dk_35[k]
                   + pb_z[k] * fk_87[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pa_x, pb_y, dk_52, dk_100, dk_101, \
                         dk_102, dl_70, dl_71, dl_72, fk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_14 * dk_100[k]
                   + pa_x[k] * dl_70[k];

        t_204[k] = f_14 * dk_101[k]
                   + pa_x[k] * dl_71[k];

        t_205[k] = f_14 * dk_102[k]
                   + pa_x[k] * dl_72[k];

        t_206[k] = f_13 * dk_52[k]
                   + pb_y[k] * fk_88[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_y, pa_z, pb_x, dk_104, dk_105, \
                         dk_106, dl_32, dl_39, fk_89, fk_90, fk_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_y[k] * dl_39[k];

        t_208[k] = pa_z[k] * dl_32[k];

        t_209[k] = f_13 * dk_104[k]
                   + pb_x[k] * fk_89[k];

        t_210[k] = f_13 * dk_105[k]
                   + pb_x[k] * fk_90[k];

        t_211[k] = f_13 * dk_106[k]
                   + pb_x[k] * fk_91[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pa_y, pb_x, dk_107, dk_108, \
                         dk_109, dl_40, dl_73, fk_92, fk_93, fk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_13 * dk_107[k]
                   + pb_x[k] * fk_92[k];

        t_213[k] = f_13 * dk_108[k]
                   + pb_x[k] * fk_93[k];

        t_214[k] = f_13 * dk_109[k]
                   + pb_x[k] * fk_94[k];

        t_215[k] = pa_y[k] * dl_40[k];

        t_216[k] = pa_x[k] * dl_73[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, t_222, t_223, pa_x, dl_74, dl_75, \
                         dl_76, dl_77, dl_78, dl_79, dl_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_x[k] * dl_74[k];

        t_218[k] = pa_x[k] * dl_75[k];

        t_219[k] = pa_x[k] * dl_76[k];

        t_220[k] = pa_x[k] * dl_77[k];

        t_221[k] = pa_x[k] * dl_78[k];

        t_222[k] = pa_x[k] * dl_79[k];

        t_223[k] = pa_x[k] * dl_80[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_x, pb_y, pb_z, dk_43, dk_111, \
                         dk_114, dl_81, dl_82, dl_84, fk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_x[k] * dl_81[k];

        t_225[k] = f_18 * dk_111[k]
                   + pa_x[k] * dl_82[k];

        t_226[k] = pb_y[k] * fk_95[k];

        t_227[k] = f_14 * dk_43[k]
                   + pb_z[k] * fk_95[k];

        t_228[k] = f_17 * dk_114[k]
                   + pa_x[k] * dl_84[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_x, pb_y, pb_z, dk_45, dk_115, \
                         dk_116, dl_85, dl_86, fk_96, fk_97, fk_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * fk_96[k];

        t_230[k] = f_17 * dk_115[k]
                   + pa_x[k] * dl_85[k];

        t_231[k] = f_16 * dk_116[k]
                   + pa_x[k] * dl_86[k];

        t_232[k] = f_14 * dk_45[k]
                   + pb_z[k] * fk_97[k];

        t_233[k] = pb_y[k] * fk_98[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_x, pb_z, dk_47, dk_118, dk_119, \
                         dk_120, dl_87, dl_88, dl_89, fk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_16 * dk_118[k]
                   + pa_x[k] * dl_87[k];

        t_235[k] = f_15 * dk_119[k]
                   + pa_x[k] * dl_88[k];

        t_236[k] = f_14 * dk_47[k]
                   + pb_z[k] * fk_99[k];

        t_237[k] = f_15 * dk_120[k]
                   + pa_x[k] * dl_89[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_x, pb_y, pb_z, dk_49, dk_122, dk_123, \
                         dl_90, dl_91, fk_100, fk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pb_y[k] * fk_100[k];

        t_239[k] = f_15 * dk_122[k]
                   + pa_x[k] * dl_90[k];

        t_240[k] = f_0 * dk_123[k]
                   + pa_x[k] * dl_91[k];

        t_241[k] = f_14 * dk_49[k]
                   + pb_z[k] * fk_101[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pa_x, pb_y, dk_124, dk_125, \
                         dk_127, dk_128, dl_92, dl_93, dl_94, dl_95, \
                         fk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * dk_124[k]
                   + pa_x[k] * dl_92[k];

        t_243[k] = f_0 * dk_125[k]
                   + pa_x[k] * dl_93[k];

        t_244[k] = pb_y[k] * fk_102[k];

        t_245[k] = f_0 * dk_127[k]
                   + pa_x[k] * dl_94[k];

        t_246[k] = f_14 * dk_128[k]
                   + pa_x[k] * dl_95[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_x, pb_z, dk_51, dk_129, dk_130, \
                         dk_131, dl_96, dl_97, dl_98, fk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_14 * dk_51[k]
                   + pb_z[k] * fk_103[k];

        t_248[k] = f_14 * dk_129[k]
                   + pa_x[k] * dl_96[k];

        t_249[k] = f_14 * dk_130[k]
                   + pa_x[k] * dl_97[k];

        t_250[k] = f_14 * dk_131[k]
                   + pa_x[k] * dl_98[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_x, pb_x, pb_y, dk_132, dk_133, dk_134, \
                         dl_99, fk_104, fk_106, fk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * fk_104[k];

        t_252[k] = f_14 * dk_132[k]
                   + pa_x[k] * dl_99[k];

        t_253[k] = f_13 * dk_133[k]
                   + pb_x[k] * fk_106[k];

        t_254[k] = f_13 * dk_134[k]
                   + pb_x[k] * fk_107[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pb_x, pb_y, dk_135, dk_136, \
                         dk_137, dk_138, fk_105, fk_108, fk_109, fk_110, \
                         fk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_13 * dk_135[k]
                   + pb_x[k] * fk_108[k];

        t_256[k] = f_13 * dk_136[k]
                   + pb_x[k] * fk_109[k];

        t_257[k] = f_13 * dk_137[k]
                   + pb_x[k] * fk_110[k];

        t_258[k] = f_13 * dk_138[k]
                   + pb_x[k] * fk_111[k];

        t_259[k] = pb_y[k] * fk_105[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, t_265, pa_x, pb_x, dk_140, dl_100, \
                         dl_101, dl_102, dl_103, dl_104, fk_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_13 * dk_140[k]
                   + pb_x[k] * fk_112[k];

        t_261[k] = pa_x[k] * dl_100[k];

        t_262[k] = pa_x[k] * dl_101[k];

        t_263[k] = pa_x[k] * dl_102[k];

        t_264[k] = pa_x[k] * dl_103[k];

        t_265[k] = pa_x[k] * dl_104[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, pa_x, pb_x, pb_y, dl_105, dl_106, \
                         dl_107, fi0_18, fi1_18, fk_112, fk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_x[k] * dl_105[k];

        t_267[k] = pa_x[k] * dl_106[k];

        t_268[k] = pb_y[k] * fk_112[k];

        t_269[k] = pa_x[k] * dl_107[k];

        t_270[k] = f_1 * fi0_18[k]
                   - f_2 * fi1_18[k]
                   + pb_x[k] * fk_113[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_x, pb_y, pb_z, dk_59, fi0_19, fi1_19, \
                         fk_113, fk_114, fk_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_0 * dk_59[k]
                   + pb_y[k] * fk_113[k];

        t_272[k] = pb_z[k] * fk_113[k];

        t_273[k] = f_11 * fi0_19[k]
                   - f_12 * fi1_19[k]
                   + pb_x[k] * fk_115[k];

        t_274[k] = pb_z[k] * fk_114[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_x, pb_y, pb_z, dk_62, fi0_20, fi0_21, \
                         fi1_20, fi1_21, fk_115, fk_116, fk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * fi0_20[k]
                   - f_12 * fi1_20[k]
                   + pb_x[k] * fk_116[k];

        t_276[k] = f_9 * fi0_21[k]
                   - f_10 * fi1_21[k]
                   + pb_x[k] * fk_117[k];

        t_277[k] = pb_z[k] * fk_115[k];

        t_278[k] = f_0 * dk_62[k]
                   + pb_y[k] * fk_116[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_x, pb_z, fi0_22, fi0_23, fi0_24, \
                         fi1_22, fi1_23, fi1_24, fk_117, fk_118, fk_119, \
                         fk_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_9 * fi0_22[k]
                   - f_10 * fi1_22[k]
                   + pb_x[k] * fk_118[k];

        t_280[k] = f_7 * fi0_23[k]
                   - f_8 * fi1_23[k]
                   + pb_x[k] * fk_119[k];

        t_281[k] = pb_z[k] * fk_117[k];

        t_282[k] = f_7 * fi0_24[k]
                   - f_8 * fi1_24[k]
                   + pb_x[k] * fk_120[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_x, pb_y, pb_z, dk_65, fi0_25, fi0_26, \
                         fi1_25, fi1_26, fk_118, fk_119, fk_121, \
                         fk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_0 * dk_65[k]
                   + pb_y[k] * fk_118[k];

        t_284[k] = f_7 * fi0_25[k]
                   - f_8 * fi1_25[k]
                   + pb_x[k] * fk_121[k];

        t_285[k] = f_5 * fi0_26[k]
                   - f_6 * fi1_26[k]
                   + pb_x[k] * fk_122[k];

        t_286[k] = pb_z[k] * fk_119[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pb_x, pb_y, dk_69, fi0_27, fi0_28, fi1_27, \
                         fi1_28, fk_121, fk_123, fk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_5 * fi0_27[k]
                   - f_6 * fi1_27[k]
                   + pb_x[k] * fk_123[k];

        t_288[k] = f_5 * fi0_28[k]
                   - f_6 * fi1_28[k]
                   + pb_x[k] * fk_124[k];

        t_289[k] = f_0 * dk_69[k]
                   + pb_y[k] * fk_121[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_x, pb_z, fi0_29, fi0_30, fi0_32, \
                         fi1_29, fi1_30, fi1_32, fk_122, fk_125, fk_126, \
                         fk_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_5 * fi0_29[k]
                   - f_6 * fi1_29[k]
                   + pb_x[k] * fk_125[k];

        t_291[k] = f_3 * fi0_30[k]
                   - f_4 * fi1_30[k]
                   + pb_x[k] * fk_126[k];

        t_292[k] = pb_z[k] * fk_122[k];

        t_293[k] = f_3 * fi0_32[k]
                   - f_4 * fi1_32[k]
                   + pb_x[k] * fk_127[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pb_x, pb_y, dk_74, fi0_33, fi0_34, fi1_33, \
                         fi1_34, fk_125, fk_128, fk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_3 * fi0_33[k]
                   - f_4 * fi1_33[k]
                   + pb_x[k] * fk_128[k];

        t_295[k] = f_3 * fi0_34[k]
                   - f_4 * fi1_34[k]
                   + pb_x[k] * fk_129[k];

        t_296[k] = f_0 * dk_74[k]
                   + pb_y[k] * fk_125[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, t_302, pb_x, fi0_35, fi1_35, \
                         fk_130, fk_131, fk_132, fk_133, fk_134, \
                         fk_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_3 * fi0_35[k]
                   - f_4 * fi1_35[k]
                   + pb_x[k] * fk_130[k];

        t_298[k] = pb_x[k] * fk_131[k];

        t_299[k] = pb_x[k] * fk_132[k];

        t_300[k] = pb_x[k] * fk_133[k];

        t_301[k] = pb_x[k] * fk_134[k];

        t_302[k] = pb_x[k] * fk_135[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, pb_x, pb_y, pb_z, dk_80, fi0_30, \
                         fi1_30, fk_131, fk_136, fk_137, fk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pb_x[k] * fk_136[k];

        t_304[k] = pb_x[k] * fk_137[k];

        t_305[k] = pb_x[k] * fk_138[k];

        t_306[k] = f_0 * dk_80[k]
                   + f_1 * fi0_30[k]
                   - f_2 * fi1_30[k]
                   + pb_y[k] * fk_131[k];

        t_307[k] = pb_z[k] * fk_131[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pb_z, fi0_30, fi0_31, fi0_32, fi1_30, fi1_31, \
                         fi1_32, fk_132, fk_133, fk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_3 * fi0_30[k]
                   - f_4 * fi1_30[k]
                   + pb_z[k] * fk_132[k];

        t_309[k] = f_5 * fi0_31[k]
                   - f_6 * fi1_31[k]
                   + pb_z[k] * fk_133[k];

        t_310[k] = f_7 * fi0_32[k]
                   - f_8 * fi1_32[k]
                   + pb_z[k] * fk_134[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_y, pb_z, dk_87, fi0_33, fi0_34, \
                         fi0_35, fi1_33, fi1_34, fi1_35, fk_135, fk_136, \
                         fk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_9 * fi0_33[k]
                   - f_10 * fi1_33[k]
                   + pb_z[k] * fk_135[k];

        t_312[k] = f_11 * fi0_34[k]
                   - f_12 * fi1_34[k]
                   + pb_z[k] * fk_136[k];

        t_313[k] = f_0 * dk_87[k]
                   + pb_y[k] * fk_138[k];

        t_314[k] = f_1 * fi0_35[k]
                   - f_2 * fi1_35[k]
                   + pb_z[k] * fk_138[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, pa_z, pb_y, pb_z, dk_59, dk_88, \
                         dl_41, dl_42, dl_43, fk_139, fk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * dl_41[k];

        t_316[k] = pa_z[k] * dl_42[k];

        t_317[k] = f_13 * dk_59[k]
                   + pb_z[k] * fk_139[k];

        t_318[k] = pa_z[k] * dl_43[k];

        t_319[k] = f_14 * dk_88[k]
                   + pb_y[k] * fk_140[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_z, pb_y, pb_z, dk_60, dk_61, dk_90, \
                         dl_44, dl_45, fk_141, fk_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_14 * dk_60[k]
                   + pa_z[k] * dl_44[k];

        t_321[k] = pa_z[k] * dl_45[k];

        t_322[k] = f_13 * dk_61[k]
                   + pb_z[k] * fk_141[k];

        t_323[k] = f_14 * dk_90[k]
                   + pb_y[k] * fk_142[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_z, dk_62, dk_63, dk_64, dl_46, \
                         dl_47, dl_48, fk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_0 * dk_62[k]
                   + pa_z[k] * dl_46[k];

        t_325[k] = pa_z[k] * dl_47[k];

        t_326[k] = f_13 * dk_63[k]
                   + pb_z[k] * fk_143[k];

        t_327[k] = f_14 * dk_64[k]
                   + pa_z[k] * dl_48[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_z, pb_y, pb_z, dk_65, dk_66, dk_92, \
                         dl_49, dl_50, fk_144, fk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_14 * dk_92[k]
                   + pb_y[k] * fk_144[k];

        t_329[k] = f_15 * dk_65[k]
                   + pa_z[k] * dl_49[k];

        t_330[k] = pa_z[k] * dl_50[k];

        t_331[k] = f_13 * dk_66[k]
                   + pb_z[k] * fk_145[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pa_z, pb_y, dk_67, dk_68, dk_69, \
                         dk_95, dl_51, dl_52, dl_53, dl_54, fk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_14 * dk_67[k]
                   + pa_z[k] * dl_51[k];

        t_333[k] = f_0 * dk_68[k]
                   + pa_z[k] * dl_52[k];

        t_334[k] = f_14 * dk_95[k]
                   + pb_y[k] * fk_146[k];

        t_335[k] = f_16 * dk_69[k]
                   + pa_z[k] * dl_53[k];

        t_336[k] = pa_z[k] * dl_54[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pb_z, dk_70, dk_71, dk_72, dk_73, \
                         dl_55, dl_56, dl_57, fk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_13 * dk_70[k]
                   + pb_z[k] * fk_147[k];

        t_338[k] = f_14 * dk_71[k]
                   + pa_z[k] * dl_55[k];

        t_339[k] = f_0 * dk_72[k]
                   + pa_z[k] * dl_56[k];

        t_340[k] = f_15 * dk_73[k]
                   + pa_z[k] * dl_57[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, pa_z, pb_x, pb_y, dk_74, dk_99, \
                         dl_58, fk_148, fk_149, fk_150, fk_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_14 * dk_99[k]
                   + pb_y[k] * fk_148[k];

        t_342[k] = f_17 * dk_74[k]
                   + pa_z[k] * dl_58[k];

        t_343[k] = pb_x[k] * fk_149[k];

        t_344[k] = pb_x[k] * fk_150[k];

        t_345[k] = pb_x[k] * fk_151[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, t_351, pa_z, pb_x, dl_59, fk_152, \
                         fk_153, fk_154, fk_155, fk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pb_x[k] * fk_152[k];

        t_347[k] = pb_x[k] * fk_153[k];

        t_348[k] = pb_x[k] * fk_154[k];

        t_349[k] = pb_x[k] * fk_155[k];

        t_350[k] = pb_x[k] * fk_156[k];

        t_351[k] = pa_z[k] * dl_59[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_z, pb_z, dk_80, dk_81, dk_82, dk_83, \
                         dl_60, dl_61, dl_62, fk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_13 * dk_80[k]
                   + pb_z[k] * fk_149[k];

        t_353[k] = f_14 * dk_81[k]
                   + pa_z[k] * dl_60[k];

        t_354[k] = f_0 * dk_82[k]
                   + pa_z[k] * dl_61[k];

        t_355[k] = f_15 * dk_83[k]
                   + pa_z[k] * dl_62[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_z, pb_y, dk_84, dk_85, dk_87, dk_110, \
                         dl_63, dl_64, dl_66, fk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * dk_84[k]
                   + pa_z[k] * dl_63[k];

        t_357[k] = f_17 * dk_85[k]
                   + pa_z[k] * dl_64[k];

        t_358[k] = f_14 * dk_110[k]
                   + pb_y[k] * fk_156[k];

        t_359[k] = f_18 * dk_87[k]
                   + pa_z[k] * dl_66[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, pa_y, pb_y, dk_111, dk_112, \
                         dk_113, dl_82, dl_83, dl_84, fk_157, fk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pa_y[k] * dl_82[k];

        t_361[k] = f_13 * dk_111[k]
                   + pb_y[k] * fk_157[k];

        t_362[k] = pa_y[k] * dl_83[k];

        t_363[k] = f_14 * dk_112[k]
                   + pa_y[k] * dl_84[k];

        t_364[k] = f_13 * dk_113[k]
                   + pb_y[k] * fk_158[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, pa_y, pb_y, pb_z, dk_89, dk_114, \
                         dk_115, dl_85, dl_86, dl_87, fk_159, fk_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = pa_y[k] * dl_85[k];

        t_366[k] = f_0 * dk_114[k]
                   + pa_y[k] * dl_86[k];

        t_367[k] = f_14 * dk_89[k]
                   + pb_z[k] * fk_159[k];

        t_368[k] = f_13 * dk_115[k]
                   + pb_y[k] * fk_160[k];

        t_369[k] = pa_y[k] * dl_87[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pb_y, pb_z, dk_91, dk_116, dk_117, \
                         dk_118, dl_88, dl_89, fk_161, fk_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_15 * dk_116[k]
                   + pa_y[k] * dl_88[k];

        t_371[k] = f_14 * dk_91[k]
                   + pb_z[k] * fk_161[k];

        t_372[k] = f_14 * dk_117[k]
                   + pa_y[k] * dl_89[k];

        t_373[k] = f_13 * dk_118[k]
                   + pb_y[k] * fk_162[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, pa_y, pb_z, dk_93, dk_119, dk_120, \
                         dk_121, dl_90, dl_91, dl_92, dl_93, fk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * dl_90[k];

        t_375[k] = f_16 * dk_119[k]
                   + pa_y[k] * dl_91[k];

        t_376[k] = f_14 * dk_93[k]
                   + pb_z[k] * fk_163[k];

        t_377[k] = f_0 * dk_120[k]
                   + pa_y[k] * dl_92[k];

        t_378[k] = f_14 * dk_121[k]
                   + pa_y[k] * dl_93[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pa_y, pb_y, pb_z, dk_96, dk_122, dk_123, \
                         dl_94, dl_95, fk_164, fk_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_13 * dk_122[k]
                   + pb_y[k] * fk_164[k];

        t_380[k] = pa_y[k] * dl_94[k];

        t_381[k] = f_17 * dk_123[k]
                   + pa_y[k] * dl_95[k];

        t_382[k] = f_14 * dk_96[k]
                   + pb_z[k] * fk_165[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pa_y, pb_y, dk_124, dk_125, \
                         dk_126, dk_127, dl_96, dl_97, dl_98, dl_99, \
                         fk_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_15 * dk_124[k]
                   + pa_y[k] * dl_96[k];

        t_384[k] = f_0 * dk_125[k]
                   + pa_y[k] * dl_97[k];

        t_385[k] = f_14 * dk_126[k]
                   + pa_y[k] * dl_98[k];

        t_386[k] = f_13 * dk_127[k]
                   + pb_y[k] * fk_166[k];

        t_387[k] = pa_y[k] * dl_99[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, t_393, t_394, pb_x, fk_167, \
                         fk_168, fk_169, fk_170, fk_171, fk_172, \
                         fk_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pb_x[k] * fk_167[k];

        t_389[k] = pb_x[k] * fk_168[k];

        t_390[k] = pb_x[k] * fk_169[k];

        t_391[k] = pb_x[k] * fk_170[k];

        t_392[k] = pb_x[k] * fk_171[k];

        t_393[k] = pb_x[k] * fk_172[k];

        t_394[k] = pb_x[k] * fk_173[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pb_x, pb_z, dk_103, dk_133, dk_135, \
                         dl_100, dl_102, fk_167, fk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pb_x[k] * fk_174[k];

        t_396[k] = f_18 * dk_133[k]
                   + pa_y[k] * dl_100[k];

        t_397[k] = f_14 * dk_103[k]
                   + pb_z[k] * fk_167[k];

        t_398[k] = f_17 * dk_135[k]
                   + pa_y[k] * dl_102[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, pa_y, dk_136, dk_137, dk_138, dk_139, \
                         dl_103, dl_104, dl_105, dl_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_16 * dk_136[k]
                   + pa_y[k] * dl_103[k];

        t_400[k] = f_15 * dk_137[k]
                   + pa_y[k] * dl_104[k];

        t_401[k] = f_0 * dk_138[k]
                   + pa_y[k] * dl_105[k];

        t_402[k] = f_14 * dk_139[k]
                   + pa_y[k] * dl_106[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, pa_y, pb_x, pb_y, pb_z, dk_111, \
                         dk_140, dl_107, fi0_36, fi1_36, fk_174, \
                         fk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_13 * dk_140[k]
                   + pb_y[k] * fk_174[k];

        t_404[k] = pa_y[k] * dl_107[k];

        t_405[k] = f_1 * fi0_36[k]
                   - f_2 * fi1_36[k]
                   + pb_x[k] * fk_175[k];

        t_406[k] = pb_y[k] * fk_175[k];

        t_407[k] = f_0 * dk_111[k]
                   + pb_z[k] * fk_175[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pb_x, pb_y, fi0_37, fi0_38, fi0_39, \
                         fi1_37, fi1_38, fi1_39, fk_176, fk_177, fk_178, \
                         fk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_11 * fi0_37[k]
                   - f_12 * fi1_37[k]
                   + pb_x[k] * fk_177[k];

        t_409[k] = pb_y[k] * fk_176[k];

        t_410[k] = f_11 * fi0_38[k]
                   - f_12 * fi1_38[k]
                   + pb_x[k] * fk_178[k];

        t_411[k] = f_9 * fi0_39[k]
                   - f_10 * fi1_39[k]
                   + pb_x[k] * fk_179[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pb_x, pb_y, pb_z, dk_114, fi0_40, fi0_41, \
                         fi1_40, fi1_41, fk_177, fk_178, fk_180, \
                         fk_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_0 * dk_114[k]
                   + pb_z[k] * fk_177[k];

        t_413[k] = pb_y[k] * fk_178[k];

        t_414[k] = f_9 * fi0_40[k]
                   - f_10 * fi1_40[k]
                   + pb_x[k] * fk_180[k];

        t_415[k] = f_7 * fi0_41[k]
                   - f_8 * fi1_41[k]
                   + pb_x[k] * fk_181[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pb_x, pb_y, pb_z, dk_116, fi0_42, fi0_43, \
                         fi1_42, fi1_43, fk_179, fk_180, fk_182, \
                         fk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_0 * dk_116[k]
                   + pb_z[k] * fk_179[k];

        t_417[k] = f_7 * fi0_42[k]
                   - f_8 * fi1_42[k]
                   + pb_x[k] * fk_182[k];

        t_418[k] = pb_y[k] * fk_180[k];

        t_419[k] = f_7 * fi0_43[k]
                   - f_8 * fi1_43[k]
                   + pb_x[k] * fk_183[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pb_x, pb_z, dk_119, fi0_44, fi0_45, fi1_44, \
                         fi1_45, fk_181, fk_184, fk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_5 * fi0_44[k]
                   - f_6 * fi1_44[k]
                   + pb_x[k] * fk_184[k];

        t_421[k] = f_0 * dk_119[k]
                   + pb_z[k] * fk_181[k];

        t_422[k] = f_5 * fi0_45[k]
                   - f_6 * fi1_45[k]
                   + pb_x[k] * fk_185[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, pb_x, pb_y, fi0_46, fi0_47, fi0_48, \
                         fi1_46, fi1_47, fi1_48, fk_183, fk_186, fk_187, \
                         fk_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_5 * fi0_46[k]
                   - f_6 * fi1_46[k]
                   + pb_x[k] * fk_186[k];

        t_424[k] = pb_y[k] * fk_183[k];

        t_425[k] = f_5 * fi0_47[k]
                   - f_6 * fi1_47[k]
                   + pb_x[k] * fk_187[k];

        t_426[k] = f_3 * fi0_48[k]
                   - f_4 * fi1_48[k]
                   + pb_x[k] * fk_188[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pb_x, pb_z, dk_123, fi0_49, fi0_50, fi1_49, \
                         fi1_50, fk_184, fk_189, fk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_0 * dk_123[k]
                   + pb_z[k] * fk_184[k];

        t_428[k] = f_3 * fi0_49[k]
                   - f_4 * fi1_49[k]
                   + pb_x[k] * fk_189[k];

        t_429[k] = f_3 * fi0_50[k]
                   - f_4 * fi1_50[k]
                   + pb_x[k] * fk_190[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_x, pb_y, fi0_51, fi0_53, \
                         fi1_51, fi1_53, fk_187, fk_191, fk_192, fk_193, \
                         fk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_3 * fi0_51[k]
                   - f_4 * fi1_51[k]
                   + pb_x[k] * fk_191[k];

        t_431[k] = pb_y[k] * fk_187[k];

        t_432[k] = f_3 * fi0_53[k]
                   - f_4 * fi1_53[k]
                   + pb_x[k] * fk_192[k];

        t_433[k] = pb_x[k] * fk_193[k];

        t_434[k] = pb_x[k] * fk_194[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, t_440, pb_x, fk_195, fk_196, \
                         fk_197, fk_198, fk_199, fk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pb_x[k] * fk_195[k];

        t_436[k] = pb_x[k] * fk_196[k];

        t_437[k] = pb_x[k] * fk_197[k];

        t_438[k] = pb_x[k] * fk_198[k];

        t_439[k] = pb_x[k] * fk_199[k];

        t_440[k] = pb_x[k] * fk_200[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_y, pb_z, dk_133, fi0_48, fi0_49, \
                         fi0_50, fi1_48, fi1_49, fi1_50, fk_193, fk_195, \
                         fk_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * fi0_48[k]
                   - f_2 * fi1_48[k]
                   + pb_y[k] * fk_193[k];

        t_442[k] = f_0 * dk_133[k]
                   + pb_z[k] * fk_193[k];

        t_443[k] = f_11 * fi0_49[k]
                   - f_12 * fi1_49[k]
                   + pb_y[k] * fk_195[k];

        t_444[k] = f_9 * fi0_50[k]
                   - f_10 * fi1_50[k]
                   + pb_y[k] * fk_196[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pb_y, fi0_51, fi0_52, fi0_53, fi1_51, \
                         fi1_52, fi1_53, fk_197, fk_198, fk_199, \
                         fk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_7 * fi0_51[k]
                   - f_8 * fi1_51[k]
                   + pb_y[k] * fk_197[k];

        t_446[k] = f_5 * fi0_52[k]
                   - f_6 * fi1_52[k]
                   + pb_y[k] * fk_198[k];

        t_447[k] = f_3 * fi0_53[k]
                   - f_4 * fi1_53[k]
                   + pb_y[k] * fk_199[k];

        t_448[k] = pb_y[k] * fk_200[k];
    }

#pragma omp simd aligned(t_449, pb_z, dk_140, fi0_53, fi1_53, fk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * dk_140[k]
                   + f_1 * fi0_53[k]
                   - f_2 * fi1_53[k]
                   + pb_z[k] * fk_200[k];
    }
}

auto
compute_prim_fl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 2.0 / p;
    const auto f_16 = 2.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 4.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_1 = buffer.data(fi0 + 1);
    const auto *fi0_2 = buffer.data(fi0 + 2);
    const auto *fi0_3 = buffer.data(fi0 + 3);
    const auto *fi0_4 = buffer.data(fi0 + 4);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_6 = buffer.data(fi0 + 6);
    const auto *fi0_7 = buffer.data(fi0 + 7);
    const auto *fi0_8 = buffer.data(fi0 + 8);
    const auto *fi0_9 = buffer.data(fi0 + 9);
    const auto *fi0_10 = buffer.data(fi0 + 10);
    const auto *fi0_11 = buffer.data(fi0 + 11);
    const auto *fi0_12 = buffer.data(fi0 + 12);
    const auto *fi0_13 = buffer.data(fi0 + 13);
    const auto *fi0_14 = buffer.data(fi0 + 14);
    const auto *fi0_15 = buffer.data(fi0 + 15);
    const auto *fi0_16 = buffer.data(fi0 + 16);
    const auto *fi0_17 = buffer.data(fi0 + 17);
    const auto *fi0_22 = buffer.data(fi0 + 22);
    const auto *fi0_23 = buffer.data(fi0 + 23);
    const auto *fi0_24 = buffer.data(fi0 + 24);
    const auto *fi0_25 = buffer.data(fi0 + 25);
    const auto *fi0_26 = buffer.data(fi0 + 26);
    const auto *fi0_27 = buffer.data(fi0 + 27);
    const auto *fi0_28 = buffer.data(fi0 + 28);
    const auto *fi0_29 = buffer.data(fi0 + 29);
    const auto *fi0_30 = buffer.data(fi0 + 30);
    const auto *fi0_31 = buffer.data(fi0 + 31);
    const auto *fi0_32 = buffer.data(fi0 + 32);
    const auto *fi0_33 = buffer.data(fi0 + 33);
    const auto *fi0_34 = buffer.data(fi0 + 34);
    const auto *fi0_35 = buffer.data(fi0 + 35);
    const auto *fi0_36 = buffer.data(fi0 + 36);
    const auto *fi0_37 = buffer.data(fi0 + 37);
    const auto *fi0_38 = buffer.data(fi0 + 38);
    const auto *fi0_39 = buffer.data(fi0 + 39);
    const auto *fi0_42 = buffer.data(fi0 + 42);
    const auto *fi0_43 = buffer.data(fi0 + 43);
    const auto *fi0_44 = buffer.data(fi0 + 44);
    const auto *fi0_45 = buffer.data(fi0 + 45);
    const auto *fi0_46 = buffer.data(fi0 + 46);
    const auto *fi0_47 = buffer.data(fi0 + 47);
    const auto *fi0_48 = buffer.data(fi0 + 48);
    const auto *fi0_49 = buffer.data(fi0 + 49);
    const auto *fi0_50 = buffer.data(fi0 + 50);
    const auto *fi0_51 = buffer.data(fi0 + 51);
    const auto *fi0_52 = buffer.data(fi0 + 52);
    const auto *fi0_53 = buffer.data(fi0 + 53);
    const auto *fi0_54 = buffer.data(fi0 + 54);
    const auto *fi0_55 = buffer.data(fi0 + 55);
    const auto *fi0_56 = buffer.data(fi0 + 56);
    const auto *fi0_57 = buffer.data(fi0 + 57);
    const auto *fi0_58 = buffer.data(fi0 + 58);
    const auto *fi0_59 = buffer.data(fi0 + 59);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_1 = buffer.data(fi1 + 1);
    const auto *fi1_2 = buffer.data(fi1 + 2);
    const auto *fi1_3 = buffer.data(fi1 + 3);
    const auto *fi1_4 = buffer.data(fi1 + 4);
    const auto *fi1_5 = buffer.data(fi1 + 5);
    const auto *fi1_6 = buffer.data(fi1 + 6);
    const auto *fi1_7 = buffer.data(fi1 + 7);
    const auto *fi1_8 = buffer.data(fi1 + 8);
    const auto *fi1_9 = buffer.data(fi1 + 9);
    const auto *fi1_10 = buffer.data(fi1 + 10);
    const auto *fi1_11 = buffer.data(fi1 + 11);
    const auto *fi1_12 = buffer.data(fi1 + 12);
    const auto *fi1_14 = buffer.data(fi1 + 14);
    const auto *fi1_15 = buffer.data(fi1 + 15);
    const auto *fi1_16 = buffer.data(fi1 + 16);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_18 = buffer.data(fi1 + 18);
    const auto *fi1_44 = buffer.data(fi1 + 44);
    const auto *fi1_46 = buffer.data(fi1 + 46);
    const auto *fi1_47 = buffer.data(fi1 + 47);
    const auto *fi1_48 = buffer.data(fi1 + 48);
    const auto *fi1_49 = buffer.data(fi1 + 49);
    const auto *fi1_50 = buffer.data(fi1 + 50);
    const auto *fi1_51 = buffer.data(fi1 + 51);
    const auto *fi1_52 = buffer.data(fi1 + 52);
    const auto *fi1_53 = buffer.data(fi1 + 53);
    const auto *fi1_54 = buffer.data(fi1 + 54);
    const auto *fi1_55 = buffer.data(fi1 + 55);
    const auto *fi1_56 = buffer.data(fi1 + 56);
    const auto *fi1_57 = buffer.data(fi1 + 57);
    const auto *fi1_58 = buffer.data(fi1 + 58);
    const auto *fi1_59 = buffer.data(fi1 + 59);
    const auto *fi1_60 = buffer.data(fi1 + 60);
    const auto *fi1_61 = buffer.data(fi1 + 61);
    const auto *fi1_62 = buffer.data(fi1 + 62);
    const auto *fi1_83 = buffer.data(fi1 + 83);
    const auto *fi1_85 = buffer.data(fi1 + 85);
    const auto *fi1_86 = buffer.data(fi1 + 86);
    const auto *fi1_87 = buffer.data(fi1 + 87);
    const auto *fi1_88 = buffer.data(fi1 + 88);
    const auto *fi1_89 = buffer.data(fi1 + 89);
    const auto *fi1_90 = buffer.data(fi1 + 90);
    const auto *fi1_91 = buffer.data(fi1 + 91);
    const auto *fi1_92 = buffer.data(fi1 + 92);
    const auto *fi1_93 = buffer.data(fi1 + 93);
    const auto *fi1_94 = buffer.data(fi1 + 94);
    const auto *fi1_95 = buffer.data(fi1 + 95);
    const auto *fi1_96 = buffer.data(fi1 + 96);
    const auto *fi1_97 = buffer.data(fi1 + 97);
    const auto *fi1_98 = buffer.data(fi1 + 98);
    const auto *fi1_99 = buffer.data(fi1 + 99);
    const auto *fi1_100 = buffer.data(fi1 + 100);
    const auto *fi1_101 = buffer.data(fi1 + 101);

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
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
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
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
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
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dk_0, fi0_0, fi0_1, fi1_0, \
                         fi1_1, fk_0, fk_1, fk_2, fk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = f_3 * fi0_0[k]
                 - f_4 * fi1_0[k]
                 + pb_y[k] * fk_1[k];

        t_2[k] = f_3 * fi0_0[k]
                 - f_4 * fi1_0[k]
                 + pb_z[k] * fk_2[k];

        t_3[k] = f_5 * fi0_1[k]
                 - f_6 * fi1_1[k]
                 + pb_y[k] * fk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, fi0_2, fi0_3, fi0_4, fi1_2, fi1_3, \
                         fi1_4, fk_4, fk_5, fk_6, fk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi0_2[k]
                 - f_6 * fi1_2[k]
                 + pb_z[k] * fk_4[k];

        t_5[k] = f_7 * fi0_3[k]
                 - f_8 * fi1_3[k]
                 + pb_y[k] * fk_5[k];

        t_6[k] = f_3 * fi0_4[k]
                 - f_4 * fi1_4[k]
                 + pb_y[k] * fk_6[k];

        t_7[k] = f_7 * fi0_4[k]
                 - f_8 * fi1_4[k]
                 + pb_z[k] * fk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, fi0_5, fi0_6, fi0_7, fi1_5, fi1_6, \
                         fi1_7, fk_8, fk_9, fk_10, fk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * fi0_5[k]
                 - f_10 * fi1_5[k]
                 + pb_y[k] * fk_8[k];

        t_9[k] = f_5 * fi0_6[k]
                 - f_6 * fi1_6[k]
                 + pb_y[k] * fk_9[k];

        t_10[k] = f_3 * fi0_7[k]
                  - f_4 * fi1_7[k]
                  + pb_y[k] * fk_10[k];

        t_11[k] = f_9 * fi0_7[k]
                  - f_10 * fi1_7[k]
                  + pb_z[k] * fk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, fi0_8, fi0_9, fi0_10, fi1_8, fi1_9, fi1_10, \
                         fk_12, fk_13, fk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * fi0_8[k]
                  - f_12 * fi1_8[k]
                  + pb_y[k] * fk_12[k];

        t_13[k] = f_7 * fi0_9[k]
                  - f_8 * fi1_9[k]
                  + pb_y[k] * fk_13[k];

        t_14[k] = f_5 * fi0_10[k]
                  - f_6 * fi1_10[k]
                  + pb_y[k] * fk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, dk_17, dk_23, fi0_11, \
                         fi1_11, fk_15, fk_16, fk_17, fk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * fi0_11[k]
                  - f_4 * fi1_11[k]
                  + pb_y[k] * fk_15[k];

        t_16[k] = f_11 * fi0_11[k]
                  - f_12 * fi1_11[k]
                  + pb_z[k] * fk_16[k];

        t_17[k] = f_0 * dk_17[k]
                  + pb_x[k] * fk_17[k];

        t_18[k] = f_0 * dk_23[k]
                  + pb_x[k] * fk_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, fi0_12, fi0_13, fi0_14, fi1_12, fi1_14, \
                         fi1_15, fk_17, fk_18, fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * fi0_12[k]
                  - f_2 * fi1_12[k]
                  + pb_y[k] * fk_17[k];

        t_20[k] = f_11 * fi0_13[k]
                  - f_12 * fi1_14[k]
                  + pb_y[k] * fk_18[k];

        t_21[k] = f_9 * fi0_14[k]
                  - f_10 * fi1_15[k]
                  + pb_y[k] * fk_19[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_y, pb_z, fi0_15, fi0_16, fi0_17, fi1_16, \
                         fi1_17, fi1_18, fk_20, fk_21, fk_22, fk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * fi0_15[k]
                  - f_8 * fi1_16[k]
                  + pb_y[k] * fk_20[k];

        t_23[k] = f_5 * fi0_16[k]
                  - f_6 * fi1_17[k]
                  + pb_y[k] * fk_21[k];

        t_24[k] = f_3 * fi0_17[k]
                  - f_4 * fi1_18[k]
                  + pb_y[k] * fk_22[k];

        t_25[k] = f_1 * fi0_17[k]
                  - f_2 * fi1_18[k]
                  + pb_z[k] * fk_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pb_y, dk_0, dk_1, dk_3, dk_5, \
                         dl_0, dl_1, dl_3, dl_5, fk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * dl_0[k];

        t_27[k] = f_13 * dk_0[k]
                  + pb_y[k] * fk_24[k];

        t_28[k] = f_14 * dk_1[k]
                  + pa_y[k] * dl_1[k];

        t_29[k] = f_0 * dk_3[k]
                  + pa_y[k] * dl_3[k];

        t_30[k] = f_15 * dk_5[k]
                  + pa_y[k] * dl_5[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pb_x, dk_8, dk_12, dk_17, dk_29, dl_8, \
                         dl_12, dl_17, fk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_16 * dk_8[k]
                  + pa_y[k] * dl_8[k];

        t_32[k] = f_17 * dk_12[k]
                  + pa_y[k] * dl_12[k];

        t_33[k] = f_14 * dk_29[k]
                  + pb_x[k] * fk_29[k];

        t_34[k] = f_18 * dk_17[k]
                  + pa_y[k] * dl_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, pb_z, dk_0, dk_2, dk_4, dk_6, \
                         dl_0, dl_2, dl_4, dl_6, fk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * dl_0[k];

        t_36[k] = f_13 * dk_0[k]
                  + pb_z[k] * fk_30[k];

        t_37[k] = f_14 * dk_2[k]
                  + pa_z[k] * dl_2[k];

        t_38[k] = f_0 * dk_4[k]
                  + pa_z[k] * dl_4[k];

        t_39[k] = f_14 * dk_6[k]
                  + pa_z[k] * dl_6[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_z, dk_7, dk_9, dk_10, dk_11, dk_13, \
                         dl_7, dl_9, dl_10, dl_11, dl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_15 * dk_7[k]
                  + pa_z[k] * dl_7[k];

        t_41[k] = f_14 * dk_9[k]
                  + pa_z[k] * dl_9[k];

        t_42[k] = f_0 * dk_10[k]
                  + pa_z[k] * dl_10[k];

        t_43[k] = f_16 * dk_11[k]
                  + pa_z[k] * dl_11[k];

        t_44[k] = f_14 * dk_13[k]
                  + pa_z[k] * dl_13[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_z, pb_x, dk_14, dk_15, dk_16, dk_36, \
                         dl_14, dl_15, dl_16, fk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * dk_14[k]
                  + pa_z[k] * dl_14[k];

        t_46[k] = f_15 * dk_15[k]
                  + pa_z[k] * dl_15[k];

        t_47[k] = f_17 * dk_16[k]
                  + pa_z[k] * dl_16[k];

        t_48[k] = f_14 * dk_36[k]
                  + pb_x[k] * fk_40[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_z, dk_18, dk_19, dk_20, dk_21, \
                         dk_22, dl_18, dl_19, dl_20, dl_21, dl_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_14 * dk_18[k]
                  + pa_z[k] * dl_18[k];

        t_50[k] = f_0 * dk_19[k]
                  + pa_z[k] * dl_19[k];

        t_51[k] = f_15 * dk_20[k]
                  + pa_z[k] * dl_20[k];

        t_52[k] = f_16 * dk_21[k]
                  + pa_z[k] * dl_21[k];

        t_53[k] = f_17 * dk_22[k]
                  + pa_z[k] * dl_22[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pa_z, pb_y, dk_23, dk_24, dk_37, dk_39, \
                         dl_23, dl_24, dl_25, fk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_18 * dk_23[k]
                  + pa_z[k] * dl_23[k];

        t_55[k] = f_18 * dk_37[k]
                  + pa_x[k] * dl_24[k];

        t_56[k] = f_14 * dk_24[k]
                  + pb_y[k] * fk_41[k];

        t_57[k] = f_17 * dk_39[k]
                  + pa_x[k] * dl_25[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, dk_41, dk_44, dk_48, dk_53, dl_27, \
                         dl_29, dl_32, dl_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_16 * dk_41[k]
                  + pa_x[k] * dl_27[k];

        t_59[k] = f_15 * dk_44[k]
                  + pa_x[k] * dl_29[k];

        t_60[k] = f_0 * dk_48[k]
                  + pa_x[k] * dl_32[k];

        t_61[k] = f_14 * dk_53[k]
                  + pa_x[k] * dl_36[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pb_x, pb_z, dk_30, dk_54, dk_70, dl_41, \
                         dl_48, fk_46, fk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_13 * dk_54[k]
                  + pb_x[k] * fk_46[k];

        t_63[k] = pa_x[k] * dl_41[k];

        t_64[k] = f_18 * dk_70[k]
                  + pa_x[k] * dl_48[k];

        t_65[k] = f_14 * dk_30[k]
                  + pb_z[k] * fk_47[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_x, dk_74, dk_77, dk_81, dk_86, \
                         dk_87, dl_50, dl_52, dl_55, dl_59, dl_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_17 * dk_74[k]
                  + pa_x[k] * dl_50[k];

        t_67[k] = f_16 * dk_77[k]
                  + pa_x[k] * dl_52[k];

        t_68[k] = f_15 * dk_81[k]
                  + pa_x[k] * dl_55[k];

        t_69[k] = f_0 * dk_86[k]
                  + pa_x[k] * dl_59[k];

        t_70[k] = f_14 * dk_87[k]
                  + pa_x[k] * dl_64[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_y, dk_37, dk_95, dl_71, \
                         fi0_22, fi1_44, fk_53, fk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_13 * dk_95[k]
                  + pb_x[k] * fk_53[k];

        t_72[k] = pa_x[k] * dl_71[k];

        t_73[k] = f_1 * fi0_22[k]
                  - f_2 * fi1_44[k]
                  + pb_x[k] * fk_54[k];

        t_74[k] = f_0 * dk_37[k]
                  + pb_y[k] * fk_54[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, fi0_23, fi0_24, fi0_25, fi1_46, fi1_47, \
                         fi1_48, fk_55, fk_56, fk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_11 * fi0_23[k]
                  - f_12 * fi1_46[k]
                  + pb_x[k] * fk_55[k];

        t_76[k] = f_11 * fi0_24[k]
                  - f_12 * fi1_47[k]
                  + pb_x[k] * fk_56[k];

        t_77[k] = f_9 * fi0_25[k]
                  - f_10 * fi1_48[k]
                  + pb_x[k] * fk_57[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_x, fi0_26, fi0_27, fi0_28, fi1_49, fi1_50, \
                         fi1_51, fk_58, fk_59, fk_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * fi0_26[k]
                  - f_10 * fi1_49[k]
                  + pb_x[k] * fk_58[k];

        t_79[k] = f_7 * fi0_27[k]
                  - f_8 * fi1_50[k]
                  + pb_x[k] * fk_59[k];

        t_80[k] = f_7 * fi0_28[k]
                  - f_8 * fi1_51[k]
                  + pb_x[k] * fk_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, fi0_29, fi0_30, fi0_31, fi1_52, fi1_53, \
                         fi1_54, fk_61, fk_62, fk_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_7 * fi0_29[k]
                  - f_8 * fi1_52[k]
                  + pb_x[k] * fk_61[k];

        t_82[k] = f_5 * fi0_30[k]
                  - f_6 * fi1_53[k]
                  + pb_x[k] * fk_62[k];

        t_83[k] = f_5 * fi0_31[k]
                  - f_6 * fi1_54[k]
                  + pb_x[k] * fk_63[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, fi0_32, fi0_33, fi0_34, fi1_55, fi1_56, \
                         fi1_57, fk_64, fk_65, fk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * fi0_32[k]
                  - f_6 * fi1_55[k]
                  + pb_x[k] * fk_64[k];

        t_85[k] = f_5 * fi0_33[k]
                  - f_6 * fi1_56[k]
                  + pb_x[k] * fk_65[k];

        t_86[k] = f_3 * fi0_34[k]
                  - f_4 * fi1_57[k]
                  + pb_x[k] * fk_66[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, fi0_36, fi0_37, fi0_38, fi1_59, fi1_60, \
                         fi1_61, fk_67, fk_68, fk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * fi0_36[k]
                  - f_4 * fi1_59[k]
                  + pb_x[k] * fk_67[k];

        t_88[k] = f_3 * fi0_37[k]
                  - f_4 * fi1_60[k]
                  + pb_x[k] * fk_68[k];

        t_89[k] = f_3 * fi0_38[k]
                  - f_4 * fi1_61[k]
                  + pb_x[k] * fk_69[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_x, pb_y, pb_z, dk_54, fi0_34, fi0_39, fi1_57, \
                         fi1_62, fk_70, fk_71, fk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * fi0_39[k]
                  - f_4 * fi1_62[k]
                  + pb_x[k] * fk_70[k];

        t_91[k] = f_0 * dk_54[k]
                  + f_1 * fi0_34[k]
                  - f_2 * fi1_57[k]
                  + pb_y[k] * fk_71[k];

        t_92[k] = f_3 * fi0_34[k]
                  - f_4 * fi1_57[k]
                  + pb_z[k] * fk_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_z, fi0_35, fi0_36, fi0_37, fi1_58, fi1_59, \
                         fi1_60, fk_73, fk_74, fk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_5 * fi0_35[k]
                  - f_6 * fi1_58[k]
                  + pb_z[k] * fk_73[k];

        t_94[k] = f_7 * fi0_36[k]
                  - f_8 * fi1_59[k]
                  + pb_z[k] * fk_74[k];

        t_95[k] = f_9 * fi0_37[k]
                  - f_10 * fi1_60[k]
                  + pb_z[k] * fk_75[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_y, pb_z, dk_38, dk_61, dl_26, \
                         fi0_38, fi0_39, fi1_61, fi1_62, fk_76, fk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_11 * fi0_38[k]
                  - f_12 * fi1_61[k]
                  + pb_z[k] * fk_76[k];

        t_97[k] = f_0 * dk_61[k]
                  + pb_y[k] * fk_78[k];

        t_98[k] = f_1 * fi0_39[k]
                  - f_2 * fi1_62[k]
                  + pb_z[k] * fk_78[k];

        t_99[k] = f_14 * dk_38[k]
                  + pa_z[k] * dl_26[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pa_z, dk_40, dk_42, dk_43, dk_45, \
                         dk_46, dl_28, dl_30, dl_31, dl_33, dl_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_0 * dk_40[k]
                   + pa_z[k] * dl_28[k];

        t_101[k] = f_14 * dk_42[k]
                   + pa_z[k] * dl_30[k];

        t_102[k] = f_15 * dk_43[k]
                   + pa_z[k] * dl_31[k];

        t_103[k] = f_14 * dk_45[k]
                   + pa_z[k] * dl_33[k];

        t_104[k] = f_0 * dk_46[k]
                   + pa_z[k] * dl_34[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_z, dk_47, dk_49, dk_50, dk_51, \
                         dk_52, dl_35, dl_37, dl_38, dl_39, dl_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_16 * dk_47[k]
                   + pa_z[k] * dl_35[k];

        t_106[k] = f_14 * dk_49[k]
                   + pa_z[k] * dl_37[k];

        t_107[k] = f_0 * dk_50[k]
                   + pa_z[k] * dl_38[k];

        t_108[k] = f_15 * dk_51[k]
                   + pa_z[k] * dl_39[k];

        t_109[k] = f_17 * dk_52[k]
                   + pa_z[k] * dl_40[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_z, pb_z, dk_54, dk_55, dk_56, \
                         dk_57, dl_41, dl_42, dl_43, dl_44, fk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * dl_41[k];

        t_111[k] = f_13 * dk_54[k]
                   + pb_z[k] * fk_83[k];

        t_112[k] = f_14 * dk_55[k]
                   + pa_z[k] * dl_42[k];

        t_113[k] = f_0 * dk_56[k]
                   + pa_z[k] * dl_43[k];

        t_114[k] = f_15 * dk_57[k]
                   + pa_z[k] * dl_44[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_z, pb_y, dk_58, dk_59, dk_61, dk_69, \
                         dl_45, dl_46, dl_47, fk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_16 * dk_58[k]
                   + pa_z[k] * dl_45[k];

        t_116[k] = f_17 * dk_59[k]
                   + pa_z[k] * dl_46[k];

        t_117[k] = f_14 * dk_69[k]
                   + pb_y[k] * fk_90[k];

        t_118[k] = f_18 * dk_61[k]
                   + pa_z[k] * dl_47[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pa_y, dk_71, dk_73, dk_75, dk_76, \
                         dk_78, dl_49, dl_51, dl_53, dl_54, dl_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_14 * dk_71[k]
                   + pa_y[k] * dl_49[k];

        t_120[k] = f_0 * dk_73[k]
                   + pa_y[k] * dl_51[k];

        t_121[k] = f_15 * dk_75[k]
                   + pa_y[k] * dl_53[k];

        t_122[k] = f_14 * dk_76[k]
                   + pa_y[k] * dl_54[k];

        t_123[k] = f_16 * dk_78[k]
                   + pa_y[k] * dl_56[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_y, dk_79, dk_80, dk_82, dk_83, \
                         dk_84, dl_57, dl_58, dl_60, dl_61, dl_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_0 * dk_79[k]
                   + pa_y[k] * dl_57[k];

        t_125[k] = f_14 * dk_80[k]
                   + pa_y[k] * dl_58[k];

        t_126[k] = f_17 * dk_82[k]
                   + pa_y[k] * dl_60[k];

        t_127[k] = f_15 * dk_83[k]
                   + pa_y[k] * dl_61[k];

        t_128[k] = f_0 * dk_84[k]
                   + pa_y[k] * dl_62[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pb_z, dk_62, dk_85, dk_88, dk_90, \
                         dl_63, dl_65, dl_66, fk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_14 * dk_85[k]
                   + pa_y[k] * dl_63[k];

        t_130[k] = f_18 * dk_88[k]
                   + pa_y[k] * dl_65[k];

        t_131[k] = f_14 * dk_62[k]
                   + pb_z[k] * fk_95[k];

        t_132[k] = f_17 * dk_90[k]
                   + pa_y[k] * dl_66[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, dk_91, dk_92, dk_93, dk_94, dl_67, \
                         dl_68, dl_69, dl_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_16 * dk_91[k]
                   + pa_y[k] * dl_67[k];

        t_134[k] = f_15 * dk_92[k]
                   + pa_y[k] * dl_68[k];

        t_135[k] = f_0 * dk_93[k]
                   + pa_y[k] * dl_69[k];

        t_136[k] = f_14 * dk_94[k]
                   + pa_y[k] * dl_70[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_y, pb_x, pb_y, pb_z, dk_70, dk_95, \
                         dl_71, fi0_42, fi1_83, fk_102, fk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_13 * dk_95[k]
                   + pb_y[k] * fk_102[k];

        t_138[k] = pa_y[k] * dl_71[k];

        t_139[k] = f_1 * fi0_42[k]
                   - f_2 * fi1_83[k]
                   + pb_x[k] * fk_103[k];

        t_140[k] = f_0 * dk_70[k]
                   + pb_z[k] * fk_103[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pb_x, fi0_43, fi0_44, fi0_45, fi1_85, fi1_86, \
                         fi1_87, fk_105, fk_106, fk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * fi0_43[k]
                   - f_12 * fi1_85[k]
                   + pb_x[k] * fk_105[k];

        t_142[k] = f_11 * fi0_44[k]
                   - f_12 * fi1_86[k]
                   + pb_x[k] * fk_106[k];

        t_143[k] = f_9 * fi0_45[k]
                   - f_10 * fi1_87[k]
                   + pb_x[k] * fk_107[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pb_x, fi0_46, fi0_47, fi0_48, fi1_88, fi1_89, \
                         fi1_90, fk_108, fk_109, fk_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_9 * fi0_46[k]
                   - f_10 * fi1_88[k]
                   + pb_x[k] * fk_108[k];

        t_145[k] = f_7 * fi0_47[k]
                   - f_8 * fi1_89[k]
                   + pb_x[k] * fk_109[k];

        t_146[k] = f_7 * fi0_48[k]
                   - f_8 * fi1_90[k]
                   + pb_x[k] * fk_110[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_x, fi0_49, fi0_50, fi0_51, fi1_91, fi1_92, \
                         fi1_93, fk_111, fk_112, fk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_7 * fi0_49[k]
                   - f_8 * fi1_91[k]
                   + pb_x[k] * fk_111[k];

        t_148[k] = f_5 * fi0_50[k]
                   - f_6 * fi1_92[k]
                   + pb_x[k] * fk_112[k];

        t_149[k] = f_5 * fi0_51[k]
                   - f_6 * fi1_93[k]
                   + pb_x[k] * fk_113[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_x, fi0_52, fi0_53, fi0_54, fi1_94, fi1_95, \
                         fi1_96, fk_114, fk_115, fk_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * fi0_52[k]
                   - f_6 * fi1_94[k]
                   + pb_x[k] * fk_114[k];

        t_151[k] = f_5 * fi0_53[k]
                   - f_6 * fi1_95[k]
                   + pb_x[k] * fk_115[k];

        t_152[k] = f_3 * fi0_54[k]
                   - f_4 * fi1_96[k]
                   + pb_x[k] * fk_116[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_x, fi0_55, fi0_56, fi0_57, fi1_97, fi1_98, \
                         fi1_99, fk_117, fk_118, fk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_3 * fi0_55[k]
                   - f_4 * fi1_97[k]
                   + pb_x[k] * fk_117[k];

        t_154[k] = f_3 * fi0_56[k]
                   - f_4 * fi1_98[k]
                   + pb_x[k] * fk_118[k];

        t_155[k] = f_3 * fi0_57[k]
                   - f_4 * fi1_99[k]
                   + pb_x[k] * fk_119[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_y, pb_z, dk_88, fi0_54, fi0_59, fi1_96, \
                         fi1_101, fk_120, fk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_3 * fi0_59[k]
                   - f_4 * fi1_101[k]
                   + pb_x[k] * fk_120[k];

        t_157[k] = f_1 * fi0_54[k]
                   - f_2 * fi1_96[k]
                   + pb_y[k] * fk_121[k];

        t_158[k] = f_0 * dk_88[k]
                   + pb_z[k] * fk_121[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_y, fi0_55, fi0_56, fi0_57, fi1_97, fi1_98, \
                         fi1_99, fk_123, fk_124, fk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_11 * fi0_55[k]
                   - f_12 * fi1_97[k]
                   + pb_y[k] * fk_123[k];

        t_160[k] = f_9 * fi0_56[k]
                   - f_10 * fi1_98[k]
                   + pb_y[k] * fk_124[k];

        t_161[k] = f_7 * fi0_57[k]
                   - f_8 * fi1_99[k]
                   + pb_y[k] * fk_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_y, pb_z, dk_95, fi0_58, fi0_59, fi1_100, \
                         fi1_101, fk_126, fk_127, fk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * fi0_58[k]
                   - f_6 * fi1_100[k]
                   + pb_y[k] * fk_126[k];

        t_163[k] = f_3 * fi0_59[k]
                   - f_4 * fi1_101[k]
                   + pb_y[k] * fk_127[k];

        t_164[k] = f_0 * dk_95[k]
                   + f_1 * fi0_59[k]
                   - f_2 * fi1_101[k]
                   + pb_z[k] * fk_128[k];
    }
}

auto
compute_prim_fl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_8 = buffer.data(fi0 + 8);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_5 = buffer.data(fi1 + 5);
    const auto *fi1_8 = buffer.data(fi1 + 8);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_8 = buffer.data(fk + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_1, dl_1, dl_2, fi0_5, fi1_5, \
                         fk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_1[k]
                 + f_1 * fi0_5[k]
                 - f_2 * fi1_5[k]
                 + pb_y[k] * fk_5[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_2, fi0_8, fi1_8, fk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_2[k]
                 + f_1 * fi0_8[k]
                 - f_2 * fi1_8[k]
                 + pb_z[k] * fk_8[k];
    }
}

auto
compute_prim_fl_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_8 = buffer.data(fi0 + 8);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_58 = buffer.data(fi1 + 58);
    const auto *fi1_101 = buffer.data(fi1 + 101);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_130 = buffer.data(fk + 130);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_1, dl_1, dl_2, fi0_5, fi1_58, \
                         fk_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_1[k]
                 + f_1 * fi0_5[k]
                 - f_2 * fi1_58[k]
                 + pb_y[k] * fk_73[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_2, fi0_8, fi1_101, fk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_2[k]
                 + f_1 * fi0_8[k]
                 - f_2 * fi1_101[k]
                 + pb_z[k] * fk_130[k];
    }
}

auto
compute_prim_fl_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);

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

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_44 = buffer.data(dk + 44);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_1 = buffer.data(fi0 + 1);
    const auto *fi0_2 = buffer.data(fi0 + 2);
    const auto *fi0_3 = buffer.data(fi0 + 3);
    const auto *fi0_4 = buffer.data(fi0 + 4);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_7 = buffer.data(fi0 + 7);
    const auto *fi0_8 = buffer.data(fi0 + 8);
    const auto *fi0_9 = buffer.data(fi0 + 9);
    const auto *fi0_11 = buffer.data(fi0 + 11);
    const auto *fi0_12 = buffer.data(fi0 + 12);
    const auto *fi0_13 = buffer.data(fi0 + 13);
    const auto *fi0_14 = buffer.data(fi0 + 14);
    const auto *fi0_16 = buffer.data(fi0 + 16);
    const auto *fi0_17 = buffer.data(fi0 + 17);
    const auto *fi0_18 = buffer.data(fi0 + 18);
    const auto *fi0_19 = buffer.data(fi0 + 19);
    const auto *fi0_20 = buffer.data(fi0 + 20);
    const auto *fi0_43 = buffer.data(fi0 + 43);
    const auto *fi0_45 = buffer.data(fi0 + 45);
    const auto *fi0_46 = buffer.data(fi0 + 46);
    const auto *fi0_47 = buffer.data(fi0 + 47);
    const auto *fi0_49 = buffer.data(fi0 + 49);
    const auto *fi0_50 = buffer.data(fi0 + 50);
    const auto *fi0_52 = buffer.data(fi0 + 52);
    const auto *fi0_53 = buffer.data(fi0 + 53);
    const auto *fi0_54 = buffer.data(fi0 + 54);
    const auto *fi0_55 = buffer.data(fi0 + 55);
    const auto *fi0_56 = buffer.data(fi0 + 56);
    const auto *fi0_57 = buffer.data(fi0 + 57);
    const auto *fi0_58 = buffer.data(fi0 + 58);
    const auto *fi0_59 = buffer.data(fi0 + 59);
    const auto *fi0_60 = buffer.data(fi0 + 60);
    const auto *fi0_61 = buffer.data(fi0 + 61);
    const auto *fi0_62 = buffer.data(fi0 + 62);
    const auto *fi0_63 = buffer.data(fi0 + 63);
    const auto *fi0_81 = buffer.data(fi0 + 81);
    const auto *fi0_83 = buffer.data(fi0 + 83);
    const auto *fi0_84 = buffer.data(fi0 + 84);
    const auto *fi0_85 = buffer.data(fi0 + 85);
    const auto *fi0_87 = buffer.data(fi0 + 87);
    const auto *fi0_88 = buffer.data(fi0 + 88);
    const auto *fi0_89 = buffer.data(fi0 + 89);
    const auto *fi0_91 = buffer.data(fi0 + 91);
    const auto *fi0_92 = buffer.data(fi0 + 92);
    const auto *fi0_93 = buffer.data(fi0 + 93);
    const auto *fi0_94 = buffer.data(fi0 + 94);
    const auto *fi0_95 = buffer.data(fi0 + 95);
    const auto *fi0_96 = buffer.data(fi0 + 96);
    const auto *fi0_97 = buffer.data(fi0 + 97);
    const auto *fi0_98 = buffer.data(fi0 + 98);
    const auto *fi0_99 = buffer.data(fi0 + 99);
    const auto *fi0_100 = buffer.data(fi0 + 100);
    const auto *fi0_101 = buffer.data(fi0 + 101);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_1 = buffer.data(fi1 + 1);
    const auto *fi1_2 = buffer.data(fi1 + 2);
    const auto *fi1_3 = buffer.data(fi1 + 3);
    const auto *fi1_4 = buffer.data(fi1 + 4);
    const auto *fi1_5 = buffer.data(fi1 + 5);
    const auto *fi1_6 = buffer.data(fi1 + 6);
    const auto *fi1_7 = buffer.data(fi1 + 7);
    const auto *fi1_8 = buffer.data(fi1 + 8);
    const auto *fi1_9 = buffer.data(fi1 + 9);
    const auto *fi1_10 = buffer.data(fi1 + 10);
    const auto *fi1_11 = buffer.data(fi1 + 11);
    const auto *fi1_12 = buffer.data(fi1 + 12);
    const auto *fi1_14 = buffer.data(fi1 + 14);
    const auto *fi1_15 = buffer.data(fi1 + 15);
    const auto *fi1_16 = buffer.data(fi1 + 16);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_18 = buffer.data(fi1 + 18);
    const auto *fi1_35 = buffer.data(fi1 + 35);
    const auto *fi1_37 = buffer.data(fi1 + 37);
    const auto *fi1_38 = buffer.data(fi1 + 38);
    const auto *fi1_39 = buffer.data(fi1 + 39);
    const auto *fi1_40 = buffer.data(fi1 + 40);
    const auto *fi1_41 = buffer.data(fi1 + 41);
    const auto *fi1_42 = buffer.data(fi1 + 42);
    const auto *fi1_43 = buffer.data(fi1 + 43);
    const auto *fi1_44 = buffer.data(fi1 + 44);
    const auto *fi1_45 = buffer.data(fi1 + 45);
    const auto *fi1_46 = buffer.data(fi1 + 46);
    const auto *fi1_47 = buffer.data(fi1 + 47);
    const auto *fi1_48 = buffer.data(fi1 + 48);
    const auto *fi1_49 = buffer.data(fi1 + 49);
    const auto *fi1_50 = buffer.data(fi1 + 50);
    const auto *fi1_51 = buffer.data(fi1 + 51);
    const auto *fi1_52 = buffer.data(fi1 + 52);
    const auto *fi1_53 = buffer.data(fi1 + 53);
    const auto *fi1_62 = buffer.data(fi1 + 62);
    const auto *fi1_64 = buffer.data(fi1 + 64);
    const auto *fi1_65 = buffer.data(fi1 + 65);
    const auto *fi1_66 = buffer.data(fi1 + 66);
    const auto *fi1_67 = buffer.data(fi1 + 67);
    const auto *fi1_68 = buffer.data(fi1 + 68);
    const auto *fi1_69 = buffer.data(fi1 + 69);
    const auto *fi1_70 = buffer.data(fi1 + 70);
    const auto *fi1_71 = buffer.data(fi1 + 71);
    const auto *fi1_72 = buffer.data(fi1 + 72);
    const auto *fi1_73 = buffer.data(fi1 + 73);
    const auto *fi1_74 = buffer.data(fi1 + 74);
    const auto *fi1_75 = buffer.data(fi1 + 75);
    const auto *fi1_76 = buffer.data(fi1 + 76);
    const auto *fi1_77 = buffer.data(fi1 + 77);
    const auto *fi1_78 = buffer.data(fi1 + 78);
    const auto *fi1_79 = buffer.data(fi1 + 79);
    const auto *fi1_80 = buffer.data(fi1 + 80);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dk_0, fi0_0, fi0_1, fi1_0, \
                         fi1_1, fk_0, fk_1, fk_2, fk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = f_3 * fi0_0[k]
                 - f_4 * fi1_0[k]
                 + pb_y[k] * fk_1[k];

        t_2[k] = f_3 * fi0_0[k]
                 - f_4 * fi1_0[k]
                 + pb_z[k] * fk_2[k];

        t_3[k] = f_5 * fi0_1[k]
                 - f_6 * fi1_1[k]
                 + pb_y[k] * fk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, fi0_2, fi0_3, fi0_4, fi1_2, fi1_3, \
                         fi1_4, fk_4, fk_5, fk_6, fk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fi0_2[k]
                 - f_6 * fi1_2[k]
                 + pb_z[k] * fk_4[k];

        t_5[k] = f_7 * fi0_3[k]
                 - f_8 * fi1_3[k]
                 + pb_y[k] * fk_5[k];

        t_6[k] = f_3 * fi0_4[k]
                 - f_4 * fi1_4[k]
                 + pb_y[k] * fk_6[k];

        t_7[k] = f_7 * fi0_4[k]
                 - f_8 * fi1_4[k]
                 + pb_z[k] * fk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, fi0_5, fi0_7, fi0_8, fi1_5, fi1_6, \
                         fi1_7, fk_8, fk_9, fk_10, fk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * fi0_5[k]
                 - f_10 * fi1_5[k]
                 + pb_y[k] * fk_8[k];

        t_9[k] = f_5 * fi0_7[k]
                 - f_6 * fi1_6[k]
                 + pb_y[k] * fk_9[k];

        t_10[k] = f_3 * fi0_8[k]
                  - f_4 * fi1_7[k]
                  + pb_y[k] * fk_10[k];

        t_11[k] = f_9 * fi0_8[k]
                  - f_10 * fi1_7[k]
                  + pb_z[k] * fk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, fi0_9, fi0_11, fi0_12, fi1_8, fi1_9, fi1_10, \
                         fk_12, fk_13, fk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * fi0_9[k]
                  - f_12 * fi1_8[k]
                  + pb_y[k] * fk_12[k];

        t_13[k] = f_7 * fi0_11[k]
                  - f_8 * fi1_9[k]
                  + pb_y[k] * fk_13[k];

        t_14[k] = f_5 * fi0_12[k]
                  - f_6 * fi1_10[k]
                  + pb_y[k] * fk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, fi0_13, fi0_14, fi0_16, fi1_11, \
                         fi1_12, fi1_14, fk_15, fk_16, fk_17, fk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * fi0_13[k]
                  - f_4 * fi1_11[k]
                  + pb_y[k] * fk_15[k];

        t_16[k] = f_11 * fi0_13[k]
                  - f_12 * fi1_11[k]
                  + pb_z[k] * fk_16[k];

        t_17[k] = f_1 * fi0_14[k]
                  - f_2 * fi1_12[k]
                  + pb_y[k] * fk_17[k];

        t_18[k] = f_11 * fi0_16[k]
                  - f_12 * fi1_14[k]
                  + pb_y[k] * fk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, fi0_17, fi0_18, fi0_19, fi1_15, fi1_16, \
                         fi1_17, fk_19, fk_20, fk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_9 * fi0_17[k]
                  - f_10 * fi1_15[k]
                  + pb_y[k] * fk_19[k];

        t_20[k] = f_7 * fi0_18[k]
                  - f_8 * fi1_16[k]
                  + pb_y[k] * fk_20[k];

        t_21[k] = f_5 * fi0_19[k]
                  - f_6 * fi1_17[k]
                  + pb_y[k] * fk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_x, pa_y, pa_z, pb_y, pb_z, dl_0, \
                         dl_1, fi0_20, fi1_18, fk_22, fk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * fi0_20[k]
                  - f_4 * fi1_18[k]
                  + pb_y[k] * fk_22[k];

        t_23[k] = f_1 * fi0_20[k]
                  - f_2 * fi1_18[k]
                  + pb_z[k] * fk_23[k];

        t_24[k] = pa_y[k] * dl_0[k];

        t_25[k] = pa_z[k] * dl_0[k];

        t_26[k] = pa_x[k] * dl_1[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_x, dl_2, fi0_43, fi0_45, fi0_46, \
                         fi1_35, fi1_37, fi1_38, fk_28, fk_29, fk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * dl_2[k];

        t_28[k] = f_1 * fi0_43[k]
                  - f_2 * fi1_35[k]
                  + pb_x[k] * fk_28[k];

        t_29[k] = f_11 * fi0_45[k]
                  - f_12 * fi1_37[k]
                  + pb_x[k] * fk_29[k];

        t_30[k] = f_11 * fi0_46[k]
                  - f_12 * fi1_38[k]
                  + pb_x[k] * fk_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, fi0_47, fi0_49, fi0_50, fi1_39, fi1_40, \
                         fi1_41, fk_31, fk_32, fk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * fi0_47[k]
                  - f_10 * fi1_39[k]
                  + pb_x[k] * fk_31[k];

        t_32[k] = f_9 * fi0_49[k]
                  - f_10 * fi1_40[k]
                  + pb_x[k] * fk_32[k];

        t_33[k] = f_7 * fi0_50[k]
                  - f_8 * fi1_41[k]
                  + pb_x[k] * fk_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, fi0_52, fi0_53, fi0_54, fi1_42, fi1_43, \
                         fi1_44, fk_34, fk_35, fk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * fi0_52[k]
                  - f_8 * fi1_42[k]
                  + pb_x[k] * fk_34[k];

        t_35[k] = f_7 * fi0_53[k]
                  - f_8 * fi1_43[k]
                  + pb_x[k] * fk_35[k];

        t_36[k] = f_5 * fi0_54[k]
                  - f_6 * fi1_44[k]
                  + pb_x[k] * fk_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, fi0_55, fi0_56, fi0_57, fi1_45, fi1_46, \
                         fi1_47, fk_37, fk_38, fk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * fi0_55[k]
                  - f_6 * fi1_45[k]
                  + pb_x[k] * fk_37[k];

        t_38[k] = f_5 * fi0_56[k]
                  - f_6 * fi1_46[k]
                  + pb_x[k] * fk_38[k];

        t_39[k] = f_5 * fi0_57[k]
                  - f_6 * fi1_47[k]
                  + pb_x[k] * fk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, fi0_58, fi0_60, fi0_61, fi1_48, fi1_50, \
                         fi1_51, fk_40, fk_41, fk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_3 * fi0_58[k]
                  - f_4 * fi1_48[k]
                  + pb_x[k] * fk_40[k];

        t_41[k] = f_3 * fi0_60[k]
                  - f_4 * fi1_50[k]
                  + pb_x[k] * fk_41[k];

        t_42[k] = f_3 * fi0_61[k]
                  - f_4 * fi1_51[k]
                  + pb_x[k] * fk_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_y, dk_24, fi0_58, fi0_62, fi0_63, fi1_48, \
                         fi1_52, fi1_53, fk_43, fk_44, fk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * fi0_62[k]
                  - f_4 * fi1_52[k]
                  + pb_x[k] * fk_43[k];

        t_44[k] = f_3 * fi0_63[k]
                  - f_4 * fi1_53[k]
                  + pb_x[k] * fk_44[k];

        t_45[k] = f_0 * dk_24[k]
                  + f_1 * fi0_58[k]
                  - f_2 * fi1_48[k]
                  + pb_y[k] * fk_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_z, fi0_58, fi0_59, fi0_60, fi1_48, fi1_49, \
                         fi1_50, fk_46, fk_47, fk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * fi0_58[k]
                  - f_4 * fi1_48[k]
                  + pb_z[k] * fk_46[k];

        t_47[k] = f_5 * fi0_59[k]
                  - f_6 * fi1_49[k]
                  + pb_z[k] * fk_47[k];

        t_48[k] = f_7 * fi0_60[k]
                  - f_8 * fi1_50[k]
                  + pb_z[k] * fk_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_z, dl_1, fi0_61, fi0_62, fi0_63, \
                         fi1_51, fi1_52, fi1_53, fk_49, fk_50, fk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * fi0_61[k]
                  - f_10 * fi1_51[k]
                  + pb_z[k] * fk_49[k];

        t_50[k] = f_11 * fi0_62[k]
                  - f_12 * fi1_52[k]
                  + pb_z[k] * fk_50[k];

        t_51[k] = f_1 * fi0_63[k]
                  - f_2 * fi1_53[k]
                  + pb_z[k] * fk_51[k];

        t_52[k] = pa_z[k] * dl_1[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pb_x, dl_2, fi0_81, fi0_83, fi0_84, \
                         fi1_62, fi1_64, fi1_65, fk_54, fk_55, fk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_y[k] * dl_2[k];

        t_54[k] = f_1 * fi0_81[k]
                  - f_2 * fi1_62[k]
                  + pb_x[k] * fk_54[k];

        t_55[k] = f_11 * fi0_83[k]
                  - f_12 * fi1_64[k]
                  + pb_x[k] * fk_55[k];

        t_56[k] = f_11 * fi0_84[k]
                  - f_12 * fi1_65[k]
                  + pb_x[k] * fk_56[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, fi0_85, fi0_87, fi0_88, fi1_66, fi1_67, \
                         fi1_68, fk_57, fk_58, fk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_9 * fi0_85[k]
                  - f_10 * fi1_66[k]
                  + pb_x[k] * fk_57[k];

        t_58[k] = f_9 * fi0_87[k]
                  - f_10 * fi1_67[k]
                  + pb_x[k] * fk_58[k];

        t_59[k] = f_7 * fi0_88[k]
                  - f_8 * fi1_68[k]
                  + pb_x[k] * fk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_x, fi0_89, fi0_91, fi0_92, fi1_69, fi1_70, \
                         fi1_71, fk_60, fk_61, fk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_7 * fi0_89[k]
                  - f_8 * fi1_69[k]
                  + pb_x[k] * fk_60[k];

        t_61[k] = f_7 * fi0_91[k]
                  - f_8 * fi1_70[k]
                  + pb_x[k] * fk_61[k];

        t_62[k] = f_5 * fi0_92[k]
                  - f_6 * fi1_71[k]
                  + pb_x[k] * fk_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, fi0_93, fi0_94, fi0_95, fi1_72, fi1_73, \
                         fi1_74, fk_63, fk_64, fk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_5 * fi0_93[k]
                  - f_6 * fi1_72[k]
                  + pb_x[k] * fk_63[k];

        t_64[k] = f_5 * fi0_94[k]
                  - f_6 * fi1_73[k]
                  + pb_x[k] * fk_64[k];

        t_65[k] = f_5 * fi0_95[k]
                  - f_6 * fi1_74[k]
                  + pb_x[k] * fk_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, fi0_96, fi0_97, fi0_98, fi1_75, fi1_76, \
                         fi1_77, fk_66, fk_67, fk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * fi0_96[k]
                  - f_4 * fi1_75[k]
                  + pb_x[k] * fk_66[k];

        t_67[k] = f_3 * fi0_97[k]
                  - f_4 * fi1_76[k]
                  + pb_x[k] * fk_67[k];

        t_68[k] = f_3 * fi0_98[k]
                  - f_4 * fi1_77[k]
                  + pb_x[k] * fk_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_y, fi0_96, fi0_99, fi0_101, fi1_75, \
                         fi1_78, fi1_80, fk_69, fk_70, fk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * fi0_99[k]
                  - f_4 * fi1_78[k]
                  + pb_x[k] * fk_69[k];

        t_70[k] = f_3 * fi0_101[k]
                  - f_4 * fi1_80[k]
                  + pb_x[k] * fk_70[k];

        t_71[k] = f_1 * fi0_96[k]
                  - f_2 * fi1_75[k]
                  + pb_y[k] * fk_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, fi0_97, fi0_98, fi0_99, fi1_76, fi1_77, \
                         fi1_78, fk_72, fk_73, fk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * fi0_97[k]
                  - f_12 * fi1_76[k]
                  + pb_y[k] * fk_72[k];

        t_73[k] = f_9 * fi0_98[k]
                  - f_10 * fi1_77[k]
                  + pb_y[k] * fk_73[k];

        t_74[k] = f_7 * fi0_99[k]
                  - f_8 * fi1_78[k]
                  + pb_y[k] * fk_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_y, pb_z, dk_44, fi0_100, fi0_101, fi1_79, \
                         fi1_80, fk_75, fk_76, fk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * fi0_100[k]
                  - f_6 * fi1_79[k]
                  + pb_y[k] * fk_75[k];

        t_76[k] = f_3 * fi0_101[k]
                  - f_4 * fi1_80[k]
                  + pb_y[k] * fk_76[k];

        t_77[k] = f_0 * dk_44[k]
                  + f_1 * fi0_101[k]
                  - f_2 * fi1_80[k]
                  + pb_z[k] * fk_77[k];
    }
}

auto
compute_prim_fl_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_8 = buffer.data(fi0 + 8);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_13 = buffer.data(fi1 + 13);
    const auto *fi1_20 = buffer.data(fi1 + 20);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_23 = buffer.data(fk + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_5, fi1_13, \
                         fk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_5[k]
                 - f_2 * fi1_13[k]
                 + pb_y[k] * fk_15[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_8, fi1_20, fk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_8[k]
                 - f_2 * fi1_20[k]
                 + pb_z[k] * fk_23[k];
    }
}

auto
compute_prim_fl_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_13 = buffer.data(fi0 + 13);
    const auto *fi0_20 = buffer.data(fi0 + 20);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_40 = buffer.data(fi1 + 40);
    const auto *fi1_68 = buffer.data(fi1 + 68);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_65 = buffer.data(fk + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_13, fi1_40, \
                         fk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_13[k]
                 - f_2 * fi1_40[k]
                 + pb_y[k] * fk_38[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_20, fi1_68, fk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_20[k]
                 - f_2 * fi1_68[k]
                 + pb_z[k] * fk_65[k];
    }
}

auto
compute_prim_fl_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_40 = buffer.data(fi0 + 40);
    const auto *fi0_68 = buffer.data(fi0 + 68);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_36 = buffer.data(fi1 + 36);
    const auto *fi1_62 = buffer.data(fi1 + 62);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_8 = buffer.data(fk + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_40, fi1_36, \
                         fk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_40[k]
                 - f_2 * fi1_36[k]
                 + pb_y[k] * fk_5[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_68, fi1_62, fk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_68[k]
                 - f_2 * fi1_62[k]
                 + pb_z[k] * fk_8[k];
    }
}

auto
compute_prim_fl_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_8 = buffer.data(fi0 + 8);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_32 = buffer.data(fi1 + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_35 = buffer.data(fk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_1, dl_1, dl_2, fi0_5, fi1_17, \
                         fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_1[k]
                 + f_1 * fi0_5[k]
                 - f_2 * fi1_17[k]
                 + pb_y[k] * fk_19[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_2, fi0_8, fi1_32, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_2[k]
                 + f_1 * fi0_8[k]
                 - f_2 * fi1_32[k]
                 + pb_z[k] * fk_35[k];
    }
}

auto
compute_prim_fl_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk, const size_t dl,
                                     const size_t fi0, const size_t fi1, const size_t fk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_17 = buffer.data(fi0 + 17);
    const auto *fi0_32 = buffer.data(fi0 + 32);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_32 = buffer.data(fi1 + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_35 = buffer.data(fk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_17, fi1_17, \
                         fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_17[k]
                 - f_2 * fi1_17[k]
                 + pb_y[k] * fk_19[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_32, fi1_32, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_32[k]
                 - f_2 * fi1_32[k]
                 + pb_z[k] * fk_35[k];
    }
}

auto
compute_prim_fl_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk, const size_t dl,
                                      const size_t fi0, const size_t fi1, const size_t fk,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_17 = buffer.data(fi0 + 17);
    const auto *fi0_32 = buffer.data(fi0 + 32);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_20 = buffer.data(fi1 + 20);
    const auto *fi1_38 = buffer.data(fi1 + 38);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_35 = buffer.data(fk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_17, fi1_20, \
                         fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_17[k]
                 - f_2 * fi1_20[k]
                 + pb_y[k] * fk_19[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_32, fi1_38, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_32[k]
                 - f_2 * fi1_38[k]
                 + pb_z[k] * fk_35[k];
    }
}

auto
compute_prim_fl_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk, const size_t dl,
                                      const size_t fi0, const size_t fi1, const size_t fk,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_20 = buffer.data(fi0 + 20);
    const auto *fi0_38 = buffer.data(fi0 + 38);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_32 = buffer.data(fi1 + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_8 = buffer.data(fk + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_20, fi1_17, \
                         fk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_20[k]
                 - f_2 * fi1_17[k]
                 + pb_y[k] * fk_5[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_38, fi1_32, fk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_38[k]
                 - f_2 * fi1_32[k]
                 + pb_z[k] * fk_8[k];
    }
}

auto
compute_prim_fl_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk, const size_t dl,
                                      const size_t fi0, const size_t fi1, const size_t fk,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_13 = buffer.data(fi0 + 13);
    const auto *fi0_20 = buffer.data(fi0 + 20);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_32 = buffer.data(fi1 + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_35 = buffer.data(fk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_13, fi1_17, \
                         fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_13[k]
                 - f_2 * fi1_17[k]
                 + pb_y[k] * fk_19[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_20, fi1_32, fk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_20[k]
                 - f_2 * fi1_32[k]
                 + pb_z[k] * fk_35[k];
    }
}

auto
compute_prim_fl_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dk, const size_t dl,
                                      const size_t fi0, const size_t fi1, const size_t fk,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_17 = buffer.data(fi0 + 17);
    const auto *fi0_32 = buffer.data(fi0 + 32);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_17 = buffer.data(fi1 + 17);
    const auto *fi1_32 = buffer.data(fi1 + 32);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_8 = buffer.data(fk + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dk_0, dl_0, dl_1, \
                         dl_2, fi0_0, fi1_0, fk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_0[k]
                 + f_1 * fi0_0[k]
                 - f_2 * fi1_0[k]
                 + pb_x[k] * fk_0[k];

        t_1[k] = pa_y[k] * dl_0[k];

        t_2[k] = pa_z[k] * dl_0[k];

        t_3[k] = pa_x[k] * dl_1[k];

        t_4[k] = pa_x[k] * dl_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dk_5, dl_1, dl_2, fi0_17, fi1_17, \
                         fk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dk_5[k]
                 + f_1 * fi0_17[k]
                 - f_2 * fi1_17[k]
                 + pb_y[k] * fk_5[k];

        t_6[k] = pa_z[k] * dl_1[k];

        t_7[k] = pa_y[k] * dl_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dk_14, fi0_32, fi1_32, fk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_14[k]
                 + f_1 * fi0_32[k]
                 - f_2 * fi1_32[k]
                 + pb_z[k] * fk_8[k];
    }
}

}  // namespace simdt2ceri
