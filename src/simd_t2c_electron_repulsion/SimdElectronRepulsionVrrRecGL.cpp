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


#include "SimdElectronRepulsionVrrRecGL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_gl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

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
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);
    const auto *fl_9 = buffer.data(fl + 9);
    const auto *fl_10 = buffer.data(fl + 10);
    const auto *fl_11 = buffer.data(fl + 11);
    const auto *fl_12 = buffer.data(fl + 12);
    const auto *fl_13 = buffer.data(fl + 13);
    const auto *fl_14 = buffer.data(fl + 14);
    const auto *fl_15 = buffer.data(fl + 15);
    const auto *fl_16 = buffer.data(fl + 16);
    const auto *fl_17 = buffer.data(fl + 17);
    const auto *fl_18 = buffer.data(fl + 18);
    const auto *fl_19 = buffer.data(fl + 19);
    const auto *fl_20 = buffer.data(fl + 20);
    const auto *fl_21 = buffer.data(fl + 21);
    const auto *fl_22 = buffer.data(fl + 22);
    const auto *fl_23 = buffer.data(fl + 23);
    const auto *fl_24 = buffer.data(fl + 24);
    const auto *fl_25 = buffer.data(fl + 25);
    const auto *fl_26 = buffer.data(fl + 26);
    const auto *fl_27 = buffer.data(fl + 27);
    const auto *fl_28 = buffer.data(fl + 28);
    const auto *fl_29 = buffer.data(fl + 29);
    const auto *fl_30 = buffer.data(fl + 30);
    const auto *fl_31 = buffer.data(fl + 31);
    const auto *fl_32 = buffer.data(fl + 32);
    const auto *fl_33 = buffer.data(fl + 33);
    const auto *fl_34 = buffer.data(fl + 34);
    const auto *fl_35 = buffer.data(fl + 35);
    const auto *fl_36 = buffer.data(fl + 36);
    const auto *fl_37 = buffer.data(fl + 37);
    const auto *fl_38 = buffer.data(fl + 38);
    const auto *fl_39 = buffer.data(fl + 39);
    const auto *fl_40 = buffer.data(fl + 40);
    const auto *fl_41 = buffer.data(fl + 41);
    const auto *fl_42 = buffer.data(fl + 42);
    const auto *fl_43 = buffer.data(fl + 43);
    const auto *fl_44 = buffer.data(fl + 44);
    const auto *fl_45 = buffer.data(fl + 45);
    const auto *fl_46 = buffer.data(fl + 46);
    const auto *fl_47 = buffer.data(fl + 47);
    const auto *fl_48 = buffer.data(fl + 48);
    const auto *fl_49 = buffer.data(fl + 49);
    const auto *fl_50 = buffer.data(fl + 50);
    const auto *fl_51 = buffer.data(fl + 51);
    const auto *fl_52 = buffer.data(fl + 52);
    const auto *fl_53 = buffer.data(fl + 53);
    const auto *fl_54 = buffer.data(fl + 54);
    const auto *fl_55 = buffer.data(fl + 55);
    const auto *fl_56 = buffer.data(fl + 56);
    const auto *fl_57 = buffer.data(fl + 57);
    const auto *fl_58 = buffer.data(fl + 58);
    const auto *fl_59 = buffer.data(fl + 59);
    const auto *fl_60 = buffer.data(fl + 60);
    const auto *fl_61 = buffer.data(fl + 61);
    const auto *fl_62 = buffer.data(fl + 62);
    const auto *fl_63 = buffer.data(fl + 63);
    const auto *fl_64 = buffer.data(fl + 64);
    const auto *fl_65 = buffer.data(fl + 65);
    const auto *fl_66 = buffer.data(fl + 66);
    const auto *fl_67 = buffer.data(fl + 67);
    const auto *fl_68 = buffer.data(fl + 68);
    const auto *fl_69 = buffer.data(fl + 69);
    const auto *fl_70 = buffer.data(fl + 70);
    const auto *fl_71 = buffer.data(fl + 71);
    const auto *fl_72 = buffer.data(fl + 72);
    const auto *fl_73 = buffer.data(fl + 73);
    const auto *fl_74 = buffer.data(fl + 74);
    const auto *fl_75 = buffer.data(fl + 75);
    const auto *fl_76 = buffer.data(fl + 76);
    const auto *fl_77 = buffer.data(fl + 77);
    const auto *fl_78 = buffer.data(fl + 78);
    const auto *fl_79 = buffer.data(fl + 79);
    const auto *fl_80 = buffer.data(fl + 80);
    const auto *fl_81 = buffer.data(fl + 81);
    const auto *fl_82 = buffer.data(fl + 82);
    const auto *fl_83 = buffer.data(fl + 83);
    const auto *fl_84 = buffer.data(fl + 84);
    const auto *fl_85 = buffer.data(fl + 85);
    const auto *fl_86 = buffer.data(fl + 86);
    const auto *fl_87 = buffer.data(fl + 87);
    const auto *fl_88 = buffer.data(fl + 88);
    const auto *fl_89 = buffer.data(fl + 89);
    const auto *fl_90 = buffer.data(fl + 90);
    const auto *fl_91 = buffer.data(fl + 91);
    const auto *fl_92 = buffer.data(fl + 92);
    const auto *fl_93 = buffer.data(fl + 93);
    const auto *fl_94 = buffer.data(fl + 94);
    const auto *fl_95 = buffer.data(fl + 95);
    const auto *fl_96 = buffer.data(fl + 96);
    const auto *fl_97 = buffer.data(fl + 97);
    const auto *fl_98 = buffer.data(fl + 98);
    const auto *fl_99 = buffer.data(fl + 99);
    const auto *fl_100 = buffer.data(fl + 100);
    const auto *fl_101 = buffer.data(fl + 101);
    const auto *fl_102 = buffer.data(fl + 102);
    const auto *fl_103 = buffer.data(fl + 103);
    const auto *fl_104 = buffer.data(fl + 104);
    const auto *fl_105 = buffer.data(fl + 105);
    const auto *fl_106 = buffer.data(fl + 106);
    const auto *fl_107 = buffer.data(fl + 107);
    const auto *fl_108 = buffer.data(fl + 108);
    const auto *fl_109 = buffer.data(fl + 109);
    const auto *fl_110 = buffer.data(fl + 110);
    const auto *fl_111 = buffer.data(fl + 111);
    const auto *fl_112 = buffer.data(fl + 112);
    const auto *fl_113 = buffer.data(fl + 113);
    const auto *fl_114 = buffer.data(fl + 114);
    const auto *fl_115 = buffer.data(fl + 115);
    const auto *fl_116 = buffer.data(fl + 116);
    const auto *fl_117 = buffer.data(fl + 117);
    const auto *fl_118 = buffer.data(fl + 118);
    const auto *fl_119 = buffer.data(fl + 119);
    const auto *fl_120 = buffer.data(fl + 120);
    const auto *fl_121 = buffer.data(fl + 121);
    const auto *fl_122 = buffer.data(fl + 122);
    const auto *fl_123 = buffer.data(fl + 123);
    const auto *fl_124 = buffer.data(fl + 124);
    const auto *fl_125 = buffer.data(fl + 125);
    const auto *fl_126 = buffer.data(fl + 126);
    const auto *fl_127 = buffer.data(fl + 127);
    const auto *fl_128 = buffer.data(fl + 128);
    const auto *fl_129 = buffer.data(fl + 129);
    const auto *fl_130 = buffer.data(fl + 130);
    const auto *fl_131 = buffer.data(fl + 131);
    const auto *fl_132 = buffer.data(fl + 132);
    const auto *fl_133 = buffer.data(fl + 133);
    const auto *fl_134 = buffer.data(fl + 134);
    const auto *fl_135 = buffer.data(fl + 135);
    const auto *fl_136 = buffer.data(fl + 136);
    const auto *fl_137 = buffer.data(fl + 137);
    const auto *fl_138 = buffer.data(fl + 138);
    const auto *fl_139 = buffer.data(fl + 139);
    const auto *fl_140 = buffer.data(fl + 140);
    const auto *fl_141 = buffer.data(fl + 141);
    const auto *fl_142 = buffer.data(fl + 142);
    const auto *fl_143 = buffer.data(fl + 143);
    const auto *fl_144 = buffer.data(fl + 144);
    const auto *fl_145 = buffer.data(fl + 145);
    const auto *fl_146 = buffer.data(fl + 146);
    const auto *fl_147 = buffer.data(fl + 147);
    const auto *fl_148 = buffer.data(fl + 148);
    const auto *fl_149 = buffer.data(fl + 149);
    const auto *fl_150 = buffer.data(fl + 150);
    const auto *fl_151 = buffer.data(fl + 151);
    const auto *fl_152 = buffer.data(fl + 152);
    const auto *fl_153 = buffer.data(fl + 153);
    const auto *fl_154 = buffer.data(fl + 154);
    const auto *fl_155 = buffer.data(fl + 155);
    const auto *fl_156 = buffer.data(fl + 156);
    const auto *fl_157 = buffer.data(fl + 157);
    const auto *fl_158 = buffer.data(fl + 158);
    const auto *fl_159 = buffer.data(fl + 159);
    const auto *fl_160 = buffer.data(fl + 160);
    const auto *fl_161 = buffer.data(fl + 161);
    const auto *fl_162 = buffer.data(fl + 162);
    const auto *fl_163 = buffer.data(fl + 163);
    const auto *fl_164 = buffer.data(fl + 164);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_1 = buffer.data(gi0 + 1);
    const auto *gi0_2 = buffer.data(gi0 + 2);
    const auto *gi0_3 = buffer.data(gi0 + 3);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_9 = buffer.data(gi0 + 9);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_15 = buffer.data(gi0 + 15);
    const auto *gi0_16 = buffer.data(gi0 + 16);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_18 = buffer.data(gi0 + 18);
    const auto *gi0_19 = buffer.data(gi0 + 19);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_25 = buffer.data(gi0 + 25);
    const auto *gi0_26 = buffer.data(gi0 + 26);
    const auto *gi0_27 = buffer.data(gi0 + 27);
    const auto *gi0_28 = buffer.data(gi0 + 28);
    const auto *gi0_29 = buffer.data(gi0 + 29);
    const auto *gi0_30 = buffer.data(gi0 + 30);
    const auto *gi0_31 = buffer.data(gi0 + 31);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_33 = buffer.data(gi0 + 33);
    const auto *gi0_34 = buffer.data(gi0 + 34);
    const auto *gi0_35 = buffer.data(gi0 + 35);
    const auto *gi0_36 = buffer.data(gi0 + 36);
    const auto *gi0_37 = buffer.data(gi0 + 37);
    const auto *gi0_38 = buffer.data(gi0 + 38);
    const auto *gi0_39 = buffer.data(gi0 + 39);
    const auto *gi0_40 = buffer.data(gi0 + 40);
    const auto *gi0_41 = buffer.data(gi0 + 41);
    const auto *gi0_42 = buffer.data(gi0 + 42);
    const auto *gi0_43 = buffer.data(gi0 + 43);
    const auto *gi0_44 = buffer.data(gi0 + 44);
    const auto *gi0_45 = buffer.data(gi0 + 45);
    const auto *gi0_46 = buffer.data(gi0 + 46);
    const auto *gi0_47 = buffer.data(gi0 + 47);
    const auto *gi0_48 = buffer.data(gi0 + 48);
    const auto *gi0_49 = buffer.data(gi0 + 49);
    const auto *gi0_50 = buffer.data(gi0 + 50);
    const auto *gi0_51 = buffer.data(gi0 + 51);
    const auto *gi0_52 = buffer.data(gi0 + 52);
    const auto *gi0_53 = buffer.data(gi0 + 53);
    const auto *gi0_54 = buffer.data(gi0 + 54);
    const auto *gi0_55 = buffer.data(gi0 + 55);
    const auto *gi0_56 = buffer.data(gi0 + 56);
    const auto *gi0_57 = buffer.data(gi0 + 57);
    const auto *gi0_58 = buffer.data(gi0 + 58);
    const auto *gi0_59 = buffer.data(gi0 + 59);
    const auto *gi0_60 = buffer.data(gi0 + 60);
    const auto *gi0_61 = buffer.data(gi0 + 61);
    const auto *gi0_62 = buffer.data(gi0 + 62);
    const auto *gi0_63 = buffer.data(gi0 + 63);
    const auto *gi0_64 = buffer.data(gi0 + 64);
    const auto *gi0_65 = buffer.data(gi0 + 65);
    const auto *gi0_66 = buffer.data(gi0 + 66);
    const auto *gi0_67 = buffer.data(gi0 + 67);
    const auto *gi0_68 = buffer.data(gi0 + 68);
    const auto *gi0_69 = buffer.data(gi0 + 69);
    const auto *gi0_70 = buffer.data(gi0 + 70);
    const auto *gi0_71 = buffer.data(gi0 + 71);
    const auto *gi0_72 = buffer.data(gi0 + 72);
    const auto *gi0_73 = buffer.data(gi0 + 73);
    const auto *gi0_74 = buffer.data(gi0 + 74);
    const auto *gi0_75 = buffer.data(gi0 + 75);
    const auto *gi0_76 = buffer.data(gi0 + 76);
    const auto *gi0_77 = buffer.data(gi0 + 77);
    const auto *gi0_78 = buffer.data(gi0 + 78);
    const auto *gi0_79 = buffer.data(gi0 + 79);
    const auto *gi0_80 = buffer.data(gi0 + 80);
    const auto *gi0_81 = buffer.data(gi0 + 81);
    const auto *gi0_82 = buffer.data(gi0 + 82);
    const auto *gi0_83 = buffer.data(gi0 + 83);
    const auto *gi0_84 = buffer.data(gi0 + 84);
    const auto *gi0_85 = buffer.data(gi0 + 85);
    const auto *gi0_86 = buffer.data(gi0 + 86);
    const auto *gi0_87 = buffer.data(gi0 + 87);
    const auto *gi0_88 = buffer.data(gi0 + 88);
    const auto *gi0_89 = buffer.data(gi0 + 89);
    const auto *gi0_90 = buffer.data(gi0 + 90);
    const auto *gi0_91 = buffer.data(gi0 + 91);
    const auto *gi0_92 = buffer.data(gi0 + 92);
    const auto *gi0_93 = buffer.data(gi0 + 93);
    const auto *gi0_94 = buffer.data(gi0 + 94);
    const auto *gi0_95 = buffer.data(gi0 + 95);
    const auto *gi0_96 = buffer.data(gi0 + 96);
    const auto *gi0_97 = buffer.data(gi0 + 97);
    const auto *gi0_98 = buffer.data(gi0 + 98);
    const auto *gi0_99 = buffer.data(gi0 + 99);
    const auto *gi0_100 = buffer.data(gi0 + 100);
    const auto *gi0_101 = buffer.data(gi0 + 101);
    const auto *gi0_102 = buffer.data(gi0 + 102);
    const auto *gi0_103 = buffer.data(gi0 + 103);
    const auto *gi0_104 = buffer.data(gi0 + 104);
    const auto *gi0_105 = buffer.data(gi0 + 105);
    const auto *gi0_106 = buffer.data(gi0 + 106);
    const auto *gi0_107 = buffer.data(gi0 + 107);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_1 = buffer.data(gi1 + 1);
    const auto *gi1_2 = buffer.data(gi1 + 2);
    const auto *gi1_3 = buffer.data(gi1 + 3);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_9 = buffer.data(gi1 + 9);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_15 = buffer.data(gi1 + 15);
    const auto *gi1_16 = buffer.data(gi1 + 16);
    const auto *gi1_17 = buffer.data(gi1 + 17);
    const auto *gi1_18 = buffer.data(gi1 + 18);
    const auto *gi1_19 = buffer.data(gi1 + 19);
    const auto *gi1_20 = buffer.data(gi1 + 20);
    const auto *gi1_21 = buffer.data(gi1 + 21);
    const auto *gi1_22 = buffer.data(gi1 + 22);
    const auto *gi1_23 = buffer.data(gi1 + 23);
    const auto *gi1_24 = buffer.data(gi1 + 24);
    const auto *gi1_25 = buffer.data(gi1 + 25);
    const auto *gi1_26 = buffer.data(gi1 + 26);
    const auto *gi1_27 = buffer.data(gi1 + 27);
    const auto *gi1_28 = buffer.data(gi1 + 28);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_30 = buffer.data(gi1 + 30);
    const auto *gi1_31 = buffer.data(gi1 + 31);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_33 = buffer.data(gi1 + 33);
    const auto *gi1_34 = buffer.data(gi1 + 34);
    const auto *gi1_35 = buffer.data(gi1 + 35);
    const auto *gi1_36 = buffer.data(gi1 + 36);
    const auto *gi1_37 = buffer.data(gi1 + 37);
    const auto *gi1_38 = buffer.data(gi1 + 38);
    const auto *gi1_39 = buffer.data(gi1 + 39);
    const auto *gi1_40 = buffer.data(gi1 + 40);
    const auto *gi1_41 = buffer.data(gi1 + 41);
    const auto *gi1_42 = buffer.data(gi1 + 42);
    const auto *gi1_43 = buffer.data(gi1 + 43);
    const auto *gi1_44 = buffer.data(gi1 + 44);
    const auto *gi1_45 = buffer.data(gi1 + 45);
    const auto *gi1_46 = buffer.data(gi1 + 46);
    const auto *gi1_47 = buffer.data(gi1 + 47);
    const auto *gi1_48 = buffer.data(gi1 + 48);
    const auto *gi1_49 = buffer.data(gi1 + 49);
    const auto *gi1_50 = buffer.data(gi1 + 50);
    const auto *gi1_51 = buffer.data(gi1 + 51);
    const auto *gi1_52 = buffer.data(gi1 + 52);
    const auto *gi1_53 = buffer.data(gi1 + 53);
    const auto *gi1_54 = buffer.data(gi1 + 54);
    const auto *gi1_55 = buffer.data(gi1 + 55);
    const auto *gi1_56 = buffer.data(gi1 + 56);
    const auto *gi1_57 = buffer.data(gi1 + 57);
    const auto *gi1_58 = buffer.data(gi1 + 58);
    const auto *gi1_59 = buffer.data(gi1 + 59);
    const auto *gi1_60 = buffer.data(gi1 + 60);
    const auto *gi1_61 = buffer.data(gi1 + 61);
    const auto *gi1_62 = buffer.data(gi1 + 62);
    const auto *gi1_63 = buffer.data(gi1 + 63);
    const auto *gi1_64 = buffer.data(gi1 + 64);
    const auto *gi1_65 = buffer.data(gi1 + 65);
    const auto *gi1_66 = buffer.data(gi1 + 66);
    const auto *gi1_67 = buffer.data(gi1 + 67);
    const auto *gi1_68 = buffer.data(gi1 + 68);
    const auto *gi1_69 = buffer.data(gi1 + 69);
    const auto *gi1_70 = buffer.data(gi1 + 70);
    const auto *gi1_71 = buffer.data(gi1 + 71);
    const auto *gi1_72 = buffer.data(gi1 + 72);
    const auto *gi1_73 = buffer.data(gi1 + 73);
    const auto *gi1_74 = buffer.data(gi1 + 74);
    const auto *gi1_75 = buffer.data(gi1 + 75);
    const auto *gi1_76 = buffer.data(gi1 + 76);
    const auto *gi1_77 = buffer.data(gi1 + 77);
    const auto *gi1_78 = buffer.data(gi1 + 78);
    const auto *gi1_79 = buffer.data(gi1 + 79);
    const auto *gi1_80 = buffer.data(gi1 + 80);
    const auto *gi1_81 = buffer.data(gi1 + 81);
    const auto *gi1_82 = buffer.data(gi1 + 82);
    const auto *gi1_83 = buffer.data(gi1 + 83);
    const auto *gi1_84 = buffer.data(gi1 + 84);
    const auto *gi1_85 = buffer.data(gi1 + 85);
    const auto *gi1_86 = buffer.data(gi1 + 86);
    const auto *gi1_87 = buffer.data(gi1 + 87);
    const auto *gi1_88 = buffer.data(gi1 + 88);
    const auto *gi1_89 = buffer.data(gi1 + 89);
    const auto *gi1_90 = buffer.data(gi1 + 90);
    const auto *gi1_91 = buffer.data(gi1 + 91);
    const auto *gi1_92 = buffer.data(gi1 + 92);
    const auto *gi1_93 = buffer.data(gi1 + 93);
    const auto *gi1_94 = buffer.data(gi1 + 94);
    const auto *gi1_95 = buffer.data(gi1 + 95);
    const auto *gi1_96 = buffer.data(gi1 + 96);
    const auto *gi1_97 = buffer.data(gi1 + 97);
    const auto *gi1_98 = buffer.data(gi1 + 98);
    const auto *gi1_99 = buffer.data(gi1 + 99);
    const auto *gi1_100 = buffer.data(gi1 + 100);
    const auto *gi1_101 = buffer.data(gi1 + 101);
    const auto *gi1_102 = buffer.data(gi1 + 102);
    const auto *gi1_103 = buffer.data(gi1 + 103);
    const auto *gi1_104 = buffer.data(gi1 + 104);
    const auto *gi1_105 = buffer.data(gi1 + 105);
    const auto *gi1_106 = buffer.data(gi1 + 106);
    const auto *gi1_107 = buffer.data(gi1 + 107);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_1 = buffer.data(gk + 1);
    const auto *gk_2 = buffer.data(gk + 2);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_16 = buffer.data(gk + 16);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_47 = buffer.data(gk + 47);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_49 = buffer.data(gk + 49);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_53 = buffer.data(gk + 53);
    const auto *gk_54 = buffer.data(gk + 54);
    const auto *gk_55 = buffer.data(gk + 55);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_58 = buffer.data(gk + 58);
    const auto *gk_59 = buffer.data(gk + 59);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_65 = buffer.data(gk + 65);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_73 = buffer.data(gk + 73);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_83 = buffer.data(gk + 83);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_85 = buffer.data(gk + 85);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_88 = buffer.data(gk + 88);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_91 = buffer.data(gk + 91);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_93 = buffer.data(gk + 93);
    const auto *gk_94 = buffer.data(gk + 94);
    const auto *gk_95 = buffer.data(gk + 95);
    const auto *gk_96 = buffer.data(gk + 96);
    const auto *gk_97 = buffer.data(gk + 97);
    const auto *gk_98 = buffer.data(gk + 98);
    const auto *gk_99 = buffer.data(gk + 99);
    const auto *gk_100 = buffer.data(gk + 100);
    const auto *gk_101 = buffer.data(gk + 101);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_106 = buffer.data(gk + 106);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_110 = buffer.data(gk + 110);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_112 = buffer.data(gk + 112);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_115 = buffer.data(gk + 115);
    const auto *gk_116 = buffer.data(gk + 116);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_119 = buffer.data(gk + 119);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_121 = buffer.data(gk + 121);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_124 = buffer.data(gk + 124);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_126 = buffer.data(gk + 126);
    const auto *gk_127 = buffer.data(gk + 127);
    const auto *gk_128 = buffer.data(gk + 128);
    const auto *gk_129 = buffer.data(gk + 129);
    const auto *gk_130 = buffer.data(gk + 130);
    const auto *gk_131 = buffer.data(gk + 131);
    const auto *gk_132 = buffer.data(gk + 132);
    const auto *gk_133 = buffer.data(gk + 133);
    const auto *gk_134 = buffer.data(gk + 134);
    const auto *gk_135 = buffer.data(gk + 135);
    const auto *gk_136 = buffer.data(gk + 136);
    const auto *gk_137 = buffer.data(gk + 137);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_148 = buffer.data(gk + 148);
    const auto *gk_149 = buffer.data(gk + 149);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_151 = buffer.data(gk + 151);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_160 = buffer.data(gk + 160);
    const auto *gk_161 = buffer.data(gk + 161);
    const auto *gk_162 = buffer.data(gk + 162);
    const auto *gk_163 = buffer.data(gk + 163);
    const auto *gk_164 = buffer.data(gk + 164);
    const auto *gk_165 = buffer.data(gk + 165);
    const auto *gk_166 = buffer.data(gk + 166);
    const auto *gk_167 = buffer.data(gk + 167);
    const auto *gk_168 = buffer.data(gk + 168);
    const auto *gk_169 = buffer.data(gk + 169);
    const auto *gk_170 = buffer.data(gk + 170);
    const auto *gk_171 = buffer.data(gk + 171);
    const auto *gk_172 = buffer.data(gk + 172);
    const auto *gk_173 = buffer.data(gk + 173);
    const auto *gk_174 = buffer.data(gk + 174);
    const auto *gk_175 = buffer.data(gk + 175);
    const auto *gk_176 = buffer.data(gk + 176);
    const auto *gk_177 = buffer.data(gk + 177);
    const auto *gk_178 = buffer.data(gk + 178);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_181 = buffer.data(gk + 181);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);
    const auto *gk_184 = buffer.data(gk + 184);
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_187 = buffer.data(gk + 187);
    const auto *gk_188 = buffer.data(gk + 188);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_191 = buffer.data(gk + 191);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_193 = buffer.data(gk + 193);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_196 = buffer.data(gk + 196);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_199 = buffer.data(gk + 199);
    const auto *gk_200 = buffer.data(gk + 200);
    const auto *gk_201 = buffer.data(gk + 201);
    const auto *gk_202 = buffer.data(gk + 202);
    const auto *gk_203 = buffer.data(gk + 203);
    const auto *gk_204 = buffer.data(gk + 204);
    const auto *gk_205 = buffer.data(gk + 205);
    const auto *gk_206 = buffer.data(gk + 206);
    const auto *gk_207 = buffer.data(gk + 207);
    const auto *gk_208 = buffer.data(gk + 208);
    const auto *gk_209 = buffer.data(gk + 209);
    const auto *gk_210 = buffer.data(gk + 210);
    const auto *gk_211 = buffer.data(gk + 211);
    const auto *gk_212 = buffer.data(gk + 212);
    const auto *gk_213 = buffer.data(gk + 213);
    const auto *gk_214 = buffer.data(gk + 214);
    const auto *gk_215 = buffer.data(gk + 215);
    const auto *gk_216 = buffer.data(gk + 216);
    const auto *gk_217 = buffer.data(gk + 217);
    const auto *gk_218 = buffer.data(gk + 218);
    const auto *gk_219 = buffer.data(gk + 219);
    const auto *gk_220 = buffer.data(gk + 220);
    const auto *gk_221 = buffer.data(gk + 221);
    const auto *gk_222 = buffer.data(gk + 222);
    const auto *gk_223 = buffer.data(gk + 223);
    const auto *gk_224 = buffer.data(gk + 224);
    const auto *gk_225 = buffer.data(gk + 225);
    const auto *gk_226 = buffer.data(gk + 226);
    const auto *gk_227 = buffer.data(gk + 227);
    const auto *gk_228 = buffer.data(gk + 228);
    const auto *gk_229 = buffer.data(gk + 229);
    const auto *gk_230 = buffer.data(gk + 230);
    const auto *gk_231 = buffer.data(gk + 231);
    const auto *gk_232 = buffer.data(gk + 232);
    const auto *gk_233 = buffer.data(gk + 233);
    const auto *gk_234 = buffer.data(gk + 234);
    const auto *gk_235 = buffer.data(gk + 235);
    const auto *gk_236 = buffer.data(gk + 236);
    const auto *gk_237 = buffer.data(gk + 237);
    const auto *gk_238 = buffer.data(gk + 238);
    const auto *gk_239 = buffer.data(gk + 239);
    const auto *gk_240 = buffer.data(gk + 240);
    const auto *gk_241 = buffer.data(gk + 241);
    const auto *gk_242 = buffer.data(gk + 242);
    const auto *gk_243 = buffer.data(gk + 243);
    const auto *gk_244 = buffer.data(gk + 244);
    const auto *gk_245 = buffer.data(gk + 245);
    const auto *gk_246 = buffer.data(gk + 246);
    const auto *gk_247 = buffer.data(gk + 247);
    const auto *gk_248 = buffer.data(gk + 248);
    const auto *gk_249 = buffer.data(gk + 249);
    const auto *gk_250 = buffer.data(gk + 250);
    const auto *gk_251 = buffer.data(gk + 251);
    const auto *gk_252 = buffer.data(gk + 252);
    const auto *gk_253 = buffer.data(gk + 253);
    const auto *gk_254 = buffer.data(gk + 254);
    const auto *gk_255 = buffer.data(gk + 255);
    const auto *gk_256 = buffer.data(gk + 256);
    const auto *gk_257 = buffer.data(gk + 257);
    const auto *gk_258 = buffer.data(gk + 258);
    const auto *gk_259 = buffer.data(gk + 259);
    const auto *gk_260 = buffer.data(gk + 260);
    const auto *gk_261 = buffer.data(gk + 261);
    const auto *gk_262 = buffer.data(gk + 262);
    const auto *gk_263 = buffer.data(gk + 263);
    const auto *gk_264 = buffer.data(gk + 264);
    const auto *gk_265 = buffer.data(gk + 265);
    const auto *gk_266 = buffer.data(gk + 266);
    const auto *gk_267 = buffer.data(gk + 267);
    const auto *gk_268 = buffer.data(gk + 268);
    const auto *gk_269 = buffer.data(gk + 269);
    const auto *gk_270 = buffer.data(gk + 270);
    const auto *gk_271 = buffer.data(gk + 271);
    const auto *gk_272 = buffer.data(gk + 272);
    const auto *gk_273 = buffer.data(gk + 273);
    const auto *gk_274 = buffer.data(gk + 274);
    const auto *gk_275 = buffer.data(gk + 275);
    const auto *gk_276 = buffer.data(gk + 276);
    const auto *gk_277 = buffer.data(gk + 277);
    const auto *gk_278 = buffer.data(gk + 278);
    const auto *gk_279 = buffer.data(gk + 279);
    const auto *gk_280 = buffer.data(gk + 280);
    const auto *gk_281 = buffer.data(gk + 281);
    const auto *gk_282 = buffer.data(gk + 282);
    const auto *gk_283 = buffer.data(gk + 283);
    const auto *gk_284 = buffer.data(gk + 284);
    const auto *gk_285 = buffer.data(gk + 285);
    const auto *gk_286 = buffer.data(gk + 286);
    const auto *gk_287 = buffer.data(gk + 287);
    const auto *gk_288 = buffer.data(gk + 288);
    const auto *gk_289 = buffer.data(gk + 289);
    const auto *gk_290 = buffer.data(gk + 290);
    const auto *gk_291 = buffer.data(gk + 291);
    const auto *gk_292 = buffer.data(gk + 292);
    const auto *gk_293 = buffer.data(gk + 293);
    const auto *gk_294 = buffer.data(gk + 294);
    const auto *gk_295 = buffer.data(gk + 295);
    const auto *gk_296 = buffer.data(gk + 296);
    const auto *gk_297 = buffer.data(gk + 297);
    const auto *gk_298 = buffer.data(gk + 298);
    const auto *gk_299 = buffer.data(gk + 299);
    const auto *gk_300 = buffer.data(gk + 300);
    const auto *gk_301 = buffer.data(gk + 301);
    const auto *gk_302 = buffer.data(gk + 302);
    const auto *gk_303 = buffer.data(gk + 303);
    const auto *gk_304 = buffer.data(gk + 304);
    const auto *gk_305 = buffer.data(gk + 305);
    const auto *gk_306 = buffer.data(gk + 306);
    const auto *gk_307 = buffer.data(gk + 307);
    const auto *gk_308 = buffer.data(gk + 308);
    const auto *gk_309 = buffer.data(gk + 309);
    const auto *gk_310 = buffer.data(gk + 310);
    const auto *gk_311 = buffer.data(gk + 311);
    const auto *gk_312 = buffer.data(gk + 312);
    const auto *gk_313 = buffer.data(gk + 313);
    const auto *gk_314 = buffer.data(gk + 314);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fk_0, gi0_0, gi1_0, \
                         gk_0, gk_1, gk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pb_y[k] * gk_0[k];

        t_2[k] = pb_z[k] * gk_0[k];

        t_3[k] = f_3 * gi0_0[k]
                 - f_4 * gi1_0[k]
                 + pb_y[k] * gk_1[k];

        t_4[k] = pb_y[k] * gk_2[k];

        t_5[k] = f_3 * gi0_0[k]
                 - f_4 * gi1_0[k]
                 + pb_z[k] * gk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, gi0_1, gi0_2, gi0_3, gi1_1, \
                         gi1_2, gi1_3, gk_3, gk_4, gk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * gi0_1[k]
                 - f_6 * gi1_1[k]
                 + pb_y[k] * gk_3[k];

        t_7[k] = pb_z[k] * gk_3[k];

        t_8[k] = pb_y[k] * gk_4[k];

        t_9[k] = f_5 * gi0_2[k]
                 - f_6 * gi1_2[k]
                 + pb_z[k] * gk_4[k];

        t_10[k] = f_7 * gi0_3[k]
                  - f_8 * gi1_3[k]
                  + pb_y[k] * gk_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, gi0_4, gi0_5, gi1_4, \
                         gi1_5, gk_5, gk_6, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gk_5[k];

        t_12[k] = f_3 * gi0_4[k]
                  - f_4 * gi1_4[k]
                  + pb_y[k] * gk_6[k];

        t_13[k] = pb_y[k] * gk_7[k];

        t_14[k] = f_7 * gi0_4[k]
                  - f_8 * gi1_4[k]
                  + pb_z[k] * gk_7[k];

        t_15[k] = f_9 * gi0_5[k]
                  - f_10 * gi1_5[k]
                  + pb_y[k] * gk_8[k];

        t_16[k] = pb_z[k] * gk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, gi0_6, gi0_7, gi1_6, gi1_7, gk_9, \
                         gk_10, gk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * gi0_6[k]
                  - f_6 * gi1_6[k]
                  + pb_y[k] * gk_9[k];

        t_18[k] = f_3 * gi0_7[k]
                  - f_4 * gi1_7[k]
                  + pb_y[k] * gk_10[k];

        t_19[k] = pb_y[k] * gk_11[k];

        t_20[k] = f_9 * gi0_7[k]
                  - f_10 * gi1_7[k]
                  + pb_z[k] * gk_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, gi0_8, gi0_9, gi0_10, gi1_8, \
                         gi1_9, gi1_10, gk_12, gk_13, gk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * gi0_8[k]
                  - f_12 * gi1_8[k]
                  + pb_y[k] * gk_12[k];

        t_22[k] = pb_z[k] * gk_12[k];

        t_23[k] = f_7 * gi0_9[k]
                  - f_8 * gi1_9[k]
                  + pb_y[k] * gk_13[k];

        t_24[k] = f_5 * gi0_10[k]
                  - f_6 * gi1_10[k]
                  + pb_y[k] * gk_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, fk_20, gi0_11, \
                         gi1_11, gk_15, gk_16, gk_17, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gi0_11[k]
                  - f_4 * gi1_11[k]
                  + pb_y[k] * gk_15[k];

        t_26[k] = pb_y[k] * gk_16[k];

        t_27[k] = f_11 * gi0_11[k]
                  - f_12 * gi1_11[k]
                  + pb_z[k] * gk_16[k];

        t_28[k] = f_0 * fk_20[k]
                  + pb_x[k] * gk_19[k];

        t_29[k] = pb_z[k] * gk_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, fk_22, fk_23, fk_24, fk_25, \
                         gk_18, gk_20, gk_21, gk_22, gk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * fk_22[k]
                  + pb_x[k] * gk_20[k];

        t_31[k] = f_0 * fk_23[k]
                  + pb_x[k] * gk_21[k];

        t_32[k] = f_0 * fk_24[k]
                  + pb_x[k] * gk_22[k];

        t_33[k] = f_0 * fk_25[k]
                  + pb_x[k] * gk_23[k];

        t_34[k] = pb_y[k] * gk_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, fk_27, gi0_12, gi0_13, \
                         gi1_12, gi1_13, gk_19, gk_20, gk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * fk_27[k]
                  + pb_x[k] * gk_25[k];

        t_36[k] = f_1 * gi0_12[k]
                  - f_2 * gi1_12[k]
                  + pb_y[k] * gk_19[k];

        t_37[k] = pb_z[k] * gk_19[k];

        t_38[k] = f_11 * gi0_13[k]
                  - f_12 * gi1_13[k]
                  + pb_y[k] * gk_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, gi0_14, gi0_15, gi0_16, gi1_14, gi1_15, \
                         gi1_16, gk_21, gk_22, gk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * gi0_14[k]
                  - f_10 * gi1_14[k]
                  + pb_y[k] * gk_21[k];

        t_40[k] = f_7 * gi0_15[k]
                  - f_8 * gi1_15[k]
                  + pb_y[k] * gk_22[k];

        t_41[k] = f_5 * gi0_16[k]
                  - f_6 * gi1_16[k]
                  + pb_y[k] * gk_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, fk_0, fl_0, \
                         gi0_17, gi1_17, gk_24, gk_25, gk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * gi0_17[k]
                  - f_4 * gi1_17[k]
                  + pb_y[k] * gk_24[k];

        t_43[k] = pb_y[k] * gk_25[k];

        t_44[k] = f_1 * gi0_17[k]
                  - f_2 * gi1_17[k]
                  + pb_z[k] * gk_25[k];

        t_45[k] = pa_y[k] * fl_0[k];

        t_46[k] = f_13 * fk_0[k]
                  + pb_y[k] * gk_26[k];

        t_47[k] = pb_z[k] * gk_26[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_z, fk_1, fk_3, fl_1, fl_2, \
                         fl_3, gk_27, gk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * fk_1[k]
                  + pa_y[k] * fl_1[k];

        t_49[k] = pb_z[k] * gk_27[k];

        t_50[k] = pa_y[k] * fl_2[k];

        t_51[k] = f_15 * fk_3[k]
                  + pa_y[k] * fl_3[k];

        t_52[k] = pb_z[k] * gk_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_y, pb_z, fk_4, fk_5, fk_7, \
                         fl_4, fl_5, fl_6, gk_29, gk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * fk_4[k]
                  + pb_y[k] * gk_29[k];

        t_54[k] = pa_y[k] * fl_4[k];

        t_55[k] = f_0 * fk_5[k]
                  + pa_y[k] * fl_5[k];

        t_56[k] = pb_z[k] * gk_30[k];

        t_57[k] = f_14 * fk_7[k]
                  + pa_y[k] * fl_6[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, fk_8, fk_9, fk_11, \
                         fl_7, fl_8, fl_9, gk_31, gk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_13 * fk_8[k]
                  + pb_y[k] * gk_31[k];

        t_59[k] = pa_y[k] * fl_7[k];

        t_60[k] = f_16 * fk_9[k]
                  + pa_y[k] * fl_8[k];

        t_61[k] = pb_z[k] * gk_32[k];

        t_62[k] = f_15 * fk_11[k]
                  + pa_y[k] * fl_9[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, fk_12, fk_13, fk_14, \
                         fl_10, fl_11, fl_12, gk_33, gk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * fk_12[k]
                  + pa_y[k] * fl_10[k];

        t_64[k] = f_13 * fk_13[k]
                  + pb_y[k] * gk_33[k];

        t_65[k] = pa_y[k] * fl_11[k];

        t_66[k] = f_17 * fk_14[k]
                  + pa_y[k] * fl_12[k];

        t_67[k] = pb_z[k] * gk_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pb_y, fk_16, fk_17, fk_18, fk_19, \
                         fl_13, fl_14, fl_15, fl_16, gk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * fk_16[k]
                  + pa_y[k] * fl_13[k];

        t_69[k] = f_15 * fk_17[k]
                  + pa_y[k] * fl_14[k];

        t_70[k] = f_14 * fk_18[k]
                  + pa_y[k] * fl_15[k];

        t_71[k] = f_13 * fk_19[k]
                  + pb_y[k] * gk_35[k];

        t_72[k] = pa_y[k] * fl_16[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, fk_37, fk_38, fk_39, fk_40, \
                         gk_36, gk_37, gk_38, gk_39, gk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * fk_37[k]
                  + pb_x[k] * gk_37[k];

        t_74[k] = pb_z[k] * gk_36[k];

        t_75[k] = f_15 * fk_38[k]
                  + pb_x[k] * gk_38[k];

        t_76[k] = f_15 * fk_39[k]
                  + pb_x[k] * gk_39[k];

        t_77[k] = f_15 * fk_40[k]
                  + pb_x[k] * gk_40[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, fk_20, fk_41, fk_42, \
                         fl_18, fl_19, gk_37, gk_41, gk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_15 * fk_41[k]
                  + pb_x[k] * gk_41[k];

        t_79[k] = f_15 * fk_42[k]
                  + pb_x[k] * gk_42[k];

        t_80[k] = pa_y[k] * fl_18[k];

        t_81[k] = f_18 * fk_20[k]
                  + pa_y[k] * fl_19[k];

        t_82[k] = pb_z[k] * gk_37[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, fk_22, fk_23, fk_24, fk_25, \
                         fk_26, fl_20, fl_21, fl_22, fl_23, fl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_17 * fk_22[k]
                  + pa_y[k] * fl_20[k];

        t_84[k] = f_16 * fk_23[k]
                  + pa_y[k] * fl_21[k];

        t_85[k] = f_0 * fk_24[k]
                  + pa_y[k] * fl_22[k];

        t_86[k] = f_15 * fk_25[k]
                  + pa_y[k] * fl_23[k];

        t_87[k] = f_14 * fk_26[k]
                  + pa_y[k] * fl_24[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, fk_0, fk_27, \
                         fl_0, fl_25, gk_43, gk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_13 * fk_27[k]
                  + pb_y[k] * gk_43[k];

        t_89[k] = pa_y[k] * fl_25[k];

        t_90[k] = pa_z[k] * fl_0[k];

        t_91[k] = pb_y[k] * gk_44[k];

        t_92[k] = f_13 * fk_0[k]
                  + pb_z[k] * gk_44[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, fk_2, fk_3, fl_1, \
                         fl_2, fl_3, gk_45, gk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * fl_1[k];

        t_94[k] = pb_y[k] * gk_45[k];

        t_95[k] = f_14 * fk_2[k]
                  + pa_z[k] * fl_2[k];

        t_96[k] = pa_z[k] * fl_3[k];

        t_97[k] = f_13 * fk_3[k]
                  + pb_z[k] * gk_46[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, fk_4, fk_5, fk_6, \
                         fl_4, fl_5, fl_6, gk_47, gk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * gk_47[k];

        t_99[k] = f_15 * fk_4[k]
                  + pa_z[k] * fl_4[k];

        t_100[k] = pa_z[k] * fl_5[k];

        t_101[k] = f_13 * fk_5[k]
                   + pb_z[k] * gk_48[k];

        t_102[k] = f_14 * fk_6[k]
                   + pa_z[k] * fl_6[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, fk_8, fk_9, \
                         fk_10, fl_7, fl_8, fl_9, gk_49, gk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * gk_49[k];

        t_104[k] = f_0 * fk_8[k]
                   + pa_z[k] * fl_7[k];

        t_105[k] = pa_z[k] * fl_8[k];

        t_106[k] = f_13 * fk_9[k]
                   + pb_z[k] * gk_50[k];

        t_107[k] = f_14 * fk_10[k]
                   + pa_z[k] * fl_9[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, fk_11, fk_13, \
                         fk_14, fl_10, fl_11, fl_12, gk_51, gk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * fk_11[k]
                   + pa_z[k] * fl_10[k];

        t_109[k] = pb_y[k] * gk_51[k];

        t_110[k] = f_16 * fk_13[k]
                   + pa_z[k] * fl_11[k];

        t_111[k] = pa_z[k] * fl_12[k];

        t_112[k] = f_13 * fk_14[k]
                   + pb_z[k] * gk_52[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, fk_15, fk_16, fk_17, \
                         fk_19, fl_13, fl_14, fl_15, fl_16, gk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * fk_15[k]
                   + pa_z[k] * fl_13[k];

        t_114[k] = f_15 * fk_16[k]
                   + pa_z[k] * fl_14[k];

        t_115[k] = f_0 * fk_17[k]
                   + pa_z[k] * fl_15[k];

        t_116[k] = pb_y[k] * gk_53[k];

        t_117[k] = f_17 * fk_19[k]
                   + pa_z[k] * fl_16[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, fk_61, fk_62, fk_63, \
                         fk_64, fl_17, gk_56, gk_57, gk_58, gk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * fl_17[k];

        t_119[k] = f_15 * fk_61[k]
                   + pb_x[k] * gk_56[k];

        t_120[k] = f_15 * fk_62[k]
                   + pb_x[k] * gk_57[k];

        t_121[k] = f_15 * fk_63[k]
                   + pb_x[k] * gk_58[k];

        t_122[k] = f_15 * fk_64[k]
                   + pb_x[k] * gk_59[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, pb_y, fk_65, fk_67, fl_19, \
                         gk_54, gk_60, gk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_15 * fk_65[k]
                   + pb_x[k] * gk_60[k];

        t_124[k] = pb_y[k] * gk_54[k];

        t_125[k] = f_15 * fk_67[k]
                   + pb_x[k] * gk_61[k];

        t_126[k] = pa_z[k] * fl_19[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_z, pb_z, fk_20, fk_21, fk_22, fk_23, \
                         fl_20, fl_21, fl_22, gk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * fk_20[k]
                   + pb_z[k] * gk_55[k];

        t_128[k] = f_14 * fk_21[k]
                   + pa_z[k] * fl_20[k];

        t_129[k] = f_15 * fk_22[k]
                   + pa_z[k] * fl_21[k];

        t_130[k] = f_0 * fk_23[k]
                   + pa_z[k] * fl_22[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_z, pb_y, fk_24, fk_25, fk_27, fl_23, \
                         fl_24, fl_25, gk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_16 * fk_24[k]
                   + pa_z[k] * fl_23[k];

        t_132[k] = f_17 * fk_25[k]
                   + pa_z[k] * fl_24[k];

        t_133[k] = pb_y[k] * gk_61[k];

        t_134[k] = f_18 * fk_27[k]
                   + pa_z[k] * fl_25[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pb_y, pb_z, dl0_0, dl1_0, fk_28, fl_26, \
                         gk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_19 * dl0_0[k]
                   - f_20 * dl1_0[k]
                   + pa_y[k] * fl_26[k];

        t_136[k] = f_14 * fk_28[k]
                   + pb_y[k] * gk_62[k];

        t_137[k] = pb_z[k] * gk_62[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_x, pb_z, fk_69, gi0_18, gi0_20, gi1_18, \
                         gi1_20, gk_63, gk_64, gk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_14 * fk_69[k]
                   + f_11 * gi0_20[k]
                   - f_12 * gi1_20[k]
                   + pb_x[k] * gk_65[k];

        t_139[k] = pb_z[k] * gk_63[k];

        t_140[k] = f_3 * gi0_18[k]
                   - f_4 * gi1_18[k]
                   + pb_z[k] * gk_64[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_x, pb_y, pb_z, fk_30, fk_71, gi0_19, \
                         gi0_22, gi1_19, gi1_22, gk_65, gk_66, gk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_14 * fk_71[k]
                   + f_9 * gi0_22[k]
                   - f_10 * gi1_22[k]
                   + pb_x[k] * gk_67[k];

        t_142[k] = pb_z[k] * gk_65[k];

        t_143[k] = f_14 * fk_30[k]
                   + pb_y[k] * gk_66[k];

        t_144[k] = f_5 * gi0_19[k]
                   - f_6 * gi1_19[k]
                   + pb_z[k] * gk_66[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_z, fk_73, gi0_20, gi0_25, gi1_20, \
                         gi1_25, gk_67, gk_68, gk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_14 * fk_73[k]
                   + f_7 * gi0_25[k]
                   - f_8 * gi1_25[k]
                   + pb_x[k] * gk_70[k];

        t_146[k] = pb_z[k] * gk_67[k];

        t_147[k] = f_3 * gi0_20[k]
                   - f_4 * gi1_20[k]
                   + pb_z[k] * gk_68[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pb_y, pb_z, fk_32, fk_75, gi0_21, \
                         gi0_29, gi1_21, gi1_29, gk_69, gk_70, gk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * fk_32[k]
                   + pb_y[k] * gk_69[k];

        t_149[k] = f_7 * gi0_21[k]
                   - f_8 * gi1_21[k]
                   + pb_z[k] * gk_69[k];

        t_150[k] = f_14 * fk_75[k]
                   + f_5 * gi0_29[k]
                   - f_6 * gi1_29[k]
                   + pb_x[k] * gk_74[k];

        t_151[k] = pb_z[k] * gk_70[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, fk_34, gi0_22, gi0_23, \
                         gi0_24, gi1_22, gi1_23, gi1_24, gk_71, gk_72, \
                         gk_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * gi0_22[k]
                   - f_4 * gi1_22[k]
                   + pb_z[k] * gk_71[k];

        t_153[k] = f_5 * gi0_23[k]
                   - f_6 * gi1_23[k]
                   + pb_z[k] * gk_72[k];

        t_154[k] = f_14 * fk_34[k]
                   + pb_y[k] * gk_73[k];

        t_155[k] = f_9 * gi0_24[k]
                   - f_10 * gi1_24[k]
                   + pb_z[k] * gk_73[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_z, fk_77, gi0_25, gi0_30, gi1_25, \
                         gi1_30, gk_74, gk_75, gk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_14 * fk_77[k]
                   + f_3 * gi0_30[k]
                   - f_4 * gi1_30[k]
                   + pb_x[k] * gk_79[k];

        t_157[k] = pb_z[k] * gk_74[k];

        t_158[k] = f_3 * gi0_25[k]
                   - f_4 * gi1_25[k]
                   + pb_z[k] * gk_75[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_y, pb_z, fk_36, gi0_26, gi0_27, \
                         gi0_28, gi1_26, gi1_27, gi1_28, gk_76, gk_77, \
                         gk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * gi0_26[k]
                   - f_6 * gi1_26[k]
                   + pb_z[k] * gk_76[k];

        t_160[k] = f_7 * gi0_27[k]
                   - f_8 * gi1_27[k]
                   + pb_z[k] * gk_77[k];

        t_161[k] = f_14 * fk_36[k]
                   + pb_y[k] * gk_78[k];

        t_162[k] = f_11 * gi0_28[k]
                   - f_12 * gi1_28[k]
                   + pb_z[k] * gk_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_x, pb_z, fk_78, fk_79, fk_80, \
                         fk_81, gk_79, gk_80, gk_82, gk_83, gk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_14 * fk_78[k]
                   + pb_x[k] * gk_80[k];

        t_164[k] = pb_z[k] * gk_79[k];

        t_165[k] = f_14 * fk_79[k]
                   + pb_x[k] * gk_82[k];

        t_166[k] = f_14 * fk_80[k]
                   + pb_x[k] * gk_83[k];

        t_167[k] = f_14 * fk_81[k]
                   + pb_x[k] * gk_84[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pb_x, dl0_1, dl1_1, fk_82, fk_83, \
                         fk_84, fl_63, gk_85, gk_86, gk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_14 * fk_82[k]
                   + pb_x[k] * gk_85[k];

        t_169[k] = f_14 * fk_83[k]
                   + pb_x[k] * gk_86[k];

        t_170[k] = f_14 * fk_84[k]
                   + pb_x[k] * gk_87[k];

        t_171[k] = f_19 * dl0_1[k]
                   - f_20 * dl1_1[k]
                   + pa_x[k] * fl_63[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, gi0_30, gi0_31, gi0_32, gi1_30, \
                         gi1_31, gi1_32, gk_80, gk_81, gk_82, gk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * gk_80[k];

        t_173[k] = f_3 * gi0_30[k]
                   - f_4 * gi1_30[k]
                   + pb_z[k] * gk_81[k];

        t_174[k] = f_5 * gi0_31[k]
                   - f_6 * gi1_31[k]
                   + pb_z[k] * gk_82[k];

        t_175[k] = f_7 * gi0_32[k]
                   - f_8 * gi1_32[k]
                   + pb_z[k] * gk_83[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pb_z, fk_43, gi0_33, gi0_34, \
                         gi0_35, gi1_33, gi1_34, gi1_35, gk_84, gk_85, \
                         gk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * gi0_33[k]
                   - f_10 * gi1_33[k]
                   + pb_z[k] * gk_84[k];

        t_177[k] = f_11 * gi0_34[k]
                   - f_12 * gi1_34[k]
                   + pb_z[k] * gk_85[k];

        t_178[k] = f_14 * fk_43[k]
                   + pb_y[k] * gk_87[k];

        t_179[k] = f_1 * gi0_35[k]
                   - f_2 * gi1_35[k]
                   + pb_z[k] * gk_87[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pa_y, pa_z, pb_y, fk_45, \
                         fl_27, fl_28, fl_35, fl_36, fl_37, gk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * fl_35[k];

        t_181[k] = pa_z[k] * fl_27[k];

        t_182[k] = pa_y[k] * fl_36[k];

        t_183[k] = pa_z[k] * fl_28[k];

        t_184[k] = f_13 * fk_45[k]
                   + pb_y[k] * gk_88[k];

        t_185[k] = pa_y[k] * fl_37[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, fk_29, \
                         fk_47, fl_29, fl_30, fl_38, gk_89, gk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * fl_29[k];

        t_187[k] = f_13 * fk_29[k]
                   + pb_z[k] * gk_89[k];

        t_188[k] = f_13 * fk_47[k]
                   + pb_y[k] * gk_90[k];

        t_189[k] = pa_y[k] * fl_38[k];

        t_190[k] = pa_z[k] * fl_30[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, fk_31, fk_49, fk_50, \
                         fl_39, fl_40, gk_91, gk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * fk_31[k]
                   + pb_z[k] * gk_91[k];

        t_192[k] = f_14 * fk_49[k]
                   + pa_y[k] * fl_39[k];

        t_193[k] = f_13 * fk_50[k]
                   + pb_y[k] * gk_92[k];

        t_194[k] = pa_y[k] * fl_40[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pa_z, pb_z, fk_33, fk_52, fk_53, \
                         fl_31, fl_41, fl_42, gk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * fl_31[k];

        t_196[k] = f_13 * fk_33[k]
                   + pb_z[k] * gk_93[k];

        t_197[k] = f_15 * fk_52[k]
                   + pa_y[k] * fl_41[k];

        t_198[k] = f_14 * fk_53[k]
                   + pa_y[k] * fl_42[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, fk_35, fk_54, \
                         fl_32, fl_43, gk_94, gk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * fk_54[k]
                   + pb_y[k] * gk_94[k];

        t_200[k] = pa_y[k] * fl_43[k];

        t_201[k] = pa_z[k] * fl_32[k];

        t_202[k] = f_13 * fk_35[k]
                   + pb_z[k] * gk_95[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pa_y, pb_y, fk_56, fk_57, fk_58, \
                         fk_59, fl_44, fl_45, fl_46, fl_47, gk_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * fk_56[k]
                   + pa_y[k] * fl_44[k];

        t_204[k] = f_15 * fk_57[k]
                   + pa_y[k] * fl_45[k];

        t_205[k] = f_14 * fk_58[k]
                   + pa_y[k] * fl_46[k];

        t_206[k] = f_13 * fk_59[k]
                   + pb_y[k] * gk_96[k];

        t_207[k] = pa_y[k] * fl_47[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pa_z, pb_x, fk_94, fk_95, fk_96, \
                         fk_97, fl_33, gk_98, gk_99, gk_100, gk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_z[k] * fl_33[k];

        t_209[k] = f_14 * fk_94[k]
                   + pb_x[k] * gk_98[k];

        t_210[k] = f_14 * fk_95[k]
                   + pb_x[k] * gk_99[k];

        t_211[k] = f_14 * fk_96[k]
                   + pb_x[k] * gk_100[k];

        t_212[k] = f_14 * fk_97[k]
                   + pb_x[k] * gk_101[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pa_z, pb_x, fk_98, fk_99, fl_34, \
                         fl_48, gk_102, gk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_14 * fk_98[k]
                   + pb_x[k] * gk_102[k];

        t_214[k] = f_14 * fk_99[k]
                   + pb_x[k] * gk_103[k];

        t_215[k] = pa_y[k] * fl_48[k];

        t_216[k] = pa_z[k] * fl_34[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_y, pb_z, fk_37, fk_62, fk_63, fk_64, \
                         fl_49, fl_50, fl_51, gk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_13 * fk_37[k]
                   + pb_z[k] * gk_97[k];

        t_218[k] = f_17 * fk_62[k]
                   + pa_y[k] * fl_49[k];

        t_219[k] = f_16 * fk_63[k]
                   + pa_y[k] * fl_50[k];

        t_220[k] = f_0 * fk_64[k]
                   + pa_y[k] * fl_51[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pb_y, fk_65, fk_66, fk_67, fl_52, \
                         fl_53, fl_54, gk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_15 * fk_65[k]
                   + pa_y[k] * fl_52[k];

        t_222[k] = f_14 * fk_66[k]
                   + pa_y[k] * fl_53[k];

        t_223[k] = f_13 * fk_67[k]
                   + pb_y[k] * gk_104[k];

        t_224[k] = pa_y[k] * fl_54[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_z, pb_y, pb_z, dl0_0, dl1_0, fk_44, \
                         fl_35, gi0_36, gi1_36, gk_105, gk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_19 * dl0_0[k]
                   - f_20 * dl1_0[k]
                   + pa_z[k] * fl_35[k];

        t_226[k] = pb_y[k] * gk_105[k];

        t_227[k] = f_14 * fk_44[k]
                   + pb_z[k] * gk_105[k];

        t_228[k] = f_3 * gi0_36[k]
                   - f_4 * gi1_36[k]
                   + pb_y[k] * gk_106[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, fk_46, fk_103, gi0_37, \
                         gi0_39, gi1_37, gi1_39, gk_107, gk_108, \
                         gk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * gk_107[k];

        t_230[k] = f_14 * fk_103[k]
                   + f_11 * gi0_39[k]
                   - f_12 * gi1_39[k]
                   + pb_x[k] * gk_109[k];

        t_231[k] = f_5 * gi0_37[k]
                   - f_6 * gi1_37[k]
                   + pb_y[k] * gk_108[k];

        t_232[k] = f_14 * fk_46[k]
                   + pb_z[k] * gk_108[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pb_y, pb_z, fk_48, fk_105, gi0_38, \
                         gi0_42, gi1_38, gi1_42, gk_109, gk_110, \
                         gk_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * gk_109[k];

        t_234[k] = f_14 * fk_105[k]
                   + f_9 * gi0_42[k]
                   - f_10 * gi1_42[k]
                   + pb_x[k] * gk_112[k];

        t_235[k] = f_7 * gi0_38[k]
                   - f_8 * gi1_38[k]
                   + pb_y[k] * gk_110[k];

        t_236[k] = f_14 * fk_48[k]
                   + pb_z[k] * gk_110[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_x, pb_y, fk_107, gi0_39, gi0_46, gi1_39, \
                         gi1_46, gk_111, gk_112, gk_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_3 * gi0_39[k]
                   - f_4 * gi1_39[k]
                   + pb_y[k] * gk_111[k];

        t_238[k] = pb_y[k] * gk_112[k];

        t_239[k] = f_14 * fk_107[k]
                   + f_7 * gi0_46[k]
                   - f_8 * gi1_46[k]
                   + pb_x[k] * gk_116[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_y, pb_z, fk_51, gi0_40, gi0_41, \
                         gi0_42, gi1_40, gi1_41, gi1_42, gk_113, gk_114, \
                         gk_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_9 * gi0_40[k]
                   - f_10 * gi1_40[k]
                   + pb_y[k] * gk_113[k];

        t_241[k] = f_14 * fk_51[k]
                   + pb_z[k] * gk_113[k];

        t_242[k] = f_5 * gi0_41[k]
                   - f_6 * gi1_41[k]
                   + pb_y[k] * gk_114[k];

        t_243[k] = f_3 * gi0_42[k]
                   - f_4 * gi1_42[k]
                   + pb_y[k] * gk_115[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pb_y, pb_z, fk_55, fk_109, gi0_43, \
                         gi0_47, gi1_43, gi1_47, gk_116, gk_117, \
                         gk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_y[k] * gk_116[k];

        t_245[k] = f_14 * fk_109[k]
                   + f_5 * gi0_47[k]
                   - f_6 * gi1_47[k]
                   + pb_x[k] * gk_121[k];

        t_246[k] = f_11 * gi0_43[k]
                   - f_12 * gi1_43[k]
                   + pb_y[k] * gk_117[k];

        t_247[k] = f_14 * fk_55[k]
                   + pb_z[k] * gk_117[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_y, gi0_44, gi0_45, gi0_46, gi1_44, \
                         gi1_45, gi1_46, gk_118, gk_119, gk_120, \
                         gk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * gi0_44[k]
                   - f_8 * gi1_44[k]
                   + pb_y[k] * gk_118[k];

        t_249[k] = f_5 * gi0_45[k]
                   - f_6 * gi1_45[k]
                   + pb_y[k] * gk_119[k];

        t_250[k] = f_3 * gi0_46[k]
                   - f_4 * gi1_46[k]
                   + pb_y[k] * gk_120[k];

        t_251[k] = pb_y[k] * gk_121[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, fk_110, fk_111, fk_112, fk_113, \
                         gi0_53, gi1_53, gk_122, gk_123, gk_124, \
                         gk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * fk_110[k]
                   + f_3 * gi0_53[k]
                   - f_4 * gi1_53[k]
                   + pb_x[k] * gk_122[k];

        t_253[k] = f_14 * fk_111[k]
                   + pb_x[k] * gk_123[k];

        t_254[k] = f_14 * fk_112[k]
                   + pb_x[k] * gk_124[k];

        t_255[k] = f_14 * fk_113[k]
                   + pb_x[k] * gk_125[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, fk_114, fk_115, \
                         fk_116, fk_117, gk_122, gk_126, gk_127, gk_128, \
                         gk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * fk_114[k]
                   + pb_x[k] * gk_126[k];

        t_257[k] = f_14 * fk_115[k]
                   + pb_x[k] * gk_127[k];

        t_258[k] = f_14 * fk_116[k]
                   + pb_x[k] * gk_128[k];

        t_259[k] = pb_y[k] * gk_122[k];

        t_260[k] = f_14 * fk_117[k]
                   + pb_x[k] * gk_130[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_y, pb_z, fk_60, gi0_48, gi0_49, \
                         gi0_50, gi1_48, gi1_49, gi1_50, gk_123, gk_125, \
                         gk_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * gi0_48[k]
                   - f_2 * gi1_48[k]
                   + pb_y[k] * gk_123[k];

        t_262[k] = f_14 * fk_60[k]
                   + pb_z[k] * gk_123[k];

        t_263[k] = f_11 * gi0_49[k]
                   - f_12 * gi1_49[k]
                   + pb_y[k] * gk_125[k];

        t_264[k] = f_9 * gi0_50[k]
                   - f_10 * gi1_50[k]
                   + pb_y[k] * gk_126[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_y, gi0_51, gi0_52, gi0_53, gi1_51, \
                         gi1_52, gi1_53, gk_127, gk_128, gk_129, \
                         gk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_7 * gi0_51[k]
                   - f_8 * gi1_51[k]
                   + pb_y[k] * gk_127[k];

        t_266[k] = f_5 * gi0_52[k]
                   - f_6 * gi1_52[k]
                   + pb_y[k] * gk_128[k];

        t_267[k] = f_3 * gi0_53[k]
                   - f_4 * gi1_53[k]
                   + pb_y[k] * gk_129[k];

        t_268[k] = pb_y[k] * gk_130[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pb_y, pb_z, dl0_2, dl1_2, fk_68, \
                         fk_118, fl_72, fl_73, gk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_19 * dl0_2[k]
                   - f_20 * dl1_2[k]
                   + pa_x[k] * fl_72[k];

        t_270[k] = f_18 * fk_118[k]
                   + pa_x[k] * fl_73[k];

        t_271[k] = f_15 * fk_68[k]
                   + pb_y[k] * gk_131[k];

        t_272[k] = pb_z[k] * gk_131[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, pa_x, pb_z, fk_120, fk_121, \
                         fk_122, fl_75, fl_76, fl_77, gk_132, gk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_17 * fk_120[k]
                   + pa_x[k] * fl_75[k];

        t_274[k] = pb_z[k] * gk_132[k];

        t_275[k] = f_17 * fk_121[k]
                   + pa_x[k] * fl_76[k];

        t_276[k] = f_16 * fk_122[k]
                   + pa_x[k] * fl_77[k];

        t_277[k] = pb_z[k] * gk_133[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_x, pb_y, pb_z, fk_70, fk_124, fk_125, \
                         fl_78, fl_79, gk_134, gk_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_15 * fk_70[k]
                   + pb_y[k] * gk_134[k];

        t_279[k] = f_16 * fk_124[k]
                   + pa_x[k] * fl_78[k];

        t_280[k] = f_0 * fk_125[k]
                   + pa_x[k] * fl_79[k];

        t_281[k] = pb_z[k] * gk_135[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_x, pb_y, fk_72, fk_127, fk_128, \
                         fk_129, fl_80, fl_81, fl_82, gk_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * fk_127[k]
                   + pa_x[k] * fl_80[k];

        t_283[k] = f_15 * fk_72[k]
                   + pb_y[k] * gk_136[k];

        t_284[k] = f_0 * fk_128[k]
                   + pa_x[k] * fl_81[k];

        t_285[k] = f_15 * fk_129[k]
                   + pa_x[k] * fl_82[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_x, pb_y, pb_z, fk_74, fk_131, fk_132, \
                         fl_83, fl_84, gk_137, gk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = pb_z[k] * gk_137[k];

        t_287[k] = f_15 * fk_131[k]
                   + pa_x[k] * fl_83[k];

        t_288[k] = f_15 * fk_132[k]
                   + pa_x[k] * fl_84[k];

        t_289[k] = f_15 * fk_74[k]
                   + pb_y[k] * gk_138[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pa_x, pb_z, fk_133, fk_134, \
                         fk_135, fk_136, fl_85, fl_86, fl_87, fl_88, \
                         gk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_15 * fk_133[k]
                   + pa_x[k] * fl_85[k];

        t_291[k] = f_14 * fk_134[k]
                   + pa_x[k] * fl_86[k];

        t_292[k] = pb_z[k] * gk_139[k];

        t_293[k] = f_14 * fk_135[k]
                   + pa_x[k] * fl_87[k];

        t_294[k] = f_14 * fk_136[k]
                   + pa_x[k] * fl_88[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pa_x, pb_x, pb_y, fk_76, fk_137, fk_138, \
                         fk_139, fl_89, fl_90, gk_140, gk_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_14 * fk_137[k]
                   + pa_x[k] * fl_89[k];

        t_296[k] = f_15 * fk_76[k]
                   + pb_y[k] * gk_140[k];

        t_297[k] = f_14 * fk_138[k]
                   + pa_x[k] * fl_90[k];

        t_298[k] = f_13 * fk_139[k]
                   + pb_x[k] * gk_142[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pb_x, pb_z, fk_141, fk_142, \
                         fk_143, fk_144, gk_141, gk_143, gk_144, gk_145, \
                         gk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pb_z[k] * gk_141[k];

        t_300[k] = f_13 * fk_141[k]
                   + pb_x[k] * gk_143[k];

        t_301[k] = f_13 * fk_142[k]
                   + pb_x[k] * gk_144[k];

        t_302[k] = f_13 * fk_143[k]
                   + pb_x[k] * gk_145[k];

        t_303[k] = f_13 * fk_144[k]
                   + pb_x[k] * gk_146[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pa_x, pb_x, pb_z, fk_145, fk_146, \
                         fl_91, fl_92, gk_142, gk_147, gk_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_13 * fk_145[k]
                   + pb_x[k] * gk_147[k];

        t_305[k] = f_13 * fk_146[k]
                   + pb_x[k] * gk_148[k];

        t_306[k] = pa_x[k] * fl_91[k];

        t_307[k] = pb_z[k] * gk_142[k];

        t_308[k] = pa_x[k] * fl_92[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, t_314, t_315, pa_x, pa_z, fl_55, \
                         fl_93, fl_94, fl_95, fl_96, fl_97, fl_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_x[k] * fl_93[k];

        t_310[k] = pa_x[k] * fl_94[k];

        t_311[k] = pa_x[k] * fl_95[k];

        t_312[k] = pa_x[k] * fl_96[k];

        t_313[k] = pa_x[k] * fl_97[k];

        t_314[k] = pa_x[k] * fl_98[k];

        t_315[k] = pa_z[k] * fl_55[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_z, pb_y, pb_z, fk_68, fk_85, fl_56, \
                         fl_57, gk_149, gk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_z[k] * fl_56[k];

        t_317[k] = f_13 * fk_68[k]
                   + pb_z[k] * gk_149[k];

        t_318[k] = pa_z[k] * fl_57[k];

        t_319[k] = f_14 * fk_85[k]
                   + pb_y[k] * gk_150[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_x, pa_z, pb_y, pb_z, fk_69, fk_87, \
                         fk_150, fl_58, fl_99, gk_151, gk_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_17 * fk_150[k]
                   + pa_x[k] * fl_99[k];

        t_321[k] = pa_z[k] * fl_58[k];

        t_322[k] = f_13 * fk_69[k]
                   + pb_z[k] * gk_151[k];

        t_323[k] = f_14 * fk_87[k]
                   + pb_y[k] * gk_152[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_x, pa_z, pb_z, fk_71, fk_152, fk_154, \
                         fl_59, fl_100, fl_101, gk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_16 * fk_152[k]
                   + pa_x[k] * fl_100[k];

        t_325[k] = pa_z[k] * fl_59[k];

        t_326[k] = f_13 * fk_71[k]
                   + pb_z[k] * gk_153[k];

        t_327[k] = f_0 * fk_154[k]
                   + pa_x[k] * fl_101[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_x, pa_z, pb_y, pb_z, fk_73, fk_89, \
                         fk_155, fl_60, fl_102, gk_154, gk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_14 * fk_89[k]
                   + pb_y[k] * gk_154[k];

        t_329[k] = f_0 * fk_155[k]
                   + pa_x[k] * fl_102[k];

        t_330[k] = pa_z[k] * fl_60[k];

        t_331[k] = f_13 * fk_73[k]
                   + pb_z[k] * gk_155[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pa_x, pb_y, fk_91, fk_157, fk_158, \
                         fk_159, fl_103, fl_104, fl_105, gk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_15 * fk_157[k]
                   + pa_x[k] * fl_103[k];

        t_333[k] = f_15 * fk_158[k]
                   + pa_x[k] * fl_104[k];

        t_334[k] = f_14 * fk_91[k]
                   + pb_y[k] * gk_156[k];

        t_335[k] = f_15 * fk_159[k]
                   + pa_x[k] * fl_105[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_x, pa_z, pb_z, fk_75, fk_160, fk_161, \
                         fl_61, fl_106, fl_107, gk_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_z[k] * fl_61[k];

        t_337[k] = f_13 * fk_75[k]
                   + pb_z[k] * gk_157[k];

        t_338[k] = f_14 * fk_160[k]
                   + pa_x[k] * fl_106[k];

        t_339[k] = f_14 * fk_161[k]
                   + pa_x[k] * fl_107[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_x, pa_z, pb_y, fk_93, fk_162, fk_163, \
                         fl_62, fl_108, fl_109, gk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_14 * fk_162[k]
                   + pa_x[k] * fl_108[k];

        t_341[k] = f_14 * fk_93[k]
                   + pb_y[k] * gk_158[k];

        t_342[k] = f_14 * fk_163[k]
                   + pa_x[k] * fl_109[k];

        t_343[k] = pa_z[k] * fl_62[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, pb_x, fk_165, fk_166, fk_167, \
                         fk_168, fk_169, gk_159, gk_160, gk_161, gk_162, \
                         gk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_13 * fk_165[k]
                   + pb_x[k] * gk_159[k];

        t_345[k] = f_13 * fk_166[k]
                   + pb_x[k] * gk_160[k];

        t_346[k] = f_13 * fk_167[k]
                   + pb_x[k] * gk_161[k];

        t_347[k] = f_13 * fk_168[k]
                   + pb_x[k] * gk_162[k];

        t_348[k] = f_13 * fk_169[k]
                   + pb_x[k] * gk_163[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, t_354, pa_x, pb_x, fk_170, fk_171, \
                         fl_110, fl_111, fl_112, fl_113, gk_164, \
                         gk_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_13 * fk_170[k]
                   + pb_x[k] * gk_164[k];

        t_350[k] = f_13 * fk_171[k]
                   + pb_x[k] * gk_165[k];

        t_351[k] = pa_x[k] * fl_110[k];

        t_352[k] = pa_x[k] * fl_111[k];

        t_353[k] = pa_x[k] * fl_112[k];

        t_354[k] = pa_x[k] * fl_113[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, t_360, pa_x, pa_y, fl_64, fl_114, \
                         fl_115, fl_116, fl_117, fl_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = pa_x[k] * fl_114[k];

        t_356[k] = pa_x[k] * fl_115[k];

        t_357[k] = pa_x[k] * fl_116[k];

        t_358[k] = pa_x[k] * fl_117[k];

        t_359[k] = pa_x[k] * fl_118[k];

        t_360[k] = pa_y[k] * fl_64[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, t_365, pa_x, pa_y, pb_y, fk_100, fk_101, \
                         fk_174, fl_65, fl_66, fl_119, gk_166, gk_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_13 * fk_100[k]
                   + pb_y[k] * gk_166[k];

        t_362[k] = pa_y[k] * fl_65[k];

        t_363[k] = f_17 * fk_174[k]
                   + pa_x[k] * fl_119[k];

        t_364[k] = f_13 * fk_101[k]
                   + pb_y[k] * gk_167[k];

        t_365[k] = pa_y[k] * fl_66[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_x, pa_y, pb_y, pb_z, fk_86, fk_103, \
                         fk_176, fl_67, fl_120, gk_168, gk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_16 * fk_176[k]
                   + pa_x[k] * fl_120[k];

        t_367[k] = f_14 * fk_86[k]
                   + pb_z[k] * gk_168[k];

        t_368[k] = f_13 * fk_103[k]
                   + pb_y[k] * gk_169[k];

        t_369[k] = pa_y[k] * fl_67[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pb_y, pb_z, fk_88, fk_105, fk_178, \
                         fk_179, fl_121, fl_122, gk_170, gk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_0 * fk_178[k]
                   + pa_x[k] * fl_121[k];

        t_371[k] = f_14 * fk_88[k]
                   + pb_z[k] * gk_170[k];

        t_372[k] = f_0 * fk_179[k]
                   + pa_x[k] * fl_122[k];

        t_373[k] = f_13 * fk_105[k]
                   + pb_y[k] * gk_171[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pa_x, pa_y, pb_z, fk_90, fk_181, fk_182, \
                         fl_68, fl_123, fl_124, gk_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * fl_68[k];

        t_375[k] = f_15 * fk_181[k]
                   + pa_x[k] * fl_123[k];

        t_376[k] = f_14 * fk_90[k]
                   + pb_z[k] * gk_172[k];

        t_377[k] = f_15 * fk_182[k]
                   + pa_x[k] * fl_124[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_x, pa_y, pb_y, fk_107, fk_183, fk_185, \
                         fl_69, fl_125, fl_126, gk_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_15 * fk_183[k]
                   + pa_x[k] * fl_125[k];

        t_379[k] = f_13 * fk_107[k]
                   + pb_y[k] * gk_173[k];

        t_380[k] = pa_y[k] * fl_69[k];

        t_381[k] = f_14 * fk_185[k]
                   + pa_x[k] * fl_126[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_x, pb_z, fk_92, fk_186, fk_187, \
                         fk_188, fl_127, fl_128, fl_129, gk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_14 * fk_92[k]
                   + pb_z[k] * gk_174[k];

        t_383[k] = f_14 * fk_186[k]
                   + pa_x[k] * fl_127[k];

        t_384[k] = f_14 * fk_187[k]
                   + pa_x[k] * fl_128[k];

        t_385[k] = f_14 * fk_188[k]
                   + pa_x[k] * fl_129[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_y, pb_x, pb_y, fk_109, fk_189, fk_190, \
                         fl_70, gk_175, gk_176, gk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_13 * fk_109[k]
                   + pb_y[k] * gk_175[k];

        t_387[k] = pa_y[k] * fl_70[k];

        t_388[k] = f_13 * fk_189[k]
                   + pb_x[k] * gk_176[k];

        t_389[k] = f_13 * fk_190[k]
                   + pb_x[k] * gk_177[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, pb_x, fk_191, fk_192, fk_193, \
                         fk_194, fk_195, gk_178, gk_179, gk_180, gk_181, \
                         gk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_13 * fk_191[k]
                   + pb_x[k] * gk_178[k];

        t_391[k] = f_13 * fk_192[k]
                   + pb_x[k] * gk_179[k];

        t_392[k] = f_13 * fk_193[k]
                   + pb_x[k] * gk_180[k];

        t_393[k] = f_13 * fk_194[k]
                   + pb_x[k] * gk_181[k];

        t_394[k] = f_13 * fk_195[k]
                   + pb_x[k] * gk_182[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, t_400, t_401, pa_x, pa_y, fl_71, \
                         fl_130, fl_131, fl_132, fl_133, fl_134, \
                         fl_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * fl_71[k];

        t_396[k] = pa_x[k] * fl_130[k];

        t_397[k] = pa_x[k] * fl_131[k];

        t_398[k] = pa_x[k] * fl_132[k];

        t_399[k] = pa_x[k] * fl_133[k];

        t_400[k] = pa_x[k] * fl_134[k];

        t_401[k] = pa_x[k] * fl_135[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, t_407, pa_x, pb_y, pb_z, fk_100, \
                         fk_197, fl_136, fl_137, fl_138, fl_139, \
                         gk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_x[k] * fl_136[k];

        t_403[k] = pa_x[k] * fl_137[k];

        t_404[k] = pa_x[k] * fl_138[k];

        t_405[k] = f_18 * fk_197[k]
                   + pa_x[k] * fl_139[k];

        t_406[k] = pb_y[k] * gk_183[k];

        t_407[k] = f_15 * fk_100[k]
                   + pb_z[k] * gk_183[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pa_x, pb_y, fk_200, fk_201, fk_202, \
                         fl_141, fl_142, fl_143, gk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_17 * fk_200[k]
                   + pa_x[k] * fl_141[k];

        t_409[k] = pb_y[k] * gk_184[k];

        t_410[k] = f_17 * fk_201[k]
                   + pa_x[k] * fl_142[k];

        t_411[k] = f_16 * fk_202[k]
                   + pa_x[k] * fl_143[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_x, pb_y, pb_z, fk_102, fk_204, fk_205, \
                         fl_144, fl_145, gk_185, gk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_15 * fk_102[k]
                   + pb_z[k] * gk_185[k];

        t_413[k] = pb_y[k] * gk_186[k];

        t_414[k] = f_16 * fk_204[k]
                   + pa_x[k] * fl_144[k];

        t_415[k] = f_0 * fk_205[k]
                   + pa_x[k] * fl_145[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pa_x, pb_y, pb_z, fk_104, fk_206, fk_208, \
                         fl_146, fl_147, gk_187, gk_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_15 * fk_104[k]
                   + pb_z[k] * gk_187[k];

        t_417[k] = f_0 * fk_206[k]
                   + pa_x[k] * fl_146[k];

        t_418[k] = pb_y[k] * gk_188[k];

        t_419[k] = f_0 * fk_208[k]
                   + pa_x[k] * fl_147[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pa_x, pb_z, fk_106, fk_209, fk_210, \
                         fk_211, fl_148, fl_149, fl_150, gk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_15 * fk_209[k]
                   + pa_x[k] * fl_148[k];

        t_421[k] = f_15 * fk_106[k]
                   + pb_z[k] * gk_189[k];

        t_422[k] = f_15 * fk_210[k]
                   + pa_x[k] * fl_149[k];

        t_423[k] = f_15 * fk_211[k]
                   + pa_x[k] * fl_150[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pa_x, pb_y, pb_z, fk_108, fk_213, fk_214, \
                         fl_151, fl_152, gk_190, gk_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_y[k] * gk_190[k];

        t_425[k] = f_15 * fk_213[k]
                   + pa_x[k] * fl_151[k];

        t_426[k] = f_14 * fk_214[k]
                   + pa_x[k] * fl_152[k];

        t_427[k] = f_15 * fk_108[k]
                   + pb_z[k] * gk_191[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pa_x, pb_y, fk_215, fk_216, \
                         fk_217, fk_218, fl_153, fl_154, fl_155, fl_156, \
                         gk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_14 * fk_215[k]
                   + pa_x[k] * fl_153[k];

        t_429[k] = f_14 * fk_216[k]
                   + pa_x[k] * fl_154[k];

        t_430[k] = f_14 * fk_217[k]
                   + pa_x[k] * fl_155[k];

        t_431[k] = pb_y[k] * gk_192[k];

        t_432[k] = f_14 * fk_218[k]
                   + pa_x[k] * fl_156[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, pb_x, fk_219, fk_220, fk_221, \
                         fk_222, fk_223, gk_194, gk_195, gk_196, gk_197, \
                         gk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_13 * fk_219[k]
                   + pb_x[k] * gk_194[k];

        t_434[k] = f_13 * fk_220[k]
                   + pb_x[k] * gk_195[k];

        t_435[k] = f_13 * fk_221[k]
                   + pb_x[k] * gk_196[k];

        t_436[k] = f_13 * fk_222[k]
                   + pb_x[k] * gk_197[k];

        t_437[k] = f_13 * fk_223[k]
                   + pb_x[k] * gk_198[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pa_x, pb_x, pb_y, fk_224, fk_226, \
                         fl_157, fl_158, gk_193, gk_199, gk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_13 * fk_224[k]
                   + pb_x[k] * gk_199[k];

        t_439[k] = pb_y[k] * gk_193[k];

        t_440[k] = f_13 * fk_226[k]
                   + pb_x[k] * gk_200[k];

        t_441[k] = pa_x[k] * fl_157[k];

        t_442[k] = pa_x[k] * fl_158[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, t_449, pa_x, pb_y, fl_159, \
                         fl_160, fl_161, fl_162, fl_163, fl_164, \
                         gk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pa_x[k] * fl_159[k];

        t_444[k] = pa_x[k] * fl_160[k];

        t_445[k] = pa_x[k] * fl_161[k];

        t_446[k] = pa_x[k] * fl_162[k];

        t_447[k] = pa_x[k] * fl_163[k];

        t_448[k] = pb_y[k] * gk_200[k];

        t_449[k] = pa_x[k] * fl_164[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, pb_x, pb_y, pb_z, fk_118, gi0_54, \
                         gi0_55, gi1_54, gi1_55, gk_201, gk_202, \
                         gk_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * gi0_54[k]
                   - f_2 * gi1_54[k]
                   + pb_x[k] * gk_201[k];

        t_451[k] = f_0 * fk_118[k]
                   + pb_y[k] * gk_201[k];

        t_452[k] = pb_z[k] * gk_201[k];

        t_453[k] = f_11 * gi0_55[k]
                   - f_12 * gi1_55[k]
                   + pb_x[k] * gk_203[k];

        t_454[k] = pb_z[k] * gk_202[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pb_x, pb_y, pb_z, fk_121, gi0_56, gi0_57, \
                         gi1_56, gi1_57, gk_203, gk_204, gk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * gi0_56[k]
                   - f_12 * gi1_56[k]
                   + pb_x[k] * gk_204[k];

        t_456[k] = f_9 * gi0_57[k]
                   - f_10 * gi1_57[k]
                   + pb_x[k] * gk_205[k];

        t_457[k] = pb_z[k] * gk_203[k];

        t_458[k] = f_0 * fk_121[k]
                   + pb_y[k] * gk_204[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, pb_x, pb_z, gi0_58, gi0_59, gi0_60, \
                         gi1_58, gi1_59, gi1_60, gk_205, gk_206, gk_207, \
                         gk_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_9 * gi0_58[k]
                   - f_10 * gi1_58[k]
                   + pb_x[k] * gk_206[k];

        t_460[k] = f_7 * gi0_59[k]
                   - f_8 * gi1_59[k]
                   + pb_x[k] * gk_207[k];

        t_461[k] = pb_z[k] * gk_205[k];

        t_462[k] = f_7 * gi0_60[k]
                   - f_8 * gi1_60[k]
                   + pb_x[k] * gk_208[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pb_x, pb_y, pb_z, fk_124, gi0_61, gi0_62, \
                         gi1_61, gi1_62, gk_206, gk_207, gk_209, \
                         gk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_0 * fk_124[k]
                   + pb_y[k] * gk_206[k];

        t_464[k] = f_7 * gi0_61[k]
                   - f_8 * gi1_61[k]
                   + pb_x[k] * gk_209[k];

        t_465[k] = f_5 * gi0_62[k]
                   - f_6 * gi1_62[k]
                   + pb_x[k] * gk_210[k];

        t_466[k] = pb_z[k] * gk_207[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_x, pb_y, fk_128, gi0_63, gi0_64, gi1_63, \
                         gi1_64, gk_209, gk_211, gk_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_5 * gi0_63[k]
                   - f_6 * gi1_63[k]
                   + pb_x[k] * gk_211[k];

        t_468[k] = f_5 * gi0_64[k]
                   - f_6 * gi1_64[k]
                   + pb_x[k] * gk_212[k];

        t_469[k] = f_0 * fk_128[k]
                   + pb_y[k] * gk_209[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_x, pb_z, gi0_65, gi0_66, gi0_68, \
                         gi1_65, gi1_66, gi1_68, gk_210, gk_213, gk_214, \
                         gk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_5 * gi0_65[k]
                   - f_6 * gi1_65[k]
                   + pb_x[k] * gk_213[k];

        t_471[k] = f_3 * gi0_66[k]
                   - f_4 * gi1_66[k]
                   + pb_x[k] * gk_214[k];

        t_472[k] = pb_z[k] * gk_210[k];

        t_473[k] = f_3 * gi0_68[k]
                   - f_4 * gi1_68[k]
                   + pb_x[k] * gk_215[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pb_x, pb_y, fk_133, gi0_69, gi0_70, gi1_69, \
                         gi1_70, gk_213, gk_216, gk_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_3 * gi0_69[k]
                   - f_4 * gi1_69[k]
                   + pb_x[k] * gk_216[k];

        t_475[k] = f_3 * gi0_70[k]
                   - f_4 * gi1_70[k]
                   + pb_x[k] * gk_217[k];

        t_476[k] = f_0 * fk_133[k]
                   + pb_y[k] * gk_213[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, t_482, pb_x, gi0_71, gi1_71, \
                         gk_218, gk_219, gk_220, gk_221, gk_222, \
                         gk_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_3 * gi0_71[k]
                   - f_4 * gi1_71[k]
                   + pb_x[k] * gk_218[k];

        t_478[k] = pb_x[k] * gk_219[k];

        t_479[k] = pb_x[k] * gk_220[k];

        t_480[k] = pb_x[k] * gk_221[k];

        t_481[k] = pb_x[k] * gk_222[k];

        t_482[k] = pb_x[k] * gk_223[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, pb_x, pb_y, pb_z, fk_139, gi0_66, \
                         gi1_66, gk_219, gk_224, gk_225, gk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pb_x[k] * gk_224[k];

        t_484[k] = pb_x[k] * gk_225[k];

        t_485[k] = pb_x[k] * gk_226[k];

        t_486[k] = f_0 * fk_139[k]
                   + f_1 * gi0_66[k]
                   - f_2 * gi1_66[k]
                   + pb_y[k] * gk_219[k];

        t_487[k] = pb_z[k] * gk_219[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pb_z, gi0_66, gi0_67, gi0_68, gi1_66, gi1_67, \
                         gi1_68, gk_220, gk_221, gk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_3 * gi0_66[k]
                   - f_4 * gi1_66[k]
                   + pb_z[k] * gk_220[k];

        t_489[k] = f_5 * gi0_67[k]
                   - f_6 * gi1_67[k]
                   + pb_z[k] * gk_221[k];

        t_490[k] = f_7 * gi0_68[k]
                   - f_8 * gi1_68[k]
                   + pb_z[k] * gk_222[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pb_y, pb_z, fk_146, gi0_69, gi0_70, \
                         gi0_71, gi1_69, gi1_70, gi1_71, gk_223, gk_224, \
                         gk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_9 * gi0_69[k]
                   - f_10 * gi1_69[k]
                   + pb_z[k] * gk_223[k];

        t_492[k] = f_11 * gi0_70[k]
                   - f_12 * gi1_70[k]
                   + pb_z[k] * gk_224[k];

        t_493[k] = f_0 * fk_146[k]
                   + pb_y[k] * gk_226[k];

        t_494[k] = f_1 * gi0_71[k]
                   - f_2 * gi1_71[k]
                   + pb_z[k] * gk_226[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, pa_z, pb_y, pb_z, fk_118, fk_148, \
                         fl_73, fl_74, fl_75, gk_227, gk_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * fl_73[k];

        t_496[k] = pa_z[k] * fl_74[k];

        t_497[k] = f_13 * fk_118[k]
                   + pb_z[k] * gk_227[k];

        t_498[k] = pa_z[k] * fl_75[k];

        t_499[k] = f_15 * fk_148[k]
                   + pb_y[k] * gk_228[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_z, pb_y, pb_z, fk_119, fk_120, fk_150, \
                         fl_76, fl_77, gk_229, gk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_14 * fk_119[k]
                   + pa_z[k] * fl_76[k];

        t_501[k] = pa_z[k] * fl_77[k];

        t_502[k] = f_13 * fk_120[k]
                   + pb_z[k] * gk_229[k];

        t_503[k] = f_15 * fk_150[k]
                   + pb_y[k] * gk_230[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_z, pb_z, fk_121, fk_122, fk_123, \
                         fl_78, fl_79, fl_80, gk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_15 * fk_121[k]
                   + pa_z[k] * fl_78[k];

        t_505[k] = pa_z[k] * fl_79[k];

        t_506[k] = f_13 * fk_122[k]
                   + pb_z[k] * gk_231[k];

        t_507[k] = f_14 * fk_123[k]
                   + pa_z[k] * fl_80[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_z, pb_y, pb_z, fk_124, fk_125, fk_152, \
                         fl_81, fl_82, gk_232, gk_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_15 * fk_152[k]
                   + pb_y[k] * gk_232[k];

        t_509[k] = f_0 * fk_124[k]
                   + pa_z[k] * fl_81[k];

        t_510[k] = pa_z[k] * fl_82[k];

        t_511[k] = f_13 * fk_125[k]
                   + pb_z[k] * gk_233[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pa_z, pb_y, fk_126, fk_127, \
                         fk_128, fk_155, fl_83, fl_84, fl_85, fl_86, \
                         gk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_14 * fk_126[k]
                   + pa_z[k] * fl_83[k];

        t_513[k] = f_15 * fk_127[k]
                   + pa_z[k] * fl_84[k];

        t_514[k] = f_15 * fk_155[k]
                   + pb_y[k] * gk_234[k];

        t_515[k] = f_16 * fk_128[k]
                   + pa_z[k] * fl_85[k];

        t_516[k] = pa_z[k] * fl_86[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_z, pb_z, fk_129, fk_130, fk_131, \
                         fk_132, fl_87, fl_88, fl_89, gk_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_13 * fk_129[k]
                   + pb_z[k] * gk_235[k];

        t_518[k] = f_14 * fk_130[k]
                   + pa_z[k] * fl_87[k];

        t_519[k] = f_15 * fk_131[k]
                   + pa_z[k] * fl_88[k];

        t_520[k] = f_0 * fk_132[k]
                   + pa_z[k] * fl_89[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, pa_z, pb_x, pb_y, fk_133, fk_159, \
                         fl_90, gk_236, gk_237, gk_238, gk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_15 * fk_159[k]
                   + pb_y[k] * gk_236[k];

        t_522[k] = f_17 * fk_133[k]
                   + pa_z[k] * fl_90[k];

        t_523[k] = pb_x[k] * gk_237[k];

        t_524[k] = pb_x[k] * gk_238[k];

        t_525[k] = pb_x[k] * gk_239[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, t_531, pa_z, pb_x, fl_91, gk_240, \
                         gk_241, gk_242, gk_243, gk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = pb_x[k] * gk_240[k];

        t_527[k] = pb_x[k] * gk_241[k];

        t_528[k] = pb_x[k] * gk_242[k];

        t_529[k] = pb_x[k] * gk_243[k];

        t_530[k] = pb_x[k] * gk_244[k];

        t_531[k] = pa_z[k] * fl_91[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pa_z, pb_z, fk_139, fk_140, fk_141, \
                         fk_142, fl_92, fl_93, fl_94, gk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_13 * fk_139[k]
                   + pb_z[k] * gk_237[k];

        t_533[k] = f_14 * fk_140[k]
                   + pa_z[k] * fl_92[k];

        t_534[k] = f_15 * fk_141[k]
                   + pa_z[k] * fl_93[k];

        t_535[k] = f_0 * fk_142[k]
                   + pa_z[k] * fl_94[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_z, pb_y, fk_143, fk_144, fk_146, \
                         fk_171, fl_95, fl_96, fl_98, gk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_16 * fk_143[k]
                   + pa_z[k] * fl_95[k];

        t_537[k] = f_17 * fk_144[k]
                   + pa_z[k] * fl_96[k];

        t_538[k] = f_15 * fk_171[k]
                   + pb_y[k] * gk_244[k];

        t_539[k] = f_18 * fk_146[k]
                   + pa_z[k] * fl_98[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pb_x, pb_y, pb_z, fk_147, fk_172, gi0_72, \
                         gi0_73, gi1_72, gi1_73, gk_245, gk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * gi0_72[k]
                   - f_2 * gi1_72[k]
                   + pb_x[k] * gk_245[k];

        t_541[k] = f_14 * fk_172[k]
                   + pb_y[k] * gk_245[k];

        t_542[k] = f_14 * fk_147[k]
                   + pb_z[k] * gk_245[k];

        t_543[k] = f_11 * gi0_73[k]
                   - f_12 * gi1_73[k]
                   + pb_x[k] * gk_247[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pb_x, pb_y, fk_173, gi0_74, gi0_75, gi1_74, \
                         gi1_75, gk_246, gk_248, gk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_14 * fk_173[k]
                   + pb_y[k] * gk_246[k];

        t_545[k] = f_11 * gi0_74[k]
                   - f_12 * gi1_74[k]
                   + pb_x[k] * gk_248[k];

        t_546[k] = f_9 * gi0_75[k]
                   - f_10 * gi1_75[k]
                   + pb_x[k] * gk_249[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pb_x, pb_y, pb_z, fk_149, fk_175, gi0_76, \
                         gi1_76, gk_247, gk_248, gk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_14 * fk_149[k]
                   + pb_z[k] * gk_247[k];

        t_548[k] = f_14 * fk_175[k]
                   + pb_y[k] * gk_248[k];

        t_549[k] = f_9 * gi0_76[k]
                   - f_10 * gi1_76[k]
                   + pb_x[k] * gk_250[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pb_x, pb_z, fk_151, gi0_77, gi0_78, gi1_77, \
                         gi1_78, gk_249, gk_251, gk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_7 * gi0_77[k]
                   - f_8 * gi1_77[k]
                   + pb_x[k] * gk_251[k];

        t_551[k] = f_14 * fk_151[k]
                   + pb_z[k] * gk_249[k];

        t_552[k] = f_7 * gi0_78[k]
                   - f_8 * gi1_78[k]
                   + pb_x[k] * gk_252[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pb_x, pb_y, fk_177, gi0_79, gi0_80, gi1_79, \
                         gi1_80, gk_250, gk_253, gk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_14 * fk_177[k]
                   + pb_y[k] * gk_250[k];

        t_554[k] = f_7 * gi0_79[k]
                   - f_8 * gi1_79[k]
                   + pb_x[k] * gk_253[k];

        t_555[k] = f_5 * gi0_80[k]
                   - f_6 * gi1_80[k]
                   + pb_x[k] * gk_254[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pb_x, pb_z, fk_153, gi0_81, gi0_82, gi1_81, \
                         gi1_82, gk_251, gk_255, gk_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_14 * fk_153[k]
                   + pb_z[k] * gk_251[k];

        t_557[k] = f_5 * gi0_81[k]
                   - f_6 * gi1_81[k]
                   + pb_x[k] * gk_255[k];

        t_558[k] = f_5 * gi0_82[k]
                   - f_6 * gi1_82[k]
                   + pb_x[k] * gk_256[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, pb_x, pb_y, fk_180, gi0_83, gi0_84, gi1_83, \
                         gi1_84, gk_253, gk_257, gk_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_14 * fk_180[k]
                   + pb_y[k] * gk_253[k];

        t_560[k] = f_5 * gi0_83[k]
                   - f_6 * gi1_83[k]
                   + pb_x[k] * gk_257[k];

        t_561[k] = f_3 * gi0_84[k]
                   - f_4 * gi1_84[k]
                   + pb_x[k] * gk_258[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, pb_x, pb_z, fk_156, gi0_85, gi0_86, gi1_85, \
                         gi1_86, gk_254, gk_259, gk_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_14 * fk_156[k]
                   + pb_z[k] * gk_254[k];

        t_563[k] = f_3 * gi0_85[k]
                   - f_4 * gi1_85[k]
                   + pb_x[k] * gk_259[k];

        t_564[k] = f_3 * gi0_86[k]
                   - f_4 * gi1_86[k]
                   + pb_x[k] * gk_260[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pb_x, pb_y, fk_184, gi0_87, gi0_89, \
                         gi1_87, gi1_89, gk_257, gk_261, gk_262, \
                         gk_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_3 * gi0_87[k]
                   - f_4 * gi1_87[k]
                   + pb_x[k] * gk_261[k];

        t_566[k] = f_14 * fk_184[k]
                   + pb_y[k] * gk_257[k];

        t_567[k] = f_3 * gi0_89[k]
                   - f_4 * gi1_89[k]
                   + pb_x[k] * gk_262[k];

        t_568[k] = pb_x[k] * gk_263[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, t_574, t_575, pb_x, gk_264, \
                         gk_265, gk_266, gk_267, gk_268, gk_269, \
                         gk_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = pb_x[k] * gk_264[k];

        t_570[k] = pb_x[k] * gk_265[k];

        t_571[k] = pb_x[k] * gk_266[k];

        t_572[k] = pb_x[k] * gk_267[k];

        t_573[k] = pb_x[k] * gk_268[k];

        t_574[k] = pb_x[k] * gk_269[k];

        t_575[k] = pb_x[k] * gk_270[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, pa_z, pb_y, pb_z, dl0_1, dl1_1, fk_164, fk_191, \
                         fl_110, gi0_85, gi1_85, gk_263, gk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_19 * dl0_1[k]
                   - f_20 * dl1_1[k]
                   + pa_z[k] * fl_110[k];

        t_577[k] = f_14 * fk_164[k]
                   + pb_z[k] * gk_263[k];

        t_578[k] = f_14 * fk_191[k]
                   + f_11 * gi0_85[k]
                   - f_12 * gi1_85[k]
                   + pb_y[k] * gk_265[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pb_y, fk_192, fk_193, fk_194, gi0_86, gi0_87, \
                         gi0_88, gi1_86, gi1_87, gi1_88, gk_266, gk_267, \
                         gk_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_14 * fk_192[k]
                   + f_9 * gi0_86[k]
                   - f_10 * gi1_86[k]
                   + pb_y[k] * gk_266[k];

        t_580[k] = f_14 * fk_193[k]
                   + f_7 * gi0_87[k]
                   - f_8 * gi1_87[k]
                   + pb_y[k] * gk_267[k];

        t_581[k] = f_14 * fk_194[k]
                   + f_5 * gi0_88[k]
                   - f_6 * gi1_88[k]
                   + pb_y[k] * gk_268[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pa_y, pb_y, dl0_2, dl1_2, fk_195, fk_196, \
                         fl_138, fl_139, gi0_89, gi1_89, gk_269, \
                         gk_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_14 * fk_195[k]
                   + f_3 * gi0_89[k]
                   - f_4 * gi1_89[k]
                   + pb_y[k] * gk_269[k];

        t_583[k] = f_14 * fk_196[k]
                   + pb_y[k] * gk_270[k];

        t_584[k] = f_19 * dl0_2[k]
                   - f_20 * dl1_2[k]
                   + pa_y[k] * fl_138[k];

        t_585[k] = pa_y[k] * fl_139[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, pa_y, pb_y, fk_197, fk_198, \
                         fk_199, fl_140, fl_141, fl_142, gk_271, \
                         gk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_13 * fk_197[k]
                   + pb_y[k] * gk_271[k];

        t_587[k] = pa_y[k] * fl_140[k];

        t_588[k] = f_14 * fk_198[k]
                   + pa_y[k] * fl_141[k];

        t_589[k] = f_13 * fk_199[k]
                   + pb_y[k] * gk_272[k];

        t_590[k] = pa_y[k] * fl_142[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pa_y, pb_y, pb_z, fk_174, fk_200, fk_201, \
                         fl_143, fl_144, gk_273, gk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_15 * fk_200[k]
                   + pa_y[k] * fl_143[k];

        t_592[k] = f_15 * fk_174[k]
                   + pb_z[k] * gk_273[k];

        t_593[k] = f_13 * fk_201[k]
                   + pb_y[k] * gk_274[k];

        t_594[k] = pa_y[k] * fl_144[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pa_y, pb_y, pb_z, fk_176, fk_202, fk_203, \
                         fk_204, fl_145, fl_146, gk_275, gk_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_0 * fk_202[k]
                   + pa_y[k] * fl_145[k];

        t_596[k] = f_15 * fk_176[k]
                   + pb_z[k] * gk_275[k];

        t_597[k] = f_14 * fk_203[k]
                   + pa_y[k] * fl_146[k];

        t_598[k] = f_13 * fk_204[k]
                   + pb_y[k] * gk_276[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pa_y, pb_z, fk_178, fk_205, \
                         fk_206, fk_207, fl_147, fl_148, fl_149, fl_150, \
                         gk_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pa_y[k] * fl_147[k];

        t_600[k] = f_16 * fk_205[k]
                   + pa_y[k] * fl_148[k];

        t_601[k] = f_15 * fk_178[k]
                   + pb_z[k] * gk_277[k];

        t_602[k] = f_15 * fk_206[k]
                   + pa_y[k] * fl_149[k];

        t_603[k] = f_14 * fk_207[k]
                   + pa_y[k] * fl_150[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_y, pb_y, pb_z, fk_181, fk_208, fk_209, \
                         fl_151, fl_152, gk_278, gk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_13 * fk_208[k]
                   + pb_y[k] * gk_278[k];

        t_605[k] = pa_y[k] * fl_151[k];

        t_606[k] = f_17 * fk_209[k]
                   + pa_y[k] * fl_152[k];

        t_607[k] = f_15 * fk_181[k]
                   + pb_z[k] * gk_279[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, t_612, pa_y, pb_y, fk_210, fk_211, \
                         fk_212, fk_213, fl_153, fl_154, fl_155, fl_156, \
                         gk_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_0 * fk_210[k]
                   + pa_y[k] * fl_153[k];

        t_609[k] = f_15 * fk_211[k]
                   + pa_y[k] * fl_154[k];

        t_610[k] = f_14 * fk_212[k]
                   + pa_y[k] * fl_155[k];

        t_611[k] = f_13 * fk_213[k]
                   + pb_y[k] * gk_280[k];

        t_612[k] = pa_y[k] * fl_156[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, t_617, t_618, t_619, pb_x, gk_281, \
                         gk_282, gk_283, gk_284, gk_285, gk_286, \
                         gk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = pb_x[k] * gk_281[k];

        t_614[k] = pb_x[k] * gk_282[k];

        t_615[k] = pb_x[k] * gk_283[k];

        t_616[k] = pb_x[k] * gk_284[k];

        t_617[k] = pb_x[k] * gk_285[k];

        t_618[k] = pb_x[k] * gk_286[k];

        t_619[k] = pb_x[k] * gk_287[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pa_y, pb_x, pb_z, fk_189, fk_219, fk_221, \
                         fl_157, fl_159, gk_281, gk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = pb_x[k] * gk_288[k];

        t_621[k] = f_18 * fk_219[k]
                   + pa_y[k] * fl_157[k];

        t_622[k] = f_15 * fk_189[k]
                   + pb_z[k] * gk_281[k];

        t_623[k] = f_17 * fk_221[k]
                   + pa_y[k] * fl_159[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pa_y, fk_222, fk_223, fk_224, fk_225, \
                         fl_160, fl_161, fl_162, fl_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_16 * fk_222[k]
                   + pa_y[k] * fl_160[k];

        t_625[k] = f_0 * fk_223[k]
                   + pa_y[k] * fl_161[k];

        t_626[k] = f_15 * fk_224[k]
                   + pa_y[k] * fl_162[k];

        t_627[k] = f_14 * fk_225[k]
                   + pa_y[k] * fl_163[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, t_632, pa_y, pb_x, pb_y, pb_z, fk_197, \
                         fk_226, fl_164, gi0_90, gi1_90, gk_288, \
                         gk_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_13 * fk_226[k]
                   + pb_y[k] * gk_288[k];

        t_629[k] = pa_y[k] * fl_164[k];

        t_630[k] = f_1 * gi0_90[k]
                   - f_2 * gi1_90[k]
                   + pb_x[k] * gk_289[k];

        t_631[k] = pb_y[k] * gk_289[k];

        t_632[k] = f_0 * fk_197[k]
                   + pb_z[k] * gk_289[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pb_x, pb_y, gi0_91, gi0_92, gi0_93, \
                         gi1_91, gi1_92, gi1_93, gk_290, gk_291, gk_292, \
                         gk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_11 * gi0_91[k]
                   - f_12 * gi1_91[k]
                   + pb_x[k] * gk_291[k];

        t_634[k] = pb_y[k] * gk_290[k];

        t_635[k] = f_11 * gi0_92[k]
                   - f_12 * gi1_92[k]
                   + pb_x[k] * gk_292[k];

        t_636[k] = f_9 * gi0_93[k]
                   - f_10 * gi1_93[k]
                   + pb_x[k] * gk_293[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pb_x, pb_y, pb_z, fk_200, gi0_94, gi0_95, \
                         gi1_94, gi1_95, gk_291, gk_292, gk_294, \
                         gk_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = f_0 * fk_200[k]
                   + pb_z[k] * gk_291[k];

        t_638[k] = pb_y[k] * gk_292[k];

        t_639[k] = f_9 * gi0_94[k]
                   - f_10 * gi1_94[k]
                   + pb_x[k] * gk_294[k];

        t_640[k] = f_7 * gi0_95[k]
                   - f_8 * gi1_95[k]
                   + pb_x[k] * gk_295[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, pb_x, pb_y, pb_z, fk_202, gi0_96, gi0_97, \
                         gi1_96, gi1_97, gk_293, gk_294, gk_296, \
                         gk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_0 * fk_202[k]
                   + pb_z[k] * gk_293[k];

        t_642[k] = f_7 * gi0_96[k]
                   - f_8 * gi1_96[k]
                   + pb_x[k] * gk_296[k];

        t_643[k] = pb_y[k] * gk_294[k];

        t_644[k] = f_7 * gi0_97[k]
                   - f_8 * gi1_97[k]
                   + pb_x[k] * gk_297[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pb_x, pb_z, fk_205, gi0_98, gi0_99, gi1_98, \
                         gi1_99, gk_295, gk_298, gk_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_5 * gi0_98[k]
                   - f_6 * gi1_98[k]
                   + pb_x[k] * gk_298[k];

        t_646[k] = f_0 * fk_205[k]
                   + pb_z[k] * gk_295[k];

        t_647[k] = f_5 * gi0_99[k]
                   - f_6 * gi1_99[k]
                   + pb_x[k] * gk_299[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, t_651, pb_x, pb_y, gi0_100, gi0_101, gi0_102, \
                         gi1_100, gi1_101, gi1_102, gk_297, gk_300, gk_301, \
                         gk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_5 * gi0_100[k]
                   - f_6 * gi1_100[k]
                   + pb_x[k] * gk_300[k];

        t_649[k] = pb_y[k] * gk_297[k];

        t_650[k] = f_5 * gi0_101[k]
                   - f_6 * gi1_101[k]
                   + pb_x[k] * gk_301[k];

        t_651[k] = f_3 * gi0_102[k]
                   - f_4 * gi1_102[k]
                   + pb_x[k] * gk_302[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pb_x, pb_z, fk_209, gi0_103, gi0_104, gi1_103, \
                         gi1_104, gk_298, gk_303, gk_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_0 * fk_209[k]
                   + pb_z[k] * gk_298[k];

        t_653[k] = f_3 * gi0_103[k]
                   - f_4 * gi1_103[k]
                   + pb_x[k] * gk_303[k];

        t_654[k] = f_3 * gi0_104[k]
                   - f_4 * gi1_104[k]
                   + pb_x[k] * gk_304[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, pb_x, pb_y, gi0_105, gi0_107, \
                         gi1_105, gi1_107, gk_301, gk_305, gk_306, gk_307, \
                         gk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_3 * gi0_105[k]
                   - f_4 * gi1_105[k]
                   + pb_x[k] * gk_305[k];

        t_656[k] = pb_y[k] * gk_301[k];

        t_657[k] = f_3 * gi0_107[k]
                   - f_4 * gi1_107[k]
                   + pb_x[k] * gk_306[k];

        t_658[k] = pb_x[k] * gk_307[k];

        t_659[k] = pb_x[k] * gk_308[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, pb_x, gk_309, gk_310, \
                         gk_311, gk_312, gk_313, gk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = pb_x[k] * gk_309[k];

        t_661[k] = pb_x[k] * gk_310[k];

        t_662[k] = pb_x[k] * gk_311[k];

        t_663[k] = pb_x[k] * gk_312[k];

        t_664[k] = pb_x[k] * gk_313[k];

        t_665[k] = pb_x[k] * gk_314[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_y, pb_z, fk_219, gi0_102, gi0_103, \
                         gi0_104, gi1_102, gi1_103, gi1_104, gk_307, gk_309, \
                         gk_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * gi0_102[k]
                   - f_2 * gi1_102[k]
                   + pb_y[k] * gk_307[k];

        t_667[k] = f_0 * fk_219[k]
                   + pb_z[k] * gk_307[k];

        t_668[k] = f_11 * gi0_103[k]
                   - f_12 * gi1_103[k]
                   + pb_y[k] * gk_309[k];

        t_669[k] = f_9 * gi0_104[k]
                   - f_10 * gi1_104[k]
                   + pb_y[k] * gk_310[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_y, gi0_105, gi0_106, gi0_107, gi1_105, \
                         gi1_106, gi1_107, gk_311, gk_312, gk_313, \
                         gk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_7 * gi0_105[k]
                   - f_8 * gi1_105[k]
                   + pb_y[k] * gk_311[k];

        t_671[k] = f_5 * gi0_106[k]
                   - f_6 * gi1_106[k]
                   + pb_y[k] * gk_312[k];

        t_672[k] = f_3 * gi0_107[k]
                   - f_4 * gi1_107[k]
                   + pb_y[k] * gk_313[k];

        t_673[k] = pb_y[k] * gk_314[k];
    }

#pragma omp simd aligned(t_674, pb_z, fk_226, gi0_107, gi1_107, \
                         gk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_0 * fk_226[k]
                   + f_1 * gi0_107[k]
                   - f_2 * gi1_107[k]
                   + pb_z[k] * gk_314[k];
    }
}

auto
compute_prim_gl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

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
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
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
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
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
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);
    const auto *fl_9 = buffer.data(fl + 9);
    const auto *fl_10 = buffer.data(fl + 10);
    const auto *fl_11 = buffer.data(fl + 11);
    const auto *fl_12 = buffer.data(fl + 12);
    const auto *fl_13 = buffer.data(fl + 13);
    const auto *fl_14 = buffer.data(fl + 14);
    const auto *fl_15 = buffer.data(fl + 15);
    const auto *fl_16 = buffer.data(fl + 16);
    const auto *fl_17 = buffer.data(fl + 17);
    const auto *fl_18 = buffer.data(fl + 18);
    const auto *fl_19 = buffer.data(fl + 19);
    const auto *fl_20 = buffer.data(fl + 20);
    const auto *fl_21 = buffer.data(fl + 21);
    const auto *fl_22 = buffer.data(fl + 22);
    const auto *fl_23 = buffer.data(fl + 23);
    const auto *fl_24 = buffer.data(fl + 24);
    const auto *fl_25 = buffer.data(fl + 25);
    const auto *fl_26 = buffer.data(fl + 26);
    const auto *fl_27 = buffer.data(fl + 27);
    const auto *fl_28 = buffer.data(fl + 28);
    const auto *fl_29 = buffer.data(fl + 29);
    const auto *fl_30 = buffer.data(fl + 30);
    const auto *fl_31 = buffer.data(fl + 31);
    const auto *fl_32 = buffer.data(fl + 32);
    const auto *fl_33 = buffer.data(fl + 33);
    const auto *fl_34 = buffer.data(fl + 34);
    const auto *fl_35 = buffer.data(fl + 35);
    const auto *fl_36 = buffer.data(fl + 36);
    const auto *fl_37 = buffer.data(fl + 37);
    const auto *fl_38 = buffer.data(fl + 38);
    const auto *fl_39 = buffer.data(fl + 39);
    const auto *fl_40 = buffer.data(fl + 40);
    const auto *fl_41 = buffer.data(fl + 41);
    const auto *fl_42 = buffer.data(fl + 42);
    const auto *fl_43 = buffer.data(fl + 43);
    const auto *fl_44 = buffer.data(fl + 44);
    const auto *fl_45 = buffer.data(fl + 45);
    const auto *fl_46 = buffer.data(fl + 46);
    const auto *fl_47 = buffer.data(fl + 47);
    const auto *fl_48 = buffer.data(fl + 48);
    const auto *fl_49 = buffer.data(fl + 49);
    const auto *fl_50 = buffer.data(fl + 50);
    const auto *fl_51 = buffer.data(fl + 51);
    const auto *fl_52 = buffer.data(fl + 52);
    const auto *fl_53 = buffer.data(fl + 53);
    const auto *fl_54 = buffer.data(fl + 54);
    const auto *fl_55 = buffer.data(fl + 55);
    const auto *fl_56 = buffer.data(fl + 56);
    const auto *fl_57 = buffer.data(fl + 57);
    const auto *fl_58 = buffer.data(fl + 58);
    const auto *fl_59 = buffer.data(fl + 59);
    const auto *fl_60 = buffer.data(fl + 60);
    const auto *fl_61 = buffer.data(fl + 61);
    const auto *fl_62 = buffer.data(fl + 62);
    const auto *fl_63 = buffer.data(fl + 63);
    const auto *fl_64 = buffer.data(fl + 64);
    const auto *fl_65 = buffer.data(fl + 65);
    const auto *fl_66 = buffer.data(fl + 66);
    const auto *fl_67 = buffer.data(fl + 67);
    const auto *fl_68 = buffer.data(fl + 68);
    const auto *fl_69 = buffer.data(fl + 69);
    const auto *fl_70 = buffer.data(fl + 70);
    const auto *fl_71 = buffer.data(fl + 71);
    const auto *fl_72 = buffer.data(fl + 72);
    const auto *fl_73 = buffer.data(fl + 73);
    const auto *fl_74 = buffer.data(fl + 74);
    const auto *fl_75 = buffer.data(fl + 75);
    const auto *fl_76 = buffer.data(fl + 76);
    const auto *fl_77 = buffer.data(fl + 77);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_1 = buffer.data(gi0 + 1);
    const auto *gi0_2 = buffer.data(gi0 + 2);
    const auto *gi0_3 = buffer.data(gi0 + 3);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_9 = buffer.data(gi0 + 9);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_15 = buffer.data(gi0 + 15);
    const auto *gi0_16 = buffer.data(gi0 + 16);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_25 = buffer.data(gi0 + 25);
    const auto *gi0_26 = buffer.data(gi0 + 26);
    const auto *gi0_27 = buffer.data(gi0 + 27);
    const auto *gi0_28 = buffer.data(gi0 + 28);
    const auto *gi0_29 = buffer.data(gi0 + 29);
    const auto *gi0_30 = buffer.data(gi0 + 30);
    const auto *gi0_31 = buffer.data(gi0 + 31);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_33 = buffer.data(gi0 + 33);
    const auto *gi0_34 = buffer.data(gi0 + 34);
    const auto *gi0_35 = buffer.data(gi0 + 35);
    const auto *gi0_36 = buffer.data(gi0 + 36);
    const auto *gi0_37 = buffer.data(gi0 + 37);
    const auto *gi0_38 = buffer.data(gi0 + 38);
    const auto *gi0_39 = buffer.data(gi0 + 39);
    const auto *gi0_40 = buffer.data(gi0 + 40);
    const auto *gi0_41 = buffer.data(gi0 + 41);
    const auto *gi0_42 = buffer.data(gi0 + 42);
    const auto *gi0_43 = buffer.data(gi0 + 43);
    const auto *gi0_44 = buffer.data(gi0 + 44);
    const auto *gi0_45 = buffer.data(gi0 + 45);
    const auto *gi0_46 = buffer.data(gi0 + 46);
    const auto *gi0_47 = buffer.data(gi0 + 47);
    const auto *gi0_48 = buffer.data(gi0 + 48);
    const auto *gi0_49 = buffer.data(gi0 + 49);
    const auto *gi0_50 = buffer.data(gi0 + 50);
    const auto *gi0_51 = buffer.data(gi0 + 51);
    const auto *gi0_52 = buffer.data(gi0 + 52);
    const auto *gi0_53 = buffer.data(gi0 + 53);
    const auto *gi0_54 = buffer.data(gi0 + 54);
    const auto *gi0_55 = buffer.data(gi0 + 55);
    const auto *gi0_58 = buffer.data(gi0 + 58);
    const auto *gi0_59 = buffer.data(gi0 + 59);
    const auto *gi0_60 = buffer.data(gi0 + 60);
    const auto *gi0_61 = buffer.data(gi0 + 61);
    const auto *gi0_62 = buffer.data(gi0 + 62);
    const auto *gi0_63 = buffer.data(gi0 + 63);
    const auto *gi0_64 = buffer.data(gi0 + 64);
    const auto *gi0_65 = buffer.data(gi0 + 65);
    const auto *gi0_66 = buffer.data(gi0 + 66);
    const auto *gi0_67 = buffer.data(gi0 + 67);
    const auto *gi0_68 = buffer.data(gi0 + 68);
    const auto *gi0_69 = buffer.data(gi0 + 69);
    const auto *gi0_70 = buffer.data(gi0 + 70);
    const auto *gi0_71 = buffer.data(gi0 + 71);
    const auto *gi0_72 = buffer.data(gi0 + 72);
    const auto *gi0_73 = buffer.data(gi0 + 73);
    const auto *gi0_74 = buffer.data(gi0 + 74);
    const auto *gi0_75 = buffer.data(gi0 + 75);
    const auto *gi0_77 = buffer.data(gi0 + 77);
    const auto *gi0_78 = buffer.data(gi0 + 78);
    const auto *gi0_79 = buffer.data(gi0 + 79);
    const auto *gi0_80 = buffer.data(gi0 + 80);
    const auto *gi0_81 = buffer.data(gi0 + 81);
    const auto *gi0_82 = buffer.data(gi0 + 82);
    const auto *gi0_83 = buffer.data(gi0 + 83);
    const auto *gi0_84 = buffer.data(gi0 + 84);
    const auto *gi0_85 = buffer.data(gi0 + 85);
    const auto *gi0_86 = buffer.data(gi0 + 86);
    const auto *gi0_87 = buffer.data(gi0 + 87);
    const auto *gi0_88 = buffer.data(gi0 + 88);
    const auto *gi0_89 = buffer.data(gi0 + 89);
    const auto *gi0_90 = buffer.data(gi0 + 90);
    const auto *gi0_91 = buffer.data(gi0 + 91);
    const auto *gi0_92 = buffer.data(gi0 + 92);
    const auto *gi0_93 = buffer.data(gi0 + 93);
    const auto *gi0_94 = buffer.data(gi0 + 94);
    const auto *gi0_96 = buffer.data(gi0 + 96);
    const auto *gi0_97 = buffer.data(gi0 + 97);
    const auto *gi0_98 = buffer.data(gi0 + 98);
    const auto *gi0_99 = buffer.data(gi0 + 99);
    const auto *gi0_100 = buffer.data(gi0 + 100);
    const auto *gi0_101 = buffer.data(gi0 + 101);
    const auto *gi0_102 = buffer.data(gi0 + 102);
    const auto *gi0_103 = buffer.data(gi0 + 103);
    const auto *gi0_104 = buffer.data(gi0 + 104);
    const auto *gi0_105 = buffer.data(gi0 + 105);
    const auto *gi0_106 = buffer.data(gi0 + 106);
    const auto *gi0_107 = buffer.data(gi0 + 107);
    const auto *gi0_108 = buffer.data(gi0 + 108);
    const auto *gi0_109 = buffer.data(gi0 + 109);
    const auto *gi0_110 = buffer.data(gi0 + 110);
    const auto *gi0_111 = buffer.data(gi0 + 111);
    const auto *gi0_112 = buffer.data(gi0 + 112);
    const auto *gi0_113 = buffer.data(gi0 + 113);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_1 = buffer.data(gi1 + 1);
    const auto *gi1_2 = buffer.data(gi1 + 2);
    const auto *gi1_3 = buffer.data(gi1 + 3);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_9 = buffer.data(gi1 + 9);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_15 = buffer.data(gi1 + 15);
    const auto *gi1_16 = buffer.data(gi1 + 16);
    const auto *gi1_17 = buffer.data(gi1 + 17);
    const auto *gi1_18 = buffer.data(gi1 + 18);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_33 = buffer.data(gi1 + 33);
    const auto *gi1_34 = buffer.data(gi1 + 34);
    const auto *gi1_35 = buffer.data(gi1 + 35);
    const auto *gi1_36 = buffer.data(gi1 + 36);
    const auto *gi1_37 = buffer.data(gi1 + 37);
    const auto *gi1_38 = buffer.data(gi1 + 38);
    const auto *gi1_39 = buffer.data(gi1 + 39);
    const auto *gi1_40 = buffer.data(gi1 + 40);
    const auto *gi1_41 = buffer.data(gi1 + 41);
    const auto *gi1_42 = buffer.data(gi1 + 42);
    const auto *gi1_43 = buffer.data(gi1 + 43);
    const auto *gi1_44 = buffer.data(gi1 + 44);
    const auto *gi1_45 = buffer.data(gi1 + 45);
    const auto *gi1_46 = buffer.data(gi1 + 46);
    const auto *gi1_47 = buffer.data(gi1 + 47);
    const auto *gi1_48 = buffer.data(gi1 + 48);
    const auto *gi1_49 = buffer.data(gi1 + 49);
    const auto *gi1_52 = buffer.data(gi1 + 52);
    const auto *gi1_53 = buffer.data(gi1 + 53);
    const auto *gi1_54 = buffer.data(gi1 + 54);
    const auto *gi1_55 = buffer.data(gi1 + 55);
    const auto *gi1_56 = buffer.data(gi1 + 56);
    const auto *gi1_57 = buffer.data(gi1 + 57);
    const auto *gi1_58 = buffer.data(gi1 + 58);
    const auto *gi1_59 = buffer.data(gi1 + 59);
    const auto *gi1_60 = buffer.data(gi1 + 60);
    const auto *gi1_61 = buffer.data(gi1 + 61);
    const auto *gi1_62 = buffer.data(gi1 + 62);
    const auto *gi1_63 = buffer.data(gi1 + 63);
    const auto *gi1_64 = buffer.data(gi1 + 64);
    const auto *gi1_65 = buffer.data(gi1 + 65);
    const auto *gi1_66 = buffer.data(gi1 + 66);
    const auto *gi1_67 = buffer.data(gi1 + 67);
    const auto *gi1_68 = buffer.data(gi1 + 68);
    const auto *gi1_69 = buffer.data(gi1 + 69);
    const auto *gi1_83 = buffer.data(gi1 + 83);
    const auto *gi1_85 = buffer.data(gi1 + 85);
    const auto *gi1_86 = buffer.data(gi1 + 86);
    const auto *gi1_87 = buffer.data(gi1 + 87);
    const auto *gi1_88 = buffer.data(gi1 + 88);
    const auto *gi1_89 = buffer.data(gi1 + 89);
    const auto *gi1_90 = buffer.data(gi1 + 90);
    const auto *gi1_91 = buffer.data(gi1 + 91);
    const auto *gi1_92 = buffer.data(gi1 + 92);
    const auto *gi1_93 = buffer.data(gi1 + 93);
    const auto *gi1_94 = buffer.data(gi1 + 94);
    const auto *gi1_95 = buffer.data(gi1 + 95);
    const auto *gi1_96 = buffer.data(gi1 + 96);
    const auto *gi1_97 = buffer.data(gi1 + 97);
    const auto *gi1_98 = buffer.data(gi1 + 98);
    const auto *gi1_99 = buffer.data(gi1 + 99);
    const auto *gi1_100 = buffer.data(gi1 + 100);
    const auto *gi1_101 = buffer.data(gi1 + 101);
    const auto *gi1_112 = buffer.data(gi1 + 112);
    const auto *gi1_113 = buffer.data(gi1 + 113);
    const auto *gi1_114 = buffer.data(gi1 + 114);
    const auto *gi1_115 = buffer.data(gi1 + 115);
    const auto *gi1_116 = buffer.data(gi1 + 116);
    const auto *gi1_117 = buffer.data(gi1 + 117);
    const auto *gi1_118 = buffer.data(gi1 + 118);
    const auto *gi1_119 = buffer.data(gi1 + 119);
    const auto *gi1_120 = buffer.data(gi1 + 120);
    const auto *gi1_121 = buffer.data(gi1 + 121);
    const auto *gi1_122 = buffer.data(gi1 + 122);
    const auto *gi1_123 = buffer.data(gi1 + 123);
    const auto *gi1_124 = buffer.data(gi1 + 124);
    const auto *gi1_125 = buffer.data(gi1 + 125);
    const auto *gi1_126 = buffer.data(gi1 + 126);
    const auto *gi1_127 = buffer.data(gi1 + 127);
    const auto *gi1_128 = buffer.data(gi1 + 128);
    const auto *gi1_129 = buffer.data(gi1 + 129);
    const auto *gi1_140 = buffer.data(gi1 + 140);
    const auto *gi1_142 = buffer.data(gi1 + 142);
    const auto *gi1_143 = buffer.data(gi1 + 143);
    const auto *gi1_144 = buffer.data(gi1 + 144);
    const auto *gi1_145 = buffer.data(gi1 + 145);
    const auto *gi1_146 = buffer.data(gi1 + 146);
    const auto *gi1_147 = buffer.data(gi1 + 147);
    const auto *gi1_148 = buffer.data(gi1 + 148);
    const auto *gi1_149 = buffer.data(gi1 + 149);
    const auto *gi1_150 = buffer.data(gi1 + 150);
    const auto *gi1_151 = buffer.data(gi1 + 151);
    const auto *gi1_152 = buffer.data(gi1 + 152);
    const auto *gi1_153 = buffer.data(gi1 + 153);
    const auto *gi1_154 = buffer.data(gi1 + 154);
    const auto *gi1_155 = buffer.data(gi1 + 155);
    const auto *gi1_156 = buffer.data(gi1 + 156);
    const auto *gi1_157 = buffer.data(gi1 + 157);
    const auto *gi1_158 = buffer.data(gi1 + 158);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_1 = buffer.data(gk + 1);
    const auto *gk_2 = buffer.data(gk + 2);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_16 = buffer.data(gk + 16);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_47 = buffer.data(gk + 47);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_49 = buffer.data(gk + 49);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_53 = buffer.data(gk + 53);
    const auto *gk_54 = buffer.data(gk + 54);
    const auto *gk_55 = buffer.data(gk + 55);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_58 = buffer.data(gk + 58);
    const auto *gk_59 = buffer.data(gk + 59);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_65 = buffer.data(gk + 65);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_73 = buffer.data(gk + 73);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_83 = buffer.data(gk + 83);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_85 = buffer.data(gk + 85);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_88 = buffer.data(gk + 88);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_95 = buffer.data(gk + 95);
    const auto *gk_96 = buffer.data(gk + 96);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_106 = buffer.data(gk + 106);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_110 = buffer.data(gk + 110);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_112 = buffer.data(gk + 112);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_115 = buffer.data(gk + 115);
    const auto *gk_116 = buffer.data(gk + 116);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_119 = buffer.data(gk + 119);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_121 = buffer.data(gk + 121);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_124 = buffer.data(gk + 124);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_127 = buffer.data(gk + 127);
    const auto *gk_132 = buffer.data(gk + 132);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_148 = buffer.data(gk + 148);
    const auto *gk_149 = buffer.data(gk + 149);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_151 = buffer.data(gk + 151);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_160 = buffer.data(gk + 160);
    const auto *gk_161 = buffer.data(gk + 161);
    const auto *gk_162 = buffer.data(gk + 162);
    const auto *gk_163 = buffer.data(gk + 163);
    const auto *gk_164 = buffer.data(gk + 164);
    const auto *gk_169 = buffer.data(gk + 169);
    const auto *gk_176 = buffer.data(gk + 176);
    const auto *gk_177 = buffer.data(gk + 177);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_181 = buffer.data(gk + 181);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);
    const auto *gk_184 = buffer.data(gk + 184);
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_187 = buffer.data(gk + 187);
    const auto *gk_188 = buffer.data(gk + 188);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_191 = buffer.data(gk + 191);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_193 = buffer.data(gk + 193);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_199 = buffer.data(gk + 199);
    const auto *gk_200 = buffer.data(gk + 200);
    const auto *gk_201 = buffer.data(gk + 201);
    const auto *gk_202 = buffer.data(gk + 202);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fk_0, gi0_0, gi0_1, gi1_0, \
                         gi1_1, gk_0, gk_1, gk_2, gk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = f_3 * gi0_0[k]
                 - f_4 * gi1_0[k]
                 + pb_y[k] * gk_1[k];

        t_2[k] = f_3 * gi0_0[k]
                 - f_4 * gi1_0[k]
                 + pb_z[k] * gk_2[k];

        t_3[k] = f_5 * gi0_1[k]
                 - f_6 * gi1_1[k]
                 + pb_y[k] * gk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, gi0_2, gi0_3, gi0_4, gi1_2, gi1_3, \
                         gi1_4, gk_4, gk_5, gk_6, gk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gi0_2[k]
                 - f_6 * gi1_2[k]
                 + pb_z[k] * gk_4[k];

        t_5[k] = f_7 * gi0_3[k]
                 - f_8 * gi1_3[k]
                 + pb_y[k] * gk_5[k];

        t_6[k] = f_3 * gi0_4[k]
                 - f_4 * gi1_4[k]
                 + pb_y[k] * gk_6[k];

        t_7[k] = f_7 * gi0_4[k]
                 - f_8 * gi1_4[k]
                 + pb_z[k] * gk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, gi0_5, gi0_6, gi0_7, gi1_5, gi1_6, \
                         gi1_7, gk_8, gk_9, gk_10, gk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * gi0_5[k]
                 - f_10 * gi1_5[k]
                 + pb_y[k] * gk_8[k];

        t_9[k] = f_5 * gi0_6[k]
                 - f_6 * gi1_6[k]
                 + pb_y[k] * gk_9[k];

        t_10[k] = f_3 * gi0_7[k]
                  - f_4 * gi1_7[k]
                  + pb_y[k] * gk_10[k];

        t_11[k] = f_9 * gi0_7[k]
                  - f_10 * gi1_7[k]
                  + pb_z[k] * gk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, gi0_8, gi0_9, gi0_10, gi1_8, gi1_9, gi1_10, \
                         gk_12, gk_13, gk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * gi0_8[k]
                  - f_12 * gi1_8[k]
                  + pb_y[k] * gk_12[k];

        t_13[k] = f_7 * gi0_9[k]
                  - f_8 * gi1_9[k]
                  + pb_y[k] * gk_13[k];

        t_14[k] = f_5 * gi0_10[k]
                  - f_6 * gi1_10[k]
                  + pb_y[k] * gk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, fk_17, fk_23, gi0_11, \
                         gi1_11, gk_15, gk_16, gk_17, gk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * gi0_11[k]
                  - f_4 * gi1_11[k]
                  + pb_y[k] * gk_15[k];

        t_16[k] = f_11 * gi0_11[k]
                  - f_12 * gi1_11[k]
                  + pb_z[k] * gk_16[k];

        t_17[k] = f_0 * fk_17[k]
                  + pb_x[k] * gk_17[k];

        t_18[k] = f_0 * fk_23[k]
                  + pb_x[k] * gk_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, gi0_12, gi0_13, gi0_14, gi1_12, gi1_14, \
                         gi1_15, gk_17, gk_18, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * gi0_12[k]
                  - f_2 * gi1_12[k]
                  + pb_y[k] * gk_17[k];

        t_20[k] = f_11 * gi0_13[k]
                  - f_12 * gi1_14[k]
                  + pb_y[k] * gk_18[k];

        t_21[k] = f_9 * gi0_14[k]
                  - f_10 * gi1_15[k]
                  + pb_y[k] * gk_19[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_y, pb_z, gi0_15, gi0_16, gi0_17, gi1_16, \
                         gi1_17, gi1_18, gk_20, gk_21, gk_22, gk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * gi0_15[k]
                  - f_8 * gi1_16[k]
                  + pb_y[k] * gk_20[k];

        t_23[k] = f_5 * gi0_16[k]
                  - f_6 * gi1_17[k]
                  + pb_y[k] * gk_21[k];

        t_24[k] = f_3 * gi0_17[k]
                  - f_4 * gi1_18[k]
                  + pb_y[k] * gk_22[k];

        t_25[k] = f_1 * gi0_17[k]
                  - f_2 * gi1_18[k]
                  + pb_z[k] * gk_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pb_y, fk_0, fk_1, fk_3, fk_5, \
                         fl_0, fl_1, fl_3, fl_5, gk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * fl_0[k];

        t_27[k] = f_13 * fk_0[k]
                  + pb_y[k] * gk_24[k];

        t_28[k] = f_14 * fk_1[k]
                  + pa_y[k] * fl_1[k];

        t_29[k] = f_15 * fk_3[k]
                  + pa_y[k] * fl_3[k];

        t_30[k] = f_0 * fk_5[k]
                  + pa_y[k] * fl_5[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pb_x, fk_8, fk_12, fk_17, fk_29, fl_8, \
                         fl_12, fl_17, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_16 * fk_8[k]
                  + pa_y[k] * fl_8[k];

        t_32[k] = f_17 * fk_12[k]
                  + pa_y[k] * fl_12[k];

        t_33[k] = f_15 * fk_29[k]
                  + pb_x[k] * gk_29[k];

        t_34[k] = f_18 * fk_17[k]
                  + pa_y[k] * fl_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, pb_z, fk_0, fk_2, fk_4, fk_6, \
                         fl_0, fl_2, fl_4, fl_6, gk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * fl_0[k];

        t_36[k] = f_13 * fk_0[k]
                  + pb_z[k] * gk_30[k];

        t_37[k] = f_14 * fk_2[k]
                  + pa_z[k] * fl_2[k];

        t_38[k] = f_15 * fk_4[k]
                  + pa_z[k] * fl_4[k];

        t_39[k] = f_14 * fk_6[k]
                  + pa_z[k] * fl_6[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_z, fk_7, fk_9, fk_10, fk_11, fk_13, \
                         fl_7, fl_9, fl_10, fl_11, fl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * fk_7[k]
                  + pa_z[k] * fl_7[k];

        t_41[k] = f_14 * fk_9[k]
                  + pa_z[k] * fl_9[k];

        t_42[k] = f_15 * fk_10[k]
                  + pa_z[k] * fl_10[k];

        t_43[k] = f_16 * fk_11[k]
                  + pa_z[k] * fl_11[k];

        t_44[k] = f_14 * fk_13[k]
                  + pa_z[k] * fl_13[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_z, pb_x, fk_14, fk_15, fk_16, fk_40, \
                         fl_14, fl_15, fl_16, gk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * fk_14[k]
                  + pa_z[k] * fl_14[k];

        t_46[k] = f_0 * fk_15[k]
                  + pa_z[k] * fl_15[k];

        t_47[k] = f_17 * fk_16[k]
                  + pa_z[k] * fl_16[k];

        t_48[k] = f_15 * fk_40[k]
                  + pb_x[k] * gk_40[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_z, fk_18, fk_19, fk_20, fk_21, \
                         fk_22, fl_18, fl_19, fl_20, fl_21, fl_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_14 * fk_18[k]
                  + pa_z[k] * fl_18[k];

        t_50[k] = f_15 * fk_19[k]
                  + pa_z[k] * fl_19[k];

        t_51[k] = f_0 * fk_20[k]
                  + pa_z[k] * fl_20[k];

        t_52[k] = f_16 * fk_21[k]
                  + pa_z[k] * fl_21[k];

        t_53[k] = f_17 * fk_22[k]
                  + pa_z[k] * fl_22[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pa_z, pb_y, dl0_0, dl1_0, fk_23, fk_24, \
                         fl_23, fl_24, gk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_18 * fk_23[k]
                  + pa_z[k] * fl_23[k];

        t_55[k] = f_19 * dl0_0[k]
                  - f_20 * dl1_0[k]
                  + pa_y[k] * fl_24[k];

        t_56[k] = f_14 * fk_24[k]
                  + pb_y[k] * gk_41[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, pb_z, fk_42, fk_43, gi0_20, gi0_22, gi0_24, \
                         gi1_32, gi1_34, gi1_36, gk_42, gk_43, gk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_14 * fk_42[k]
                  + f_11 * gi0_22[k]
                  - f_12 * gi1_34[k]
                  + pb_x[k] * gk_43[k];

        t_58[k] = f_3 * gi0_20[k]
                  - f_4 * gi1_32[k]
                  + pb_z[k] * gk_42[k];

        t_59[k] = f_14 * fk_43[k]
                  + f_9 * gi0_24[k]
                  - f_10 * gi1_36[k]
                  + pb_x[k] * gk_45[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_x, pb_z, fk_44, gi0_21, gi0_22, gi0_27, gi1_33, \
                         gi1_34, gi1_39, gk_44, gk_46, gk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * gi0_21[k]
                  - f_6 * gi1_33[k]
                  + pb_z[k] * gk_44[k];

        t_61[k] = f_14 * fk_44[k]
                  + f_7 * gi0_27[k]
                  - f_8 * gi1_39[k]
                  + pb_x[k] * gk_48[k];

        t_62[k] = f_3 * gi0_22[k]
                  - f_4 * gi1_34[k]
                  + pb_z[k] * gk_46[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, pb_z, fk_45, gi0_23, gi0_24, gi0_31, gi1_35, \
                         gi1_36, gi1_43, gk_47, gk_49, gk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * gi0_23[k]
                  - f_8 * gi1_35[k]
                  + pb_z[k] * gk_47[k];

        t_64[k] = f_14 * fk_45[k]
                  + f_5 * gi0_31[k]
                  - f_6 * gi1_43[k]
                  + pb_x[k] * gk_52[k];

        t_65[k] = f_3 * gi0_24[k]
                  - f_4 * gi1_36[k]
                  + pb_z[k] * gk_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, fk_46, gi0_25, gi0_26, gi0_32, gi1_37, \
                         gi1_38, gi1_44, gk_50, gk_51, gk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * gi0_25[k]
                  - f_6 * gi1_37[k]
                  + pb_z[k] * gk_50[k];

        t_67[k] = f_9 * gi0_26[k]
                  - f_10 * gi1_38[k]
                  + pb_z[k] * gk_51[k];

        t_68[k] = f_14 * fk_46[k]
                  + f_3 * gi0_32[k]
                  - f_4 * gi1_44[k]
                  + pb_x[k] * gk_57[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, gi0_27, gi0_28, gi0_29, gi1_39, gi1_40, \
                         gi1_41, gk_53, gk_54, gk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * gi0_27[k]
                  - f_4 * gi1_39[k]
                  + pb_z[k] * gk_53[k];

        t_70[k] = f_5 * gi0_28[k]
                  - f_6 * gi1_40[k]
                  + pb_z[k] * gk_54[k];

        t_71[k] = f_7 * gi0_29[k]
                  - f_8 * gi1_41[k]
                  + pb_z[k] * gk_55[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_x, pb_x, pb_z, dl0_1, dl1_1, fk_47, fl_26, \
                         gi0_30, gi1_42, gk_56, gk_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * gi0_30[k]
                  - f_12 * gi1_42[k]
                  + pb_z[k] * gk_56[k];

        t_73[k] = f_14 * fk_47[k]
                  + pb_x[k] * gk_58[k];

        t_74[k] = f_19 * dl0_1[k]
                  - f_20 * dl1_1[k]
                  + pa_x[k] * fl_26[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, gi0_32, gi0_33, gi0_34, gi1_44, gi1_45, \
                         gi1_46, gk_59, gk_60, gk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * gi0_32[k]
                  - f_4 * gi1_44[k]
                  + pb_z[k] * gk_59[k];

        t_76[k] = f_5 * gi0_33[k]
                  - f_6 * gi1_45[k]
                  + pb_z[k] * gk_60[k];

        t_77[k] = f_7 * gi0_34[k]
                  - f_8 * gi1_46[k]
                  + pb_z[k] * gk_61[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_z, gi0_35, gi0_36, gi0_37, gi1_47, gi1_48, \
                         gi1_49, gk_62, gk_63, gk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * gi0_35[k]
                  - f_10 * gi1_47[k]
                  + pb_z[k] * gk_62[k];

        t_79[k] = f_11 * gi0_36[k]
                  - f_12 * gi1_48[k]
                  + pb_z[k] * gk_63[k];

        t_80[k] = f_1 * gi0_37[k]
                  - f_2 * gi1_49[k]
                  + pb_z[k] * gk_64[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, pb_z, dl0_0, dl1_0, fk_30, fl_25, \
                         gi0_38, gi1_52, gk_65, gk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_19 * dl0_0[k]
                  - f_20 * dl1_0[k]
                  + pa_z[k] * fl_25[k];

        t_82[k] = f_14 * fk_30[k]
                  + pb_z[k] * gk_65[k];

        t_83[k] = f_3 * gi0_38[k]
                  - f_4 * gi1_52[k]
                  + pb_y[k] * gk_66[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, fk_50, fk_51, gi0_39, gi0_41, gi0_44, \
                         gi1_53, gi1_55, gi1_58, gk_68, gk_69, gk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_14 * fk_50[k]
                  + f_11 * gi0_41[k]
                  - f_12 * gi1_55[k]
                  + pb_x[k] * gk_69[k];

        t_85[k] = f_5 * gi0_39[k]
                  - f_6 * gi1_53[k]
                  + pb_y[k] * gk_68[k];

        t_86[k] = f_14 * fk_51[k]
                  + f_9 * gi0_44[k]
                  - f_10 * gi1_58[k]
                  + pb_x[k] * gk_72[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, fk_52, gi0_40, gi0_41, gi0_48, gi1_54, \
                         gi1_55, gi1_62, gk_70, gk_71, gk_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_7 * gi0_40[k]
                  - f_8 * gi1_54[k]
                  + pb_y[k] * gk_70[k];

        t_88[k] = f_3 * gi0_41[k]
                  - f_4 * gi1_55[k]
                  + pb_y[k] * gk_71[k];

        t_89[k] = f_14 * fk_52[k]
                  + f_7 * gi0_48[k]
                  - f_8 * gi1_62[k]
                  + pb_x[k] * gk_76[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_y, gi0_42, gi0_43, gi0_44, gi1_56, gi1_57, \
                         gi1_58, gk_73, gk_74, gk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * gi0_42[k]
                  - f_10 * gi1_56[k]
                  + pb_y[k] * gk_73[k];

        t_91[k] = f_5 * gi0_43[k]
                  - f_6 * gi1_57[k]
                  + pb_y[k] * gk_74[k];

        t_92[k] = f_3 * gi0_44[k]
                  - f_4 * gi1_58[k]
                  + pb_y[k] * gk_75[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, pb_y, fk_53, gi0_45, gi0_46, gi0_49, gi1_59, \
                         gi1_60, gi1_63, gk_77, gk_78, gk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_14 * fk_53[k]
                  + f_5 * gi0_49[k]
                  - f_6 * gi1_63[k]
                  + pb_x[k] * gk_81[k];

        t_94[k] = f_11 * gi0_45[k]
                  - f_12 * gi1_59[k]
                  + pb_y[k] * gk_77[k];

        t_95[k] = f_7 * gi0_46[k]
                  - f_8 * gi1_60[k]
                  + pb_y[k] * gk_78[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_x, pb_y, fk_54, gi0_47, gi0_48, gi0_55, gi1_61, \
                         gi1_62, gi1_69, gk_79, gk_80, gk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * gi0_47[k]
                  - f_6 * gi1_61[k]
                  + pb_y[k] * gk_79[k];

        t_97[k] = f_3 * gi0_48[k]
                  - f_4 * gi1_62[k]
                  + pb_y[k] * gk_80[k];

        t_98[k] = f_14 * fk_54[k]
                  + f_3 * gi0_55[k]
                  - f_4 * gi1_69[k]
                  + pb_x[k] * gk_82[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_y, fk_55, gi0_50, gi0_51, gi1_64, \
                         gi1_65, gk_83, gk_84, gk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_14 * fk_55[k]
                  + pb_x[k] * gk_89[k];

        t_100[k] = f_1 * gi0_50[k]
                   - f_2 * gi1_64[k]
                   + pb_y[k] * gk_83[k];

        t_101[k] = f_11 * gi0_51[k]
                   - f_12 * gi1_65[k]
                   + pb_y[k] * gk_84[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_y, gi0_52, gi0_53, gi0_54, gi1_66, gi1_67, \
                         gi1_68, gk_85, gk_86, gk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * gi0_52[k]
                   - f_10 * gi1_66[k]
                   + pb_y[k] * gk_85[k];

        t_103[k] = f_7 * gi0_53[k]
                   - f_8 * gi1_67[k]
                   + pb_y[k] * gk_86[k];

        t_104[k] = f_5 * gi0_54[k]
                   - f_6 * gi1_68[k]
                   + pb_y[k] * gk_87[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_y, dl0_2, dl1_2, fk_41, fk_56, \
                         fl_27, fl_28, gi0_55, gi1_69, gk_88, gk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * gi0_55[k]
                   - f_4 * gi1_69[k]
                   + pb_y[k] * gk_88[k];

        t_106[k] = f_19 * dl0_2[k]
                   - f_20 * dl1_2[k]
                   + pa_x[k] * fl_27[k];

        t_107[k] = f_18 * fk_56[k]
                   + pa_x[k] * fl_28[k];

        t_108[k] = f_15 * fk_41[k]
                   + pb_y[k] * gk_90[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pa_x, fk_58, fk_60, fk_63, fk_67, \
                         fk_72, fl_29, fl_31, fl_33, fl_36, fl_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_17 * fk_58[k]
                   + pa_x[k] * fl_29[k];

        t_110[k] = f_16 * fk_60[k]
                   + pa_x[k] * fl_31[k];

        t_111[k] = f_0 * fk_63[k]
                   + pa_x[k] * fl_33[k];

        t_112[k] = f_15 * fk_67[k]
                   + pa_x[k] * fl_36[k];

        t_113[k] = f_14 * fk_72[k]
                   + pa_x[k] * fl_40[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_x, pb_x, pb_z, fk_48, fk_73, fk_105, \
                         fl_45, fl_54, gk_95, gk_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * fk_73[k]
                   + pb_x[k] * gk_95[k];

        t_115[k] = pa_x[k] * fl_45[k];

        t_116[k] = f_18 * fk_105[k]
                   + pa_x[k] * fl_54[k];

        t_117[k] = f_15 * fk_48[k]
                   + pb_z[k] * gk_96[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_x, fk_109, fk_112, fk_116, \
                         fk_121, fk_122, fl_56, fl_58, fl_61, fl_65, \
                         fl_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_17 * fk_109[k]
                   + pa_x[k] * fl_56[k];

        t_119[k] = f_16 * fk_112[k]
                   + pa_x[k] * fl_58[k];

        t_120[k] = f_0 * fk_116[k]
                   + pa_x[k] * fl_61[k];

        t_121[k] = f_15 * fk_121[k]
                   + pa_x[k] * fl_65[k];

        t_122[k] = f_14 * fk_122[k]
                   + pa_x[k] * fl_70[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pb_x, pb_y, fk_56, fk_130, fl_77, \
                         gi0_58, gi1_83, gk_102, gk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_13 * fk_130[k]
                   + pb_x[k] * gk_102[k];

        t_124[k] = pa_x[k] * fl_77[k];

        t_125[k] = f_1 * gi0_58[k]
                   - f_2 * gi1_83[k]
                   + pb_x[k] * gk_103[k];

        t_126[k] = f_0 * fk_56[k]
                   + pb_y[k] * gk_103[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pb_x, gi0_59, gi0_60, gi0_61, gi1_85, gi1_86, \
                         gi1_87, gk_104, gk_105, gk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_11 * gi0_59[k]
                   - f_12 * gi1_85[k]
                   + pb_x[k] * gk_104[k];

        t_128[k] = f_11 * gi0_60[k]
                   - f_12 * gi1_86[k]
                   + pb_x[k] * gk_105[k];

        t_129[k] = f_9 * gi0_61[k]
                   - f_10 * gi1_87[k]
                   + pb_x[k] * gk_106[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pb_x, gi0_62, gi0_63, gi0_64, gi1_88, gi1_89, \
                         gi1_90, gk_107, gk_108, gk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_9 * gi0_62[k]
                   - f_10 * gi1_88[k]
                   + pb_x[k] * gk_107[k];

        t_131[k] = f_7 * gi0_63[k]
                   - f_8 * gi1_89[k]
                   + pb_x[k] * gk_108[k];

        t_132[k] = f_7 * gi0_64[k]
                   - f_8 * gi1_90[k]
                   + pb_x[k] * gk_109[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_x, gi0_65, gi0_66, gi0_67, gi1_91, gi1_92, \
                         gi1_93, gk_110, gk_111, gk_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * gi0_65[k]
                   - f_8 * gi1_91[k]
                   + pb_x[k] * gk_110[k];

        t_134[k] = f_5 * gi0_66[k]
                   - f_6 * gi1_92[k]
                   + pb_x[k] * gk_111[k];

        t_135[k] = f_5 * gi0_67[k]
                   - f_6 * gi1_93[k]
                   + pb_x[k] * gk_112[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pb_x, gi0_68, gi0_69, gi0_70, gi1_94, gi1_95, \
                         gi1_96, gk_113, gk_114, gk_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_5 * gi0_68[k]
                   - f_6 * gi1_94[k]
                   + pb_x[k] * gk_113[k];

        t_137[k] = f_5 * gi0_69[k]
                   - f_6 * gi1_95[k]
                   + pb_x[k] * gk_114[k];

        t_138[k] = f_3 * gi0_70[k]
                   - f_4 * gi1_96[k]
                   + pb_x[k] * gk_115[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pb_x, gi0_72, gi0_73, gi0_74, gi1_98, gi1_99, \
                         gi1_100, gk_116, gk_117, gk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * gi0_72[k]
                   - f_4 * gi1_98[k]
                   + pb_x[k] * gk_116[k];

        t_140[k] = f_3 * gi0_73[k]
                   - f_4 * gi1_99[k]
                   + pb_x[k] * gk_117[k];

        t_141[k] = f_3 * gi0_74[k]
                   - f_4 * gi1_100[k]
                   + pb_x[k] * gk_118[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_x, pb_y, pb_z, fk_73, gi0_70, gi0_75, gi1_96, \
                         gi1_101, gk_119, gk_120, gk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * gi0_75[k]
                   - f_4 * gi1_101[k]
                   + pb_x[k] * gk_119[k];

        t_143[k] = f_0 * fk_73[k]
                   + f_1 * gi0_70[k]
                   - f_2 * gi1_96[k]
                   + pb_y[k] * gk_120[k];

        t_144[k] = f_3 * gi0_70[k]
                   - f_4 * gi1_96[k]
                   + pb_z[k] * gk_121[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_z, gi0_71, gi0_72, gi0_73, gi1_97, gi1_98, \
                         gi1_99, gk_122, gk_123, gk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_5 * gi0_71[k]
                   - f_6 * gi1_97[k]
                   + pb_z[k] * gk_122[k];

        t_146[k] = f_7 * gi0_72[k]
                   - f_8 * gi1_98[k]
                   + pb_z[k] * gk_123[k];

        t_147[k] = f_9 * gi0_73[k]
                   - f_10 * gi1_99[k]
                   + pb_z[k] * gk_124[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_z, pb_y, pb_z, fk_57, fk_80, fl_30, \
                         gi0_74, gi0_75, gi1_100, gi1_101, gk_125, \
                         gk_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_11 * gi0_74[k]
                   - f_12 * gi1_100[k]
                   + pb_z[k] * gk_125[k];

        t_149[k] = f_0 * fk_80[k]
                   + pb_y[k] * gk_127[k];

        t_150[k] = f_1 * gi0_75[k]
                   - f_2 * gi1_101[k]
                   + pb_z[k] * gk_127[k];

        t_151[k] = f_14 * fk_57[k]
                   + pa_z[k] * fl_30[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pa_z, fk_59, fk_61, fk_62, fk_64, \
                         fk_65, fl_32, fl_34, fl_35, fl_37, fl_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_15 * fk_59[k]
                   + pa_z[k] * fl_32[k];

        t_153[k] = f_14 * fk_61[k]
                   + pa_z[k] * fl_34[k];

        t_154[k] = f_0 * fk_62[k]
                   + pa_z[k] * fl_35[k];

        t_155[k] = f_14 * fk_64[k]
                   + pa_z[k] * fl_37[k];

        t_156[k] = f_15 * fk_65[k]
                   + pa_z[k] * fl_38[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_z, fk_66, fk_68, fk_69, fk_70, \
                         fk_71, fl_39, fl_41, fl_42, fl_43, fl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_16 * fk_66[k]
                   + pa_z[k] * fl_39[k];

        t_158[k] = f_14 * fk_68[k]
                   + pa_z[k] * fl_41[k];

        t_159[k] = f_15 * fk_69[k]
                   + pa_z[k] * fl_42[k];

        t_160[k] = f_0 * fk_70[k]
                   + pa_z[k] * fl_43[k];

        t_161[k] = f_17 * fk_71[k]
                   + pa_z[k] * fl_44[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pa_z, pb_z, fk_73, fk_74, fk_75, \
                         fk_76, fl_45, fl_46, fl_47, fl_48, gk_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_z[k] * fl_45[k];

        t_163[k] = f_13 * fk_73[k]
                   + pb_z[k] * gk_132[k];

        t_164[k] = f_14 * fk_74[k]
                   + pa_z[k] * fl_46[k];

        t_165[k] = f_15 * fk_75[k]
                   + pa_z[k] * fl_47[k];

        t_166[k] = f_0 * fk_76[k]
                   + pa_z[k] * fl_48[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_z, pb_y, fk_77, fk_78, fk_80, fk_92, \
                         fl_49, fl_50, fl_51, gk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_16 * fk_77[k]
                   + pa_z[k] * fl_49[k];

        t_168[k] = f_17 * fk_78[k]
                   + pa_z[k] * fl_50[k];

        t_169[k] = f_15 * fk_92[k]
                   + pb_y[k] * gk_139[k];

        t_170[k] = f_18 * fk_80[k]
                   + pa_z[k] * fl_51[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_x, gi0_77, gi0_78, gi0_79, gi1_112, gi1_113, \
                         gi1_114, gk_140, gk_141, gk_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_1 * gi0_77[k]
                   - f_2 * gi1_112[k]
                   + pb_x[k] * gk_140[k];

        t_172[k] = f_11 * gi0_78[k]
                   - f_12 * gi1_113[k]
                   + pb_x[k] * gk_141[k];

        t_173[k] = f_11 * gi0_79[k]
                   - f_12 * gi1_114[k]
                   + pb_x[k] * gk_142[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_x, gi0_80, gi0_81, gi0_82, gi1_115, gi1_116, \
                         gi1_117, gk_143, gk_144, gk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_9 * gi0_80[k]
                   - f_10 * gi1_115[k]
                   + pb_x[k] * gk_143[k];

        t_175[k] = f_9 * gi0_81[k]
                   - f_10 * gi1_116[k]
                   + pb_x[k] * gk_144[k];

        t_176[k] = f_7 * gi0_82[k]
                   - f_8 * gi1_117[k]
                   + pb_x[k] * gk_145[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_x, gi0_83, gi0_84, gi0_85, gi1_118, gi1_119, \
                         gi1_120, gk_146, gk_147, gk_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_7 * gi0_83[k]
                   - f_8 * gi1_118[k]
                   + pb_x[k] * gk_146[k];

        t_178[k] = f_7 * gi0_84[k]
                   - f_8 * gi1_119[k]
                   + pb_x[k] * gk_147[k];

        t_179[k] = f_5 * gi0_85[k]
                   - f_6 * gi1_120[k]
                   + pb_x[k] * gk_148[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pb_x, gi0_86, gi0_87, gi0_88, gi1_121, gi1_122, \
                         gi1_123, gk_149, gk_150, gk_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_5 * gi0_86[k]
                   - f_6 * gi1_121[k]
                   + pb_x[k] * gk_149[k];

        t_181[k] = f_5 * gi0_87[k]
                   - f_6 * gi1_122[k]
                   + pb_x[k] * gk_150[k];

        t_182[k] = f_5 * gi0_88[k]
                   - f_6 * gi1_123[k]
                   + pb_x[k] * gk_151[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_x, gi0_89, gi0_90, gi0_91, gi1_124, gi1_125, \
                         gi1_126, gk_152, gk_153, gk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_3 * gi0_89[k]
                   - f_4 * gi1_124[k]
                   + pb_x[k] * gk_152[k];

        t_184[k] = f_3 * gi0_90[k]
                   - f_4 * gi1_125[k]
                   + pb_x[k] * gk_153[k];

        t_185[k] = f_3 * gi0_91[k]
                   - f_4 * gi1_126[k]
                   + pb_x[k] * gk_154[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_z, pb_x, dl0_1, dl1_1, fl_52, gi0_92, gi0_94, \
                         gi1_127, gi1_129, gk_155, gk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_3 * gi0_92[k]
                   - f_4 * gi1_127[k]
                   + pb_x[k] * gk_155[k];

        t_187[k] = f_3 * gi0_94[k]
                   - f_4 * gi1_129[k]
                   + pb_x[k] * gk_156[k];

        t_188[k] = f_19 * dl0_1[k]
                   - f_20 * dl1_1[k]
                   + pa_z[k] * fl_52[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, fk_85, fk_99, fk_100, gi0_90, \
                         gi0_91, gi1_125, gi1_126, gk_157, gk_159, \
                         gk_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * fk_85[k]
                   + pb_z[k] * gk_157[k];

        t_190[k] = f_14 * fk_99[k]
                   + f_11 * gi0_90[k]
                   - f_12 * gi1_125[k]
                   + pb_y[k] * gk_159[k];

        t_191[k] = f_14 * fk_100[k]
                   + f_9 * gi0_91[k]
                   - f_10 * gi1_126[k]
                   + pb_y[k] * gk_160[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, fk_101, fk_102, fk_103, gi0_92, gi0_93, \
                         gi0_94, gi1_127, gi1_128, gi1_129, gk_161, gk_162, \
                         gk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_14 * fk_101[k]
                   + f_7 * gi0_92[k]
                   - f_8 * gi1_127[k]
                   + pb_y[k] * gk_161[k];

        t_193[k] = f_14 * fk_102[k]
                   + f_5 * gi0_93[k]
                   - f_6 * gi1_128[k]
                   + pb_y[k] * gk_162[k];

        t_194[k] = f_14 * fk_103[k]
                   + f_3 * gi0_94[k]
                   - f_4 * gi1_129[k]
                   + pb_y[k] * gk_163[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pb_y, dl0_2, dl1_2, fk_104, fk_106, \
                         fk_108, fl_53, fl_55, fl_57, gk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_14 * fk_104[k]
                   + pb_y[k] * gk_164[k];

        t_196[k] = f_19 * dl0_2[k]
                   - f_20 * dl1_2[k]
                   + pa_y[k] * fl_53[k];

        t_197[k] = f_14 * fk_106[k]
                   + pa_y[k] * fl_55[k];

        t_198[k] = f_15 * fk_108[k]
                   + pa_y[k] * fl_57[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pa_y, fk_110, fk_111, fk_113, \
                         fk_114, fk_115, fl_59, fl_60, fl_62, fl_63, \
                         fl_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * fk_110[k]
                   + pa_y[k] * fl_59[k];

        t_200[k] = f_14 * fk_111[k]
                   + pa_y[k] * fl_60[k];

        t_201[k] = f_16 * fk_113[k]
                   + pa_y[k] * fl_62[k];

        t_202[k] = f_15 * fk_114[k]
                   + pa_y[k] * fl_63[k];

        t_203[k] = f_14 * fk_115[k]
                   + pa_y[k] * fl_64[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pa_y, fk_117, fk_118, fk_119, \
                         fk_120, fk_123, fl_66, fl_67, fl_68, fl_69, \
                         fl_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_17 * fk_117[k]
                   + pa_y[k] * fl_66[k];

        t_205[k] = f_0 * fk_118[k]
                   + pa_y[k] * fl_67[k];

        t_206[k] = f_15 * fk_119[k]
                   + pa_y[k] * fl_68[k];

        t_207[k] = f_14 * fk_120[k]
                   + pa_y[k] * fl_69[k];

        t_208[k] = f_18 * fk_123[k]
                   + pa_y[k] * fl_71[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pa_y, pb_z, fk_97, fk_125, fk_126, \
                         fk_127, fl_72, fl_73, fl_74, gk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_15 * fk_97[k]
                   + pb_z[k] * gk_169[k];

        t_210[k] = f_17 * fk_125[k]
                   + pa_y[k] * fl_72[k];

        t_211[k] = f_16 * fk_126[k]
                   + pa_y[k] * fl_73[k];

        t_212[k] = f_0 * fk_127[k]
                   + pa_y[k] * fl_74[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pb_y, fk_128, fk_129, fk_130, \
                         fl_75, fl_76, fl_77, gk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_15 * fk_128[k]
                   + pa_y[k] * fl_75[k];

        t_214[k] = f_14 * fk_129[k]
                   + pa_y[k] * fl_76[k];

        t_215[k] = f_13 * fk_130[k]
                   + pb_y[k] * gk_176[k];

        t_216[k] = pa_y[k] * fl_77[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_x, pb_z, fk_105, gi0_96, gi0_97, \
                         gi0_98, gi1_140, gi1_142, gi1_143, gk_177, gk_179, \
                         gk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_1 * gi0_96[k]
                   - f_2 * gi1_140[k]
                   + pb_x[k] * gk_177[k];

        t_218[k] = f_0 * fk_105[k]
                   + pb_z[k] * gk_177[k];

        t_219[k] = f_11 * gi0_97[k]
                   - f_12 * gi1_142[k]
                   + pb_x[k] * gk_179[k];

        t_220[k] = f_11 * gi0_98[k]
                   - f_12 * gi1_143[k]
                   + pb_x[k] * gk_180[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pb_x, gi0_99, gi0_100, gi0_101, gi1_144, \
                         gi1_145, gi1_146, gk_181, gk_182, gk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * gi0_99[k]
                   - f_10 * gi1_144[k]
                   + pb_x[k] * gk_181[k];

        t_222[k] = f_9 * gi0_100[k]
                   - f_10 * gi1_145[k]
                   + pb_x[k] * gk_182[k];

        t_223[k] = f_7 * gi0_101[k]
                   - f_8 * gi1_146[k]
                   + pb_x[k] * gk_183[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, pb_x, gi0_102, gi0_103, gi0_104, gi1_147, \
                         gi1_148, gi1_149, gk_184, gk_185, gk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_7 * gi0_102[k]
                   - f_8 * gi1_147[k]
                   + pb_x[k] * gk_184[k];

        t_225[k] = f_7 * gi0_103[k]
                   - f_8 * gi1_148[k]
                   + pb_x[k] * gk_185[k];

        t_226[k] = f_5 * gi0_104[k]
                   - f_6 * gi1_149[k]
                   + pb_x[k] * gk_186[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pb_x, gi0_105, gi0_106, gi0_107, gi1_150, \
                         gi1_151, gi1_152, gk_187, gk_188, gk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_5 * gi0_105[k]
                   - f_6 * gi1_150[k]
                   + pb_x[k] * gk_187[k];

        t_228[k] = f_5 * gi0_106[k]
                   - f_6 * gi1_151[k]
                   + pb_x[k] * gk_188[k];

        t_229[k] = f_5 * gi0_107[k]
                   - f_6 * gi1_152[k]
                   + pb_x[k] * gk_189[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pb_x, gi0_108, gi0_109, gi0_110, gi1_153, \
                         gi1_154, gi1_155, gk_190, gk_191, gk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_3 * gi0_108[k]
                   - f_4 * gi1_153[k]
                   + pb_x[k] * gk_190[k];

        t_231[k] = f_3 * gi0_109[k]
                   - f_4 * gi1_154[k]
                   + pb_x[k] * gk_191[k];

        t_232[k] = f_3 * gi0_110[k]
                   - f_4 * gi1_155[k]
                   + pb_x[k] * gk_192[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_x, pb_y, gi0_108, gi0_111, gi0_113, gi1_153, \
                         gi1_156, gi1_158, gk_193, gk_194, gk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * gi0_111[k]
                   - f_4 * gi1_156[k]
                   + pb_x[k] * gk_193[k];

        t_234[k] = f_3 * gi0_113[k]
                   - f_4 * gi1_158[k]
                   + pb_x[k] * gk_194[k];

        t_235[k] = f_1 * gi0_108[k]
                   - f_2 * gi1_153[k]
                   + pb_y[k] * gk_195[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_y, pb_z, fk_123, gi0_109, gi0_110, gi1_154, \
                         gi1_155, gk_195, gk_197, gk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * fk_123[k]
                   + pb_z[k] * gk_195[k];

        t_237[k] = f_11 * gi0_109[k]
                   - f_12 * gi1_154[k]
                   + pb_y[k] * gk_197[k];

        t_238[k] = f_9 * gi0_110[k]
                   - f_10 * gi1_155[k]
                   + pb_y[k] * gk_198[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pb_y, gi0_111, gi0_112, gi0_113, gi1_156, \
                         gi1_157, gi1_158, gk_199, gk_200, gk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_7 * gi0_111[k]
                   - f_8 * gi1_156[k]
                   + pb_y[k] * gk_199[k];

        t_240[k] = f_5 * gi0_112[k]
                   - f_6 * gi1_157[k]
                   + pb_y[k] * gk_200[k];

        t_241[k] = f_3 * gi0_113[k]
                   - f_4 * gi1_158[k]
                   + pb_y[k] * gk_201[k];
    }

#pragma omp simd aligned(t_242, pb_z, fk_130, gi0_113, gi1_158, \
                         gk_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * fk_130[k]
                   + f_1 * gi0_113[k]
                   - f_2 * gi1_158[k]
                   + pb_z[k] * gk_202[k];
    }
}

auto
compute_prim_gl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_23 = buffer.data(fk + 23);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_26 = buffer.data(gi0 + 26);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_17 = buffer.data(gi1 + 17);
    const auto *gi1_20 = buffer.data(gi1 + 20);
    const auto *gi1_21 = buffer.data(gi1 + 21);
    const auto *gi1_22 = buffer.data(gi1 + 22);
    const auto *gi1_23 = buffer.data(gi1 + 23);
    const auto *gi1_24 = buffer.data(gi1 + 24);
    const auto *gi1_26 = buffer.data(gi1 + 26);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_29 = buffer.data(gk + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_4, \
                         gi1_5, gi1_6, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_4[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_5[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_6[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_7, gi1_8, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_7[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_8[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_10, gi1_11, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_10[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_11[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_12, gi1_13, gi1_14, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_12[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_13[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_14[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_15, fl_4, fl_5, \
                         fl_8, gi0_17, gi1_17, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_15[k]
                  + f_1 * gi0_17[k]
                  - f_2 * gi1_17[k]
                  + pb_y[k] * gk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_17, fl_5, fl_6, \
                         gi0_20, gi1_20, gk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_17[k]
                  + f_6 * gi0_20[k]
                  - f_7 * gi1_20[k]
                  + pb_y[k] * gk_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_18, fk_19, fk_20, gi0_21, gi0_22, gi0_23, \
                         gi1_21, gi1_22, gi1_23, gk_23, gk_24, gk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_18[k]
                  + f_8 * gi0_21[k]
                  - f_9 * gi1_21[k]
                  + pb_y[k] * gk_23[k];

        t_24[k] = f_5 * fk_19[k]
                  + f_10 * gi0_22[k]
                  - f_11 * gi1_22[k]
                  + pb_y[k] * gk_24[k];

        t_25[k] = f_5 * fk_20[k]
                  + f_12 * gi0_23[k]
                  - f_13 * gi1_23[k]
                  + pb_y[k] * gk_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_21, fl_7, fl_8, \
                         gi0_24, gi1_24, gk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_21[k]
                  + f_14 * gi0_24[k]
                  - f_15 * gi1_24[k]
                  + pb_y[k] * gk_26[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_23, gi0_26, gi1_26, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_23[k]
                  + f_1 * gi0_26[k]
                  - f_2 * gi1_26[k]
                  + pb_z[k] * gk_29[k];
    }
}

auto
compute_prim_gl_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_23 = buffer.data(fk + 23);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_26 = buffer.data(gi0 + 26);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_33 = buffer.data(gi1 + 33);
    const auto *gi1_35 = buffer.data(gi1 + 35);
    const auto *gi1_37 = buffer.data(gi1 + 37);
    const auto *gi1_39 = buffer.data(gi1 + 39);
    const auto *gi1_40 = buffer.data(gi1 + 40);
    const auto *gi1_51 = buffer.data(gi1 + 51);
    const auto *gi1_53 = buffer.data(gi1 + 53);
    const auto *gi1_55 = buffer.data(gi1 + 55);
    const auto *gi1_56 = buffer.data(gi1 + 56);
    const auto *gi1_62 = buffer.data(gi1 + 62);
    const auto *gi1_91 = buffer.data(gi1 + 91);
    const auto *gi1_117 = buffer.data(gi1 + 117);
    const auto *gi1_118 = buffer.data(gi1 + 118);
    const auto *gi1_119 = buffer.data(gi1 + 119);
    const auto *gi1_120 = buffer.data(gi1 + 120);
    const auto *gi1_121 = buffer.data(gi1 + 121);
    const auto *gi1_152 = buffer.data(gi1 + 152);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_106 = buffer.data(gk + 106);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_183 = buffer.data(gk + 183);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_33, \
                         gi1_35, gi1_37, gk_42, gk_44, gk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_33[k]
                 + pb_x[k] * gk_42[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_35[k]
                 + pb_x[k] * gk_44[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_37[k]
                 + pb_x[k] * gk_46[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_39, gi1_40, gk_48, gk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_39[k]
                 + pb_x[k] * gk_48[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_40[k]
                 + pb_x[k] * gk_50[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_51, gi1_53, gk_60, gk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_51[k]
                  + pb_x[k] * gk_60[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_53[k]
                  + pb_x[k] * gk_62[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_55, gi1_56, gi1_62, gk_64, gk_66, gk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_55[k]
                  + pb_x[k] * gk_64[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_56[k]
                  + pb_x[k] * gk_66[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_62[k]
                  + pb_x[k] * gk_67[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_15, fl_4, fl_5, \
                         fl_8, gi0_17, gi1_91, gk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_15[k]
                  + f_1 * gi0_17[k]
                  - f_2 * gi1_91[k]
                  + pb_y[k] * gk_106[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_17, fl_5, fl_6, \
                         gi0_20, gi1_117, gk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_17[k]
                  + f_6 * gi0_20[k]
                  - f_7 * gi1_117[k]
                  + pb_y[k] * gk_140[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_18, fk_19, fk_20, gi0_21, gi0_22, gi0_23, \
                         gi1_118, gi1_119, gi1_120, gk_141, gk_142, \
                         gk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_18[k]
                  + f_8 * gi0_21[k]
                  - f_9 * gi1_118[k]
                  + pb_y[k] * gk_141[k];

        t_24[k] = f_5 * fk_19[k]
                  + f_10 * gi0_22[k]
                  - f_11 * gi1_119[k]
                  + pb_y[k] * gk_142[k];

        t_25[k] = f_5 * fk_20[k]
                  + f_12 * gi0_23[k]
                  - f_13 * gi1_120[k]
                  + pb_y[k] * gk_143[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_21, fl_7, fl_8, \
                         gi0_24, gi1_121, gk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_21[k]
                  + f_14 * gi0_24[k]
                  - f_15 * gi1_121[k]
                  + pb_y[k] * gk_144[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_23, gi0_26, gi1_152, gk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_23[k]
                  + f_1 * gi0_26[k]
                  - f_2 * gi1_152[k]
                  + pb_z[k] * gk_183[k];
    }
}

auto
compute_prim_gl_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
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
    const auto f_13 = 0.5 / alpha;
    const auto f_14 = 0.5 * beta / (alpha * p);
    const auto f_15 = 1.0 / p;

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
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
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_65 = buffer.data(fk + 65);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_1 = buffer.data(gi0 + 1);
    const auto *gi0_2 = buffer.data(gi0 + 2);
    const auto *gi0_3 = buffer.data(gi0 + 3);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_9 = buffer.data(gi0 + 9);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_16 = buffer.data(gi0 + 16);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_18 = buffer.data(gi0 + 18);
    const auto *gi0_19 = buffer.data(gi0 + 19);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_33 = buffer.data(gi0 + 33);
    const auto *gi0_35 = buffer.data(gi0 + 35);
    const auto *gi0_37 = buffer.data(gi0 + 37);
    const auto *gi0_39 = buffer.data(gi0 + 39);
    const auto *gi0_40 = buffer.data(gi0 + 40);
    const auto *gi0_51 = buffer.data(gi0 + 51);
    const auto *gi0_53 = buffer.data(gi0 + 53);
    const auto *gi0_55 = buffer.data(gi0 + 55);
    const auto *gi0_56 = buffer.data(gi0 + 56);
    const auto *gi0_62 = buffer.data(gi0 + 62);
    const auto *gi0_76 = buffer.data(gi0 + 76);
    const auto *gi0_78 = buffer.data(gi0 + 78);
    const auto *gi0_79 = buffer.data(gi0 + 79);
    const auto *gi0_80 = buffer.data(gi0 + 80);
    const auto *gi0_82 = buffer.data(gi0 + 82);
    const auto *gi0_83 = buffer.data(gi0 + 83);
    const auto *gi0_85 = buffer.data(gi0 + 85);
    const auto *gi0_86 = buffer.data(gi0 + 86);
    const auto *gi0_87 = buffer.data(gi0 + 87);
    const auto *gi0_88 = buffer.data(gi0 + 88);
    const auto *gi0_89 = buffer.data(gi0 + 89);
    const auto *gi0_90 = buffer.data(gi0 + 90);
    const auto *gi0_91 = buffer.data(gi0 + 91);
    const auto *gi0_92 = buffer.data(gi0 + 92);
    const auto *gi0_93 = buffer.data(gi0 + 93);
    const auto *gi0_94 = buffer.data(gi0 + 94);
    const auto *gi0_95 = buffer.data(gi0 + 95);
    const auto *gi0_96 = buffer.data(gi0 + 96);
    const auto *gi0_117 = buffer.data(gi0 + 117);
    const auto *gi0_118 = buffer.data(gi0 + 118);
    const auto *gi0_119 = buffer.data(gi0 + 119);
    const auto *gi0_120 = buffer.data(gi0 + 120);
    const auto *gi0_121 = buffer.data(gi0 + 121);
    const auto *gi0_132 = buffer.data(gi0 + 132);
    const auto *gi0_134 = buffer.data(gi0 + 134);
    const auto *gi0_135 = buffer.data(gi0 + 135);
    const auto *gi0_136 = buffer.data(gi0 + 136);
    const auto *gi0_138 = buffer.data(gi0 + 138);
    const auto *gi0_139 = buffer.data(gi0 + 139);
    const auto *gi0_140 = buffer.data(gi0 + 140);
    const auto *gi0_142 = buffer.data(gi0 + 142);
    const auto *gi0_143 = buffer.data(gi0 + 143);
    const auto *gi0_144 = buffer.data(gi0 + 144);
    const auto *gi0_145 = buffer.data(gi0 + 145);
    const auto *gi0_146 = buffer.data(gi0 + 146);
    const auto *gi0_147 = buffer.data(gi0 + 147);
    const auto *gi0_148 = buffer.data(gi0 + 148);
    const auto *gi0_149 = buffer.data(gi0 + 149);
    const auto *gi0_150 = buffer.data(gi0 + 150);
    const auto *gi0_151 = buffer.data(gi0 + 151);
    const auto *gi0_152 = buffer.data(gi0 + 152);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_1 = buffer.data(gi1 + 1);
    const auto *gi1_2 = buffer.data(gi1 + 2);
    const auto *gi1_3 = buffer.data(gi1 + 3);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_9 = buffer.data(gi1 + 9);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_15 = buffer.data(gi1 + 15);
    const auto *gi1_16 = buffer.data(gi1 + 16);
    const auto *gi1_17 = buffer.data(gi1 + 17);
    const auto *gi1_18 = buffer.data(gi1 + 18);
    const auto *gi1_25 = buffer.data(gi1 + 25);
    const auto *gi1_27 = buffer.data(gi1 + 27);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_31 = buffer.data(gi1 + 31);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_41 = buffer.data(gi1 + 41);
    const auto *gi1_43 = buffer.data(gi1 + 43);
    const auto *gi1_45 = buffer.data(gi1 + 45);
    const auto *gi1_46 = buffer.data(gi1 + 46);
    const auto *gi1_52 = buffer.data(gi1 + 52);
    const auto *gi1_65 = buffer.data(gi1 + 65);
    const auto *gi1_67 = buffer.data(gi1 + 67);
    const auto *gi1_68 = buffer.data(gi1 + 68);
    const auto *gi1_69 = buffer.data(gi1 + 69);
    const auto *gi1_70 = buffer.data(gi1 + 70);
    const auto *gi1_71 = buffer.data(gi1 + 71);
    const auto *gi1_72 = buffer.data(gi1 + 72);
    const auto *gi1_73 = buffer.data(gi1 + 73);
    const auto *gi1_74 = buffer.data(gi1 + 74);
    const auto *gi1_75 = buffer.data(gi1 + 75);
    const auto *gi1_76 = buffer.data(gi1 + 76);
    const auto *gi1_77 = buffer.data(gi1 + 77);
    const auto *gi1_78 = buffer.data(gi1 + 78);
    const auto *gi1_79 = buffer.data(gi1 + 79);
    const auto *gi1_80 = buffer.data(gi1 + 80);
    const auto *gi1_81 = buffer.data(gi1 + 81);
    const auto *gi1_82 = buffer.data(gi1 + 82);
    const auto *gi1_83 = buffer.data(gi1 + 83);
    const auto *gi1_99 = buffer.data(gi1 + 99);
    const auto *gi1_100 = buffer.data(gi1 + 100);
    const auto *gi1_101 = buffer.data(gi1 + 101);
    const auto *gi1_102 = buffer.data(gi1 + 102);
    const auto *gi1_103 = buffer.data(gi1 + 103);
    const auto *gi1_110 = buffer.data(gi1 + 110);
    const auto *gi1_112 = buffer.data(gi1 + 112);
    const auto *gi1_113 = buffer.data(gi1 + 113);
    const auto *gi1_114 = buffer.data(gi1 + 114);
    const auto *gi1_115 = buffer.data(gi1 + 115);
    const auto *gi1_116 = buffer.data(gi1 + 116);
    const auto *gi1_117 = buffer.data(gi1 + 117);
    const auto *gi1_118 = buffer.data(gi1 + 118);
    const auto *gi1_119 = buffer.data(gi1 + 119);
    const auto *gi1_120 = buffer.data(gi1 + 120);
    const auto *gi1_121 = buffer.data(gi1 + 121);
    const auto *gi1_122 = buffer.data(gi1 + 122);
    const auto *gi1_123 = buffer.data(gi1 + 123);
    const auto *gi1_124 = buffer.data(gi1 + 124);
    const auto *gi1_125 = buffer.data(gi1 + 125);
    const auto *gi1_126 = buffer.data(gi1 + 126);
    const auto *gi1_127 = buffer.data(gi1 + 127);
    const auto *gi1_128 = buffer.data(gi1 + 128);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_1 = buffer.data(gk + 1);
    const auto *gk_2 = buffer.data(gk + 2);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_16 = buffer.data(gk + 16);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_47 = buffer.data(gk + 47);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_49 = buffer.data(gk + 49);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_53 = buffer.data(gk + 53);
    const auto *gk_54 = buffer.data(gk + 54);
    const auto *gk_55 = buffer.data(gk + 55);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_58 = buffer.data(gk + 58);
    const auto *gk_59 = buffer.data(gk + 59);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_65 = buffer.data(gk + 65);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_83 = buffer.data(gk + 83);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_85 = buffer.data(gk + 85);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_88 = buffer.data(gk + 88);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_91 = buffer.data(gk + 91);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_93 = buffer.data(gk + 93);
    const auto *gk_94 = buffer.data(gk + 94);
    const auto *gk_95 = buffer.data(gk + 95);
    const auto *gk_96 = buffer.data(gk + 96);
    const auto *gk_97 = buffer.data(gk + 97);
    const auto *gk_98 = buffer.data(gk + 98);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fk_0, gi0_0, gi0_1, gi1_0, \
                         gi1_1, gk_0, gk_1, gk_2, gk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = f_3 * gi0_0[k]
                 - f_4 * gi1_0[k]
                 + pb_y[k] * gk_1[k];

        t_2[k] = f_3 * gi0_0[k]
                 - f_4 * gi1_0[k]
                 + pb_z[k] * gk_2[k];

        t_3[k] = f_5 * gi0_1[k]
                 - f_6 * gi1_1[k]
                 + pb_y[k] * gk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, gi0_2, gi0_3, gi0_4, gi1_2, gi1_3, \
                         gi1_4, gk_4, gk_5, gk_6, gk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gi0_2[k]
                 - f_6 * gi1_2[k]
                 + pb_z[k] * gk_4[k];

        t_5[k] = f_7 * gi0_3[k]
                 - f_8 * gi1_3[k]
                 + pb_y[k] * gk_5[k];

        t_6[k] = f_3 * gi0_4[k]
                 - f_4 * gi1_4[k]
                 + pb_y[k] * gk_6[k];

        t_7[k] = f_7 * gi0_4[k]
                 - f_8 * gi1_4[k]
                 + pb_z[k] * gk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, gi0_5, gi0_7, gi0_8, gi1_5, gi1_6, \
                         gi1_7, gk_8, gk_9, gk_10, gk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * gi0_5[k]
                 - f_10 * gi1_5[k]
                 + pb_y[k] * gk_8[k];

        t_9[k] = f_5 * gi0_7[k]
                 - f_6 * gi1_6[k]
                 + pb_y[k] * gk_9[k];

        t_10[k] = f_3 * gi0_8[k]
                  - f_4 * gi1_7[k]
                  + pb_y[k] * gk_10[k];

        t_11[k] = f_9 * gi0_8[k]
                  - f_10 * gi1_7[k]
                  + pb_z[k] * gk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, gi0_9, gi0_11, gi0_12, gi1_8, gi1_9, gi1_10, \
                         gk_12, gk_13, gk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * gi0_9[k]
                  - f_12 * gi1_8[k]
                  + pb_y[k] * gk_12[k];

        t_13[k] = f_7 * gi0_11[k]
                  - f_8 * gi1_9[k]
                  + pb_y[k] * gk_13[k];

        t_14[k] = f_5 * gi0_12[k]
                  - f_6 * gi1_10[k]
                  + pb_y[k] * gk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, gi0_13, gi0_14, gi0_16, gi1_11, \
                         gi1_12, gi1_14, gk_15, gk_16, gk_17, gk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * gi0_13[k]
                  - f_4 * gi1_11[k]
                  + pb_y[k] * gk_15[k];

        t_16[k] = f_11 * gi0_13[k]
                  - f_12 * gi1_11[k]
                  + pb_z[k] * gk_16[k];

        t_17[k] = f_1 * gi0_14[k]
                  - f_2 * gi1_12[k]
                  + pb_y[k] * gk_17[k];

        t_18[k] = f_11 * gi0_16[k]
                  - f_12 * gi1_14[k]
                  + pb_y[k] * gk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, gi0_17, gi0_18, gi0_19, gi1_15, gi1_16, \
                         gi1_17, gk_19, gk_20, gk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_9 * gi0_17[k]
                  - f_10 * gi1_15[k]
                  + pb_y[k] * gk_19[k];

        t_20[k] = f_7 * gi0_18[k]
                  - f_8 * gi1_16[k]
                  + pb_y[k] * gk_20[k];

        t_21[k] = f_5 * gi0_19[k]
                  - f_6 * gi1_17[k]
                  + pb_y[k] * gk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pa_z, pb_y, pb_z, fl_0, gi0_20, gi1_18, \
                         gk_22, gk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gi0_20[k]
                  - f_4 * gi1_18[k]
                  + pb_y[k] * gk_22[k];

        t_23[k] = f_1 * gi0_20[k]
                  - f_2 * gi1_18[k]
                  + pb_z[k] * gk_23[k];

        t_24[k] = pa_y[k] * fl_0[k];

        t_25[k] = pa_z[k] * fl_0[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, dl0_0, dl1_0, fk_17, fk_18, fl_1, \
                         gi0_33, gi0_35, gi1_25, gi1_27, gk_27, gk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_13 * dl0_0[k]
                  - f_14 * dl1_0[k]
                  + pa_y[k] * fl_1[k];

        t_27[k] = f_15 * fk_17[k]
                  + f_11 * gi0_33[k]
                  - f_12 * gi1_25[k]
                  + pb_x[k] * gk_27[k];

        t_28[k] = f_15 * fk_18[k]
                  + f_9 * gi0_35[k]
                  - f_10 * gi1_27[k]
                  + pb_x[k] * gk_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, fk_19, fk_20, fk_21, gi0_37, gi0_39, gi0_40, \
                         gi1_29, gi1_31, gi1_32, gk_29, gk_30, gk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_15 * fk_19[k]
                  + f_7 * gi0_37[k]
                  - f_8 * gi1_29[k]
                  + pb_x[k] * gk_29[k];

        t_30[k] = f_15 * fk_20[k]
                  + f_5 * gi0_39[k]
                  - f_6 * gi1_31[k]
                  + pb_x[k] * gk_30[k];

        t_31[k] = f_15 * fk_21[k]
                  + f_3 * gi0_40[k]
                  - f_4 * gi1_32[k]
                  + pb_x[k] * gk_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pa_z, pb_x, dl0_0, dl0_1, dl1_0, dl1_1, \
                         fk_23, fl_2, fl_3, gi0_51, gi1_41, gk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * dl0_1[k]
                  - f_14 * dl1_1[k]
                  + pa_x[k] * fl_3[k];

        t_33[k] = f_13 * dl0_0[k]
                  - f_14 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_34[k] = f_15 * fk_23[k]
                  + f_11 * gi0_51[k]
                  - f_12 * gi1_41[k]
                  + pb_x[k] * gk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, fk_24, fk_25, fk_26, gi0_53, gi0_55, gi0_56, \
                         gi1_43, gi1_45, gi1_46, gk_35, gk_36, gk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_15 * fk_24[k]
                  + f_9 * gi0_53[k]
                  - f_10 * gi1_43[k]
                  + pb_x[k] * gk_35[k];

        t_36[k] = f_15 * fk_25[k]
                  + f_7 * gi0_55[k]
                  - f_8 * gi1_45[k]
                  + pb_x[k] * gk_36[k];

        t_37[k] = f_15 * fk_26[k]
                  + f_5 * gi0_56[k]
                  - f_6 * gi1_46[k]
                  + pb_x[k] * gk_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_x, dl0_2, dl1_2, fk_27, fl_4, fl_5, \
                         fl_8, gi0_62, gi1_52, gk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_15 * fk_27[k]
                  + f_3 * gi0_62[k]
                  - f_4 * gi1_52[k]
                  + pb_x[k] * gk_38[k];

        t_39[k] = f_13 * dl0_2[k]
                  - f_14 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_40[k] = pa_x[k] * fl_5[k];

        t_41[k] = pa_x[k] * fl_8[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, gi0_76, gi0_78, gi0_79, gi1_65, gi1_67, \
                         gi1_68, gk_42, gk_43, gk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * gi0_76[k]
                  - f_2 * gi1_65[k]
                  + pb_x[k] * gk_42[k];

        t_43[k] = f_11 * gi0_78[k]
                  - f_12 * gi1_67[k]
                  + pb_x[k] * gk_43[k];

        t_44[k] = f_11 * gi0_79[k]
                  - f_12 * gi1_68[k]
                  + pb_x[k] * gk_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, gi0_80, gi0_82, gi0_83, gi1_69, gi1_70, \
                         gi1_71, gk_45, gk_46, gk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * gi0_80[k]
                  - f_10 * gi1_69[k]
                  + pb_x[k] * gk_45[k];

        t_46[k] = f_9 * gi0_82[k]
                  - f_10 * gi1_70[k]
                  + pb_x[k] * gk_46[k];

        t_47[k] = f_7 * gi0_83[k]
                  - f_8 * gi1_71[k]
                  + pb_x[k] * gk_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, gi0_85, gi0_86, gi0_87, gi1_72, gi1_73, \
                         gi1_74, gk_48, gk_49, gk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_7 * gi0_85[k]
                  - f_8 * gi1_72[k]
                  + pb_x[k] * gk_48[k];

        t_49[k] = f_7 * gi0_86[k]
                  - f_8 * gi1_73[k]
                  + pb_x[k] * gk_49[k];

        t_50[k] = f_5 * gi0_87[k]
                  - f_6 * gi1_74[k]
                  + pb_x[k] * gk_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, gi0_88, gi0_89, gi0_90, gi1_75, gi1_76, \
                         gi1_77, gk_51, gk_52, gk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * gi0_88[k]
                  - f_6 * gi1_75[k]
                  + pb_x[k] * gk_51[k];

        t_52[k] = f_5 * gi0_89[k]
                  - f_6 * gi1_76[k]
                  + pb_x[k] * gk_52[k];

        t_53[k] = f_5 * gi0_90[k]
                  - f_6 * gi1_77[k]
                  + pb_x[k] * gk_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, gi0_91, gi0_93, gi0_94, gi1_78, gi1_80, \
                         gi1_81, gk_54, gk_55, gk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * gi0_91[k]
                  - f_4 * gi1_78[k]
                  + pb_x[k] * gk_54[k];

        t_55[k] = f_3 * gi0_93[k]
                  - f_4 * gi1_80[k]
                  + pb_x[k] * gk_55[k];

        t_56[k] = f_3 * gi0_94[k]
                  - f_4 * gi1_81[k]
                  + pb_x[k] * gk_56[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, pb_y, fk_38, gi0_91, gi0_95, gi0_96, gi1_78, \
                         gi1_82, gi1_83, gk_57, gk_58, gk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * gi0_95[k]
                  - f_4 * gi1_82[k]
                  + pb_x[k] * gk_57[k];

        t_58[k] = f_3 * gi0_96[k]
                  - f_4 * gi1_83[k]
                  + pb_x[k] * gk_58[k];

        t_59[k] = f_0 * fk_38[k]
                  + f_1 * gi0_91[k]
                  - f_2 * gi1_78[k]
                  + pb_y[k] * gk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_z, gi0_91, gi0_92, gi0_93, gi1_78, gi1_79, \
                         gi1_80, gk_60, gk_61, gk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * gi0_91[k]
                  - f_4 * gi1_78[k]
                  + pb_z[k] * gk_60[k];

        t_61[k] = f_5 * gi0_92[k]
                  - f_6 * gi1_79[k]
                  + pb_z[k] * gk_61[k];

        t_62[k] = f_7 * gi0_93[k]
                  - f_8 * gi1_80[k]
                  + pb_z[k] * gk_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_z, fl_5, gi0_94, gi0_95, gi0_96, \
                         gi1_81, gi1_82, gi1_83, gk_63, gk_64, gk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * gi0_94[k]
                  - f_10 * gi1_81[k]
                  + pb_z[k] * gk_63[k];

        t_64[k] = f_11 * gi0_95[k]
                  - f_12 * gi1_82[k]
                  + pb_z[k] * gk_64[k];

        t_65[k] = f_1 * gi0_96[k]
                  - f_2 * gi1_83[k]
                  + pb_z[k] * gk_65[k];

        t_66[k] = pa_z[k] * fl_5[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_z, pb_y, dl0_1, dl1_1, fk_45, fk_46, fl_6, \
                         gi0_117, gi0_118, gi1_99, gi1_100, gk_68, \
                         gk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * dl0_1[k]
                  - f_14 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_68[k] = f_15 * fk_45[k]
                  + f_11 * gi0_117[k]
                  - f_12 * gi1_99[k]
                  + pb_y[k] * gk_68[k];

        t_69[k] = f_15 * fk_46[k]
                  + f_9 * gi0_118[k]
                  - f_10 * gi1_100[k]
                  + pb_y[k] * gk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_y, fk_47, fk_48, fk_49, gi0_119, gi0_120, \
                         gi0_121, gi1_101, gi1_102, gi1_103, gk_70, gk_71, \
                         gk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_15 * fk_47[k]
                  + f_7 * gi0_119[k]
                  - f_8 * gi1_101[k]
                  + pb_y[k] * gk_70[k];

        t_71[k] = f_15 * fk_48[k]
                  + f_5 * gi0_120[k]
                  - f_6 * gi1_102[k]
                  + pb_y[k] * gk_71[k];

        t_72[k] = f_15 * fk_49[k]
                  + f_3 * gi0_121[k]
                  - f_4 * gi1_103[k]
                  + pb_y[k] * gk_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pb_x, dl0_2, dl1_2, fl_7, fl_8, \
                         gi0_132, gi0_134, gi1_110, gi1_112, gk_75, \
                         gk_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_13 * dl0_2[k]
                  - f_14 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_74[k] = pa_y[k] * fl_8[k];

        t_75[k] = f_1 * gi0_132[k]
                  - f_2 * gi1_110[k]
                  + pb_x[k] * gk_75[k];

        t_76[k] = f_11 * gi0_134[k]
                  - f_12 * gi1_112[k]
                  + pb_x[k] * gk_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_x, gi0_135, gi0_136, gi0_138, gi1_113, gi1_114, \
                         gi1_115, gk_77, gk_78, gk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_11 * gi0_135[k]
                  - f_12 * gi1_113[k]
                  + pb_x[k] * gk_77[k];

        t_78[k] = f_9 * gi0_136[k]
                  - f_10 * gi1_114[k]
                  + pb_x[k] * gk_78[k];

        t_79[k] = f_9 * gi0_138[k]
                  - f_10 * gi1_115[k]
                  + pb_x[k] * gk_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, gi0_139, gi0_140, gi0_142, gi1_116, gi1_117, \
                         gi1_118, gk_80, gk_81, gk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_7 * gi0_139[k]
                  - f_8 * gi1_116[k]
                  + pb_x[k] * gk_80[k];

        t_81[k] = f_7 * gi0_140[k]
                  - f_8 * gi1_117[k]
                  + pb_x[k] * gk_81[k];

        t_82[k] = f_7 * gi0_142[k]
                  - f_8 * gi1_118[k]
                  + pb_x[k] * gk_82[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, gi0_143, gi0_144, gi0_145, gi1_119, gi1_120, \
                         gi1_121, gk_83, gk_84, gk_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_5 * gi0_143[k]
                  - f_6 * gi1_119[k]
                  + pb_x[k] * gk_83[k];

        t_84[k] = f_5 * gi0_144[k]
                  - f_6 * gi1_120[k]
                  + pb_x[k] * gk_84[k];

        t_85[k] = f_5 * gi0_145[k]
                  - f_6 * gi1_121[k]
                  + pb_x[k] * gk_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_x, gi0_146, gi0_147, gi0_148, gi1_122, gi1_123, \
                         gi1_124, gk_86, gk_87, gk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_5 * gi0_146[k]
                  - f_6 * gi1_122[k]
                  + pb_x[k] * gk_86[k];

        t_87[k] = f_3 * gi0_147[k]
                  - f_4 * gi1_123[k]
                  + pb_x[k] * gk_87[k];

        t_88[k] = f_3 * gi0_148[k]
                  - f_4 * gi1_124[k]
                  + pb_x[k] * gk_88[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pb_x, gi0_149, gi0_150, gi0_152, gi1_125, gi1_126, \
                         gi1_128, gk_89, gk_90, gk_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_3 * gi0_149[k]
                  - f_4 * gi1_125[k]
                  + pb_x[k] * gk_89[k];

        t_90[k] = f_3 * gi0_150[k]
                  - f_4 * gi1_126[k]
                  + pb_x[k] * gk_90[k];

        t_91[k] = f_3 * gi0_152[k]
                  - f_4 * gi1_128[k]
                  + pb_x[k] * gk_91[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pb_y, gi0_147, gi0_148, gi0_149, gi1_123, gi1_124, \
                         gi1_125, gk_92, gk_93, gk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_1 * gi0_147[k]
                  - f_2 * gi1_123[k]
                  + pb_y[k] * gk_92[k];

        t_93[k] = f_11 * gi0_148[k]
                  - f_12 * gi1_124[k]
                  + pb_y[k] * gk_93[k];

        t_94[k] = f_9 * gi0_149[k]
                  - f_10 * gi1_125[k]
                  + pb_y[k] * gk_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_y, gi0_150, gi0_151, gi0_152, gi1_126, gi1_127, \
                         gi1_128, gk_95, gk_96, gk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_7 * gi0_150[k]
                  - f_8 * gi1_126[k]
                  + pb_y[k] * gk_95[k];

        t_96[k] = f_5 * gi0_151[k]
                  - f_6 * gi1_127[k]
                  + pb_y[k] * gk_96[k];

        t_97[k] = f_3 * gi0_152[k]
                  - f_4 * gi1_128[k]
                  + pb_y[k] * gk_97[k];
    }

#pragma omp simd aligned(t_98, pb_z, fk_65, gi0_152, gi1_128, gk_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * fk_65[k]
                  + f_1 * gi0_152[k]
                  - f_2 * gi1_128[k]
                  + pb_z[k] * gk_98[k];
    }
}

auto
compute_prim_gl_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_26 = buffer.data(gi0 + 26);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_25 = buffer.data(gi1 + 25);
    const auto *gi1_28 = buffer.data(gi1 + 28);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_30 = buffer.data(gi1 + 30);
    const auto *gi1_31 = buffer.data(gi1 + 31);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_38 = buffer.data(gi1 + 38);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_44 = buffer.data(gk + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_4, \
                         gi1_5, gi1_6, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_4[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_5[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_6[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_7, gi1_8, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_7[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_8[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_10, gi1_11, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_10[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_11[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_12, gi1_13, gi1_14, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_12[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_13[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_14[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_19, fl_4, fl_5, \
                         fl_8, gi0_17, gi1_25, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_19[k]
                  + f_1 * gi0_17[k]
                  - f_2 * gi1_25[k]
                  + pb_y[k] * gk_29[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_21, fl_5, fl_6, \
                         gi0_20, gi1_28, gk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_21[k]
                  + f_6 * gi0_20[k]
                  - f_7 * gi1_28[k]
                  + pb_y[k] * gk_32[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_22, fk_23, fk_24, gi0_21, gi0_22, gi0_23, \
                         gi1_29, gi1_30, gi1_31, gk_33, gk_34, gk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_22[k]
                  + f_8 * gi0_21[k]
                  - f_9 * gi1_29[k]
                  + pb_y[k] * gk_33[k];

        t_24[k] = f_5 * fk_23[k]
                  + f_10 * gi0_22[k]
                  - f_11 * gi1_30[k]
                  + pb_y[k] * gk_34[k];

        t_25[k] = f_5 * fk_24[k]
                  + f_12 * gi0_23[k]
                  - f_13 * gi1_31[k]
                  + pb_y[k] * gk_35[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_25, fl_7, fl_8, \
                         gi0_24, gi1_32, gk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_25[k]
                  + f_14 * gi0_24[k]
                  - f_15 * gi1_32[k]
                  + pb_y[k] * gk_36[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_35, gi0_26, gi1_38, gk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_35[k]
                  + f_1 * gi0_26[k]
                  - f_2 * gi1_38[k]
                  + pb_z[k] * gk_44[k];
    }
}

auto
compute_prim_gl_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_25 = buffer.data(gi0 + 25);
    const auto *gi0_28 = buffer.data(gi0 + 28);
    const auto *gi0_29 = buffer.data(gi0 + 29);
    const auto *gi0_30 = buffer.data(gi0 + 30);
    const auto *gi0_31 = buffer.data(gi0 + 31);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_38 = buffer.data(gi0 + 38);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_20 = buffer.data(gi1 + 20);
    const auto *gi1_21 = buffer.data(gi1 + 21);
    const auto *gi1_22 = buffer.data(gi1 + 22);
    const auto *gi1_23 = buffer.data(gi1 + 23);
    const auto *gi1_24 = buffer.data(gi1 + 24);
    const auto *gi1_27 = buffer.data(gi1 + 27);
    const auto *gi1_28 = buffer.data(gi1 + 28);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_30 = buffer.data(gi1 + 30);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_54 = buffer.data(gi1 + 54);
    const auto *gi1_65 = buffer.data(gi1 + 65);
    const auto *gi1_66 = buffer.data(gi1 + 66);
    const auto *gi1_67 = buffer.data(gi1 + 67);
    const auto *gi1_68 = buffer.data(gi1 + 68);
    const auto *gi1_69 = buffer.data(gi1 + 69);
    const auto *gi1_90 = buffer.data(gi1 + 90);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_86 = buffer.data(gk + 86);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_20, \
                         gi1_21, gi1_22, gk_18, gk_19, gk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_20[k]
                 + pb_x[k] * gk_18[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_21[k]
                 + pb_x[k] * gk_19[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_22[k]
                 + pb_x[k] * gk_20[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_23, gi1_24, gk_21, gk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_23[k]
                 + pb_x[k] * gk_21[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_24[k]
                 + pb_x[k] * gk_22[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_27, gi1_28, gk_25, gk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_27[k]
                  + pb_x[k] * gk_25[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_28[k]
                  + pb_x[k] * gk_26[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_29, gi1_30, gi1_32, gk_27, gk_28, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_29[k]
                  + pb_x[k] * gk_27[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_30[k]
                  + pb_x[k] * gk_28[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_32[k]
                  + pb_x[k] * gk_29[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_19, fl_4, fl_5, \
                         fl_8, gi0_25, gi1_54, gk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_19[k]
                  + f_1 * gi0_25[k]
                  - f_2 * gi1_54[k]
                  + pb_y[k] * gk_52[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_21, fl_5, fl_6, \
                         gi0_28, gi1_65, gk_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_21[k]
                  + f_6 * gi0_28[k]
                  - f_7 * gi1_65[k]
                  + pb_y[k] * gk_60[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_22, fk_23, fk_24, gi0_29, gi0_30, gi0_31, \
                         gi1_66, gi1_67, gi1_68, gk_61, gk_62, gk_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_22[k]
                  + f_8 * gi0_29[k]
                  - f_9 * gi1_66[k]
                  + pb_y[k] * gk_61[k];

        t_24[k] = f_5 * fk_23[k]
                  + f_10 * gi0_30[k]
                  - f_11 * gi1_67[k]
                  + pb_y[k] * gk_62[k];

        t_25[k] = f_5 * fk_24[k]
                  + f_12 * gi0_31[k]
                  - f_13 * gi1_68[k]
                  + pb_y[k] * gk_63[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_25, fl_7, fl_8, \
                         gi0_32, gi1_69, gk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_25[k]
                  + f_14 * gi0_32[k]
                  - f_15 * gi1_69[k]
                  + pb_y[k] * gk_64[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_35, gi0_38, gi1_90, gk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_35[k]
                  + f_1 * gi0_38[k]
                  - f_2 * gi1_90[k]
                  + pb_z[k] * gk_86[k];
    }
}

auto
compute_prim_gl_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_27 = buffer.data(gi0 + 27);
    const auto *gi0_28 = buffer.data(gi0 + 28);
    const auto *gi0_29 = buffer.data(gi0 + 29);
    const auto *gi0_30 = buffer.data(gi0 + 30);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_54 = buffer.data(gi0 + 54);
    const auto *gi0_65 = buffer.data(gi0 + 65);
    const auto *gi0_66 = buffer.data(gi0 + 66);
    const auto *gi0_67 = buffer.data(gi0 + 67);
    const auto *gi0_68 = buffer.data(gi0 + 68);
    const auto *gi0_69 = buffer.data(gi0 + 69);
    const auto *gi0_90 = buffer.data(gi0 + 90);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_18 = buffer.data(gi1 + 18);
    const auto *gi1_19 = buffer.data(gi1 + 19);
    const auto *gi1_20 = buffer.data(gi1 + 20);
    const auto *gi1_21 = buffer.data(gi1 + 21);
    const auto *gi1_22 = buffer.data(gi1 + 22);
    const auto *gi1_24 = buffer.data(gi1 + 24);
    const auto *gi1_25 = buffer.data(gi1 + 25);
    const auto *gi1_26 = buffer.data(gi1 + 26);
    const auto *gi1_27 = buffer.data(gi1 + 27);
    const auto *gi1_28 = buffer.data(gi1 + 28);
    const auto *gi1_48 = buffer.data(gi1 + 48);
    const auto *gi1_56 = buffer.data(gi1 + 56);
    const auto *gi1_57 = buffer.data(gi1 + 57);
    const auto *gi1_58 = buffer.data(gi1 + 58);
    const auto *gi1_59 = buffer.data(gi1 + 59);
    const auto *gi1_60 = buffer.data(gi1 + 60);
    const auto *gi1_80 = buffer.data(gi1 + 80);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_29 = buffer.data(gk + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_20, gi0_21, gi0_22, \
                         gi1_18, gi1_19, gi1_20, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_20[k]
                 - f_7 * gi1_18[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_21[k]
                 - f_9 * gi1_19[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_22[k]
                 - f_11 * gi1_20[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_23, \
                         gi0_24, gi1_21, gi1_22, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_23[k]
                 - f_13 * gi1_21[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_24[k]
                 - f_15 * gi1_22[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_27, gi0_28, gi1_24, gi1_25, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_27[k]
                  - f_7 * gi1_24[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_28[k]
                  - f_9 * gi1_25[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_29, gi0_30, gi0_32, \
                         gi1_26, gi1_27, gi1_28, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_29[k]
                  - f_11 * gi1_26[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_30[k]
                  - f_13 * gi1_27[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_32[k]
                  - f_15 * gi1_28[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_19, fl_4, fl_5, \
                         fl_8, gi0_54, gi1_48, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_19[k]
                  + f_1 * gi0_54[k]
                  - f_2 * gi1_48[k]
                  + pb_y[k] * gk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_21, fl_5, fl_6, \
                         gi0_65, gi1_56, gk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_21[k]
                  + f_6 * gi0_65[k]
                  - f_7 * gi1_56[k]
                  + pb_y[k] * gk_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_22, fk_23, fk_24, gi0_66, gi0_67, gi0_68, \
                         gi1_57, gi1_58, gi1_59, gk_23, gk_24, gk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_22[k]
                  + f_8 * gi0_66[k]
                  - f_9 * gi1_57[k]
                  + pb_y[k] * gk_23[k];

        t_24[k] = f_5 * fk_23[k]
                  + f_10 * gi0_67[k]
                  - f_11 * gi1_58[k]
                  + pb_y[k] * gk_24[k];

        t_25[k] = f_5 * fk_24[k]
                  + f_12 * gi0_68[k]
                  - f_13 * gi1_59[k]
                  + pb_y[k] * gk_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_25, fl_7, fl_8, \
                         gi0_69, gi1_60, gk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_25[k]
                  + f_14 * gi0_69[k]
                  - f_15 * gi1_60[k]
                  + pb_y[k] * gk_26[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_35, gi0_90, gi1_80, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_35[k]
                  + f_1 * gi0_90[k]
                  - f_2 * gi1_80[k]
                  + pb_z[k] * gk_29[k];
    }
}

auto
compute_prim_gl_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_23 = buffer.data(fk + 23);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_17 = buffer.data(gi0 + 17);
    const auto *gi0_20 = buffer.data(gi0 + 20);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_22 = buffer.data(gi0 + 22);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_26 = buffer.data(gi0 + 26);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_33 = buffer.data(gi1 + 33);
    const auto *gi1_34 = buffer.data(gi1 + 34);
    const auto *gi1_35 = buffer.data(gi1 + 35);
    const auto *gi1_36 = buffer.data(gi1 + 36);
    const auto *gi1_50 = buffer.data(gi1 + 50);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_56 = buffer.data(gk + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_4, \
                         gi1_5, gi1_6, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_4[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_5[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_6[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_7, gi1_8, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_7[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_8[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_10, gi1_11, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_10[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_11[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_12, gi1_13, gi1_14, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_12[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_13[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_14[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_15, fl_4, fl_5, \
                         fl_8, gi0_17, gi1_29, gk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_15[k]
                  + f_1 * gi0_17[k]
                  - f_2 * gi1_29[k]
                  + pb_y[k] * gk_33[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_17, fl_5, fl_6, \
                         gi0_20, gi1_32, gk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_17[k]
                  + f_6 * gi0_20[k]
                  - f_7 * gi1_32[k]
                  + pb_y[k] * gk_36[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_18, fk_19, fk_20, gi0_21, gi0_22, gi0_23, \
                         gi1_33, gi1_34, gi1_35, gk_37, gk_38, gk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_18[k]
                  + f_8 * gi0_21[k]
                  - f_9 * gi1_33[k]
                  + pb_y[k] * gk_37[k];

        t_24[k] = f_5 * fk_19[k]
                  + f_10 * gi0_22[k]
                  - f_11 * gi1_34[k]
                  + pb_y[k] * gk_38[k];

        t_25[k] = f_5 * fk_20[k]
                  + f_12 * gi0_23[k]
                  - f_13 * gi1_35[k]
                  + pb_y[k] * gk_39[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_21, fl_7, fl_8, \
                         gi0_24, gi1_36, gk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_21[k]
                  + f_14 * gi0_24[k]
                  - f_15 * gi1_36[k]
                  + pb_y[k] * gk_40[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_23, gi0_26, gi1_50, gk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_23[k]
                  + f_1 * gi0_26[k]
                  - f_2 * gi1_50[k]
                  + pb_z[k] * gk_56[k];
    }
}

auto
compute_prim_gl_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dl0, const size_t dl1,
                                     const size_t fk, const size_t fl, const size_t gi0,
                                     const size_t gi1, const size_t gk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_29 = buffer.data(gi0 + 29);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_33 = buffer.data(gi0 + 33);
    const auto *gi0_34 = buffer.data(gi0 + 34);
    const auto *gi0_35 = buffer.data(gi0 + 35);
    const auto *gi0_36 = buffer.data(gi0 + 36);
    const auto *gi0_50 = buffer.data(gi0 + 50);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_33 = buffer.data(gi1 + 33);
    const auto *gi1_34 = buffer.data(gi1 + 34);
    const auto *gi1_35 = buffer.data(gi1 + 35);
    const auto *gi1_36 = buffer.data(gi1 + 36);
    const auto *gi1_50 = buffer.data(gi1 + 50);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_56 = buffer.data(gk + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_4, \
                         gi1_5, gi1_6, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_4[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_5[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_6[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_7, gi1_8, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_7[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_8[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_10, gi1_11, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_10[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_11[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_12, gi1_13, gi1_14, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_12[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_13[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_14[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_19, fl_4, fl_5, \
                         fl_8, gi0_29, gi1_29, gk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_19[k]
                  + f_1 * gi0_29[k]
                  - f_2 * gi1_29[k]
                  + pb_y[k] * gk_33[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_21, fl_5, fl_6, \
                         gi0_32, gi1_32, gk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_21[k]
                  + f_6 * gi0_32[k]
                  - f_7 * gi1_32[k]
                  + pb_y[k] * gk_36[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_22, fk_23, fk_24, gi0_33, gi0_34, gi0_35, \
                         gi1_33, gi1_34, gi1_35, gk_37, gk_38, gk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_22[k]
                  + f_8 * gi0_33[k]
                  - f_9 * gi1_33[k]
                  + pb_y[k] * gk_37[k];

        t_24[k] = f_5 * fk_23[k]
                  + f_10 * gi0_34[k]
                  - f_11 * gi1_34[k]
                  + pb_y[k] * gk_38[k];

        t_25[k] = f_5 * fk_24[k]
                  + f_12 * gi0_35[k]
                  - f_13 * gi1_35[k]
                  + pb_y[k] * gk_39[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_25, fl_7, fl_8, \
                         gi0_36, gi1_36, gk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_25[k]
                  + f_14 * gi0_36[k]
                  - f_15 * gi1_36[k]
                  + pb_y[k] * gk_40[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_35, gi0_50, gi1_50, gk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_35[k]
                  + f_1 * gi0_50[k]
                  - f_2 * gi1_50[k]
                  + pb_z[k] * gk_56[k];
    }
}

auto
compute_prim_gl_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dl0, const size_t dl1,
                                      const size_t fk, const size_t fl, const size_t gi0,
                                      const size_t gi1, const size_t gk, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_4 = buffer.data(gi0 + 4);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_11 = buffer.data(gi0 + 11);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_29 = buffer.data(gi0 + 29);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_33 = buffer.data(gi0 + 33);
    const auto *gi0_34 = buffer.data(gi0 + 34);
    const auto *gi0_35 = buffer.data(gi0 + 35);
    const auto *gi0_36 = buffer.data(gi0 + 36);
    const auto *gi0_50 = buffer.data(gi0 + 50);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_9 = buffer.data(gi1 + 9);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_15 = buffer.data(gi1 + 15);
    const auto *gi1_16 = buffer.data(gi1 + 16);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_36 = buffer.data(gi1 + 36);
    const auto *gi1_37 = buffer.data(gi1 + 37);
    const auto *gi1_38 = buffer.data(gi1 + 38);
    const auto *gi1_39 = buffer.data(gi1 + 39);
    const auto *gi1_40 = buffer.data(gi1 + 40);
    const auto *gi1_56 = buffer.data(gi1 + 56);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_56 = buffer.data(gk + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_4, gi0_5, gi0_6, gi1_6, \
                         gi1_7, gi1_8, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_4[k]
                 - f_7 * gi1_6[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_5[k]
                 - f_9 * gi1_7[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_6[k]
                 - f_11 * gi1_8[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_7, \
                         gi0_8, gi1_9, gi1_10, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_7[k]
                 - f_13 * gi1_9[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_8[k]
                 - f_15 * gi1_10[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_10, gi0_11, gi1_12, gi1_13, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_10[k]
                  - f_7 * gi1_12[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_11[k]
                  - f_9 * gi1_13[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_12, gi0_13, gi0_14, \
                         gi1_14, gi1_15, gi1_16, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_12[k]
                  - f_11 * gi1_14[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_13[k]
                  - f_13 * gi1_15[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_14[k]
                  - f_15 * gi1_16[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_19, fl_4, fl_5, \
                         fl_8, gi0_29, gi1_32, gk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_19[k]
                  + f_1 * gi0_29[k]
                  - f_2 * gi1_32[k]
                  + pb_y[k] * gk_33[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_21, fl_5, fl_6, \
                         gi0_32, gi1_36, gk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_21[k]
                  + f_6 * gi0_32[k]
                  - f_7 * gi1_36[k]
                  + pb_y[k] * gk_36[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_22, fk_23, fk_24, gi0_33, gi0_34, gi0_35, \
                         gi1_37, gi1_38, gi1_39, gk_37, gk_38, gk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_22[k]
                  + f_8 * gi0_33[k]
                  - f_9 * gi1_37[k]
                  + pb_y[k] * gk_37[k];

        t_24[k] = f_5 * fk_23[k]
                  + f_10 * gi0_34[k]
                  - f_11 * gi1_38[k]
                  + pb_y[k] * gk_38[k];

        t_25[k] = f_5 * fk_24[k]
                  + f_12 * gi0_35[k]
                  - f_13 * gi1_39[k]
                  + pb_y[k] * gk_39[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_25, fl_7, fl_8, \
                         gi0_36, gi1_40, gk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_25[k]
                  + f_14 * gi0_36[k]
                  - f_15 * gi1_40[k]
                  + pb_y[k] * gk_40[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_35, gi0_50, gi1_56, gk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_35[k]
                  + f_1 * gi0_50[k]
                  - f_2 * gi1_56[k]
                  + pb_z[k] * gk_56[k];
    }
}

auto
compute_prim_gl_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dl0, const size_t dl1,
                                      const size_t fk, const size_t fl, const size_t gi0,
                                      const size_t gi1, const size_t gk, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.5 / beta;
    const auto f_7 = 2.5 * alpha / (beta * p);
    const auto f_8 = 2.0 / beta;
    const auto f_9 = 2.0 * alpha / (beta * p);
    const auto f_10 = 1.5 / beta;
    const auto f_11 = 1.5 * alpha / (beta * p);
    const auto f_12 = 1.0 / beta;
    const auto f_13 = alpha / (beta * p);
    const auto f_14 = 0.5 / beta;
    const auto f_15 = 0.5 * alpha / (beta * p);

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

    const auto *dl0_0 = buffer.data(dl0 + 0);
    const auto *dl0_1 = buffer.data(dl0 + 1);
    const auto *dl0_2 = buffer.data(dl0 + 2);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_1 = buffer.data(dl1 + 1);
    const auto *dl1_2 = buffer.data(dl1 + 2);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_35 = buffer.data(fk + 35);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_7 = buffer.data(gi0 + 7);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_9 = buffer.data(gi0 + 9);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_15 = buffer.data(gi0 + 15);
    const auto *gi0_16 = buffer.data(gi0 + 16);
    const auto *gi0_32 = buffer.data(gi0 + 32);
    const auto *gi0_36 = buffer.data(gi0 + 36);
    const auto *gi0_37 = buffer.data(gi0 + 37);
    const auto *gi0_38 = buffer.data(gi0 + 38);
    const auto *gi0_39 = buffer.data(gi0 + 39);
    const auto *gi0_40 = buffer.data(gi0 + 40);
    const auto *gi0_56 = buffer.data(gi0 + 56);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_4 = buffer.data(gi1 + 4);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_7 = buffer.data(gi1 + 7);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_11 = buffer.data(gi1 + 11);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_29 = buffer.data(gi1 + 29);
    const auto *gi1_32 = buffer.data(gi1 + 32);
    const auto *gi1_33 = buffer.data(gi1 + 33);
    const auto *gi1_34 = buffer.data(gi1 + 34);
    const auto *gi1_35 = buffer.data(gi1 + 35);
    const auto *gi1_36 = buffer.data(gi1 + 36);
    const auto *gi1_50 = buffer.data(gi1 + 50);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_29 = buffer.data(gk + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dl0_0, dl1_0, fk_0, fl_0, fl_1, \
                         gi0_0, gi1_0, gk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_0[k]
                 + f_1 * gi0_0[k]
                 - f_2 * gi1_0[k]
                 + pb_x[k] * gk_0[k];

        t_1[k] = pa_y[k] * fl_0[k];

        t_2[k] = pa_z[k] * fl_0[k];

        t_3[k] = f_3 * dl0_0[k]
                 - f_4 * dl1_0[k]
                 + pa_y[k] * fl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fk_3, fk_4, fk_5, gi0_6, gi0_7, gi0_8, gi1_4, \
                         gi1_5, gi1_6, gk_4, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fk_3[k]
                 + f_6 * gi0_6[k]
                 - f_7 * gi1_4[k]
                 + pb_x[k] * gk_4[k];

        t_5[k] = f_5 * fk_4[k]
                 + f_8 * gi0_7[k]
                 - f_9 * gi1_5[k]
                 + pb_x[k] * gk_5[k];

        t_6[k] = f_5 * fk_5[k]
                 + f_10 * gi0_8[k]
                 - f_11 * gi1_6[k]
                 + pb_x[k] * gk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, dl0_1, dl1_1, fk_6, fk_7, fl_3, gi0_9, \
                         gi0_10, gi1_7, gi1_8, gk_7, gk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * fk_6[k]
                 + f_12 * gi0_9[k]
                 - f_13 * gi1_7[k]
                 + pb_x[k] * gk_7[k];

        t_8[k] = f_5 * fk_7[k]
                 + f_14 * gi0_10[k]
                 - f_15 * gi1_8[k]
                 + pb_x[k] * gk_8[k];

        t_9[k] = f_3 * dl0_1[k]
                 - f_4 * dl1_1[k]
                 + pa_x[k] * fl_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, dl0_0, dl1_0, fk_9, fk_10, fl_2, \
                         gi0_12, gi0_13, gi1_10, gi1_11, gk_11, gk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dl0_0[k]
                  - f_4 * dl1_0[k]
                  + pa_z[k] * fl_2[k];

        t_11[k] = f_5 * fk_9[k]
                  + f_6 * gi0_12[k]
                  - f_7 * gi1_10[k]
                  + pb_x[k] * gk_11[k];

        t_12[k] = f_5 * fk_10[k]
                  + f_8 * gi0_13[k]
                  - f_9 * gi1_11[k]
                  + pb_x[k] * gk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fk_11, fk_12, fk_13, gi0_14, gi0_15, gi0_16, \
                         gi1_12, gi1_13, gi1_14, gk_13, gk_14, gk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fk_11[k]
                  + f_10 * gi0_14[k]
                  - f_11 * gi1_12[k]
                  + pb_x[k] * gk_13[k];

        t_14[k] = f_5 * fk_12[k]
                  + f_12 * gi0_15[k]
                  - f_13 * gi1_13[k]
                  + pb_x[k] * gk_14[k];

        t_15[k] = f_5 * fk_13[k]
                  + f_14 * gi0_16[k]
                  - f_15 * gi1_14[k]
                  + pb_x[k] * gk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dl0_2, dl1_2, fk_19, fl_4, fl_5, \
                         fl_8, gi0_32, gi1_29, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_x[k] * fl_4[k];

        t_17[k] = pa_x[k] * fl_5[k];

        t_18[k] = pa_x[k] * fl_8[k];

        t_19[k] = f_0 * fk_19[k]
                  + f_1 * gi0_32[k]
                  - f_2 * gi1_29[k]
                  + pb_y[k] * gk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dl0_1, dl1_1, fk_21, fl_5, fl_6, \
                         gi0_36, gi1_32, gk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fl_5[k];

        t_21[k] = f_3 * dl0_1[k]
                  - f_4 * dl1_1[k]
                  + pa_z[k] * fl_6[k];

        t_22[k] = f_5 * fk_21[k]
                  + f_6 * gi0_36[k]
                  - f_7 * gi1_32[k]
                  + pb_y[k] * gk_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, fk_22, fk_23, fk_24, gi0_37, gi0_38, gi0_39, \
                         gi1_33, gi1_34, gi1_35, gk_23, gk_24, gk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fk_22[k]
                  + f_8 * gi0_37[k]
                  - f_9 * gi1_33[k]
                  + pb_y[k] * gk_23[k];

        t_24[k] = f_5 * fk_23[k]
                  + f_10 * gi0_38[k]
                  - f_11 * gi1_34[k]
                  + pb_y[k] * gk_24[k];

        t_25[k] = f_5 * fk_24[k]
                  + f_12 * gi0_39[k]
                  - f_13 * gi1_35[k]
                  + pb_y[k] * gk_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, dl0_2, dl1_2, fk_25, fl_7, fl_8, \
                         gi0_40, gi1_36, gk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fk_25[k]
                  + f_14 * gi0_40[k]
                  - f_15 * gi1_36[k]
                  + pb_y[k] * gk_26[k];

        t_27[k] = f_3 * dl0_2[k]
                  - f_4 * dl1_2[k]
                  + pa_y[k] * fl_7[k];

        t_28[k] = pa_y[k] * fl_8[k];
    }

#pragma omp simd aligned(t_29, pb_z, fk_35, gi0_56, gi1_50, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fk_35[k]
                  + f_1 * gi0_56[k]
                  - f_2 * gi1_50[k]
                  + pb_z[k] * gk_29[k];
    }
}

}  // namespace simdt2ceri
