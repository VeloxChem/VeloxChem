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
    const auto *dl0_171 = buffer.data(dl0 + 171);
    const auto *dl0_269 = buffer.data(dl0 + 269);

    const auto *dl1_0 = buffer.data(dl1 + 0);
    const auto *dl1_171 = buffer.data(dl1 + 171);
    const auto *dl1_269 = buffer.data(dl1 + 269);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
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
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_9 = buffer.data(fl + 9);
    const auto *fl_10 = buffer.data(fl + 10);
    const auto *fl_12 = buffer.data(fl + 12);
    const auto *fl_14 = buffer.data(fl + 14);
    const auto *fl_15 = buffer.data(fl + 15);
    const auto *fl_17 = buffer.data(fl + 17);
    const auto *fl_18 = buffer.data(fl + 18);
    const auto *fl_20 = buffer.data(fl + 20);
    const auto *fl_21 = buffer.data(fl + 21);
    const auto *fl_23 = buffer.data(fl + 23);
    const auto *fl_24 = buffer.data(fl + 24);
    const auto *fl_25 = buffer.data(fl + 25);
    const auto *fl_27 = buffer.data(fl + 27);
    const auto *fl_28 = buffer.data(fl + 28);
    const auto *fl_35 = buffer.data(fl + 35);
    const auto *fl_36 = buffer.data(fl + 36);
    const auto *fl_38 = buffer.data(fl + 38);
    const auto *fl_39 = buffer.data(fl + 39);
    const auto *fl_40 = buffer.data(fl + 40);
    const auto *fl_41 = buffer.data(fl + 41);
    const auto *fl_42 = buffer.data(fl + 42);
    const auto *fl_44 = buffer.data(fl + 44);
    const auto *fl_45 = buffer.data(fl + 45);
    const auto *fl_46 = buffer.data(fl + 46);
    const auto *fl_48 = buffer.data(fl + 48);
    const auto *fl_51 = buffer.data(fl + 51);
    const auto *fl_55 = buffer.data(fl + 55);
    const auto *fl_60 = buffer.data(fl + 60);
    const auto *fl_66 = buffer.data(fl + 66);
    const auto *fl_73 = buffer.data(fl + 73);
    const auto *fl_81 = buffer.data(fl + 81);
    const auto *fl_90 = buffer.data(fl + 90);
    const auto *fl_92 = buffer.data(fl + 92);
    const auto *fl_95 = buffer.data(fl + 95);
    const auto *fl_99 = buffer.data(fl + 99);
    const auto *fl_102 = buffer.data(fl + 102);
    const auto *fl_104 = buffer.data(fl + 104);
    const auto *fl_107 = buffer.data(fl + 107);
    const auto *fl_108 = buffer.data(fl + 108);
    const auto *fl_110 = buffer.data(fl + 110);
    const auto *fl_113 = buffer.data(fl + 113);
    const auto *fl_114 = buffer.data(fl + 114);
    const auto *fl_115 = buffer.data(fl + 115);
    const auto *fl_117 = buffer.data(fl + 117);
    const auto *fl_125 = buffer.data(fl + 125);
    const auto *fl_128 = buffer.data(fl + 128);
    const auto *fl_129 = buffer.data(fl + 129);
    const auto *fl_130 = buffer.data(fl + 130);
    const auto *fl_131 = buffer.data(fl + 131);
    const auto *fl_132 = buffer.data(fl + 132);
    const auto *fl_134 = buffer.data(fl + 134);
    const auto *fl_135 = buffer.data(fl + 135);
    const auto *fl_136 = buffer.data(fl + 136);
    const auto *fl_138 = buffer.data(fl + 138);
    const auto *fl_141 = buffer.data(fl + 141);
    const auto *fl_145 = buffer.data(fl + 145);
    const auto *fl_150 = buffer.data(fl + 150);
    const auto *fl_156 = buffer.data(fl + 156);
    const auto *fl_163 = buffer.data(fl + 163);
    const auto *fl_171 = buffer.data(fl + 171);
    const auto *fl_225 = buffer.data(fl + 225);
    const auto *fl_227 = buffer.data(fl + 227);
    const auto *fl_230 = buffer.data(fl + 230);
    const auto *fl_234 = buffer.data(fl + 234);
    const auto *fl_239 = buffer.data(fl + 239);
    const auto *fl_245 = buffer.data(fl + 245);
    const auto *fl_252 = buffer.data(fl + 252);
    const auto *fl_260 = buffer.data(fl + 260);
    const auto *fl_269 = buffer.data(fl + 269);
    const auto *fl_270 = buffer.data(fl + 270);
    const auto *fl_271 = buffer.data(fl + 271);
    const auto *fl_273 = buffer.data(fl + 273);
    const auto *fl_275 = buffer.data(fl + 275);
    const auto *fl_276 = buffer.data(fl + 276);
    const auto *fl_279 = buffer.data(fl + 279);
    const auto *fl_280 = buffer.data(fl + 280);
    const auto *fl_282 = buffer.data(fl + 282);
    const auto *fl_284 = buffer.data(fl + 284);
    const auto *fl_285 = buffer.data(fl + 285);
    const auto *fl_287 = buffer.data(fl + 287);
    const auto *fl_288 = buffer.data(fl + 288);
    const auto *fl_290 = buffer.data(fl + 290);
    const auto *fl_291 = buffer.data(fl + 291);
    const auto *fl_293 = buffer.data(fl + 293);
    const auto *fl_294 = buffer.data(fl + 294);
    const auto *fl_295 = buffer.data(fl + 295);
    const auto *fl_297 = buffer.data(fl + 297);
    const auto *fl_306 = buffer.data(fl + 306);
    const auto *fl_308 = buffer.data(fl + 308);
    const auto *fl_309 = buffer.data(fl + 309);
    const auto *fl_310 = buffer.data(fl + 310);
    const auto *fl_311 = buffer.data(fl + 311);
    const auto *fl_312 = buffer.data(fl + 312);
    const auto *fl_313 = buffer.data(fl + 313);
    const auto *fl_314 = buffer.data(fl + 314);
    const auto *fl_320 = buffer.data(fl + 320);
    const auto *fl_324 = buffer.data(fl + 324);
    const auto *fl_327 = buffer.data(fl + 327);
    const auto *fl_329 = buffer.data(fl + 329);
    const auto *fl_332 = buffer.data(fl + 332);
    const auto *fl_333 = buffer.data(fl + 333);
    const auto *fl_335 = buffer.data(fl + 335);
    const auto *fl_338 = buffer.data(fl + 338);
    const auto *fl_339 = buffer.data(fl + 339);
    const auto *fl_340 = buffer.data(fl + 340);
    const auto *fl_342 = buffer.data(fl + 342);
    const auto *fl_351 = buffer.data(fl + 351);
    const auto *fl_352 = buffer.data(fl + 352);
    const auto *fl_353 = buffer.data(fl + 353);
    const auto *fl_354 = buffer.data(fl + 354);
    const auto *fl_355 = buffer.data(fl + 355);
    const auto *fl_356 = buffer.data(fl + 356);
    const auto *fl_357 = buffer.data(fl + 357);
    const auto *fl_358 = buffer.data(fl + 358);
    const auto *fl_359 = buffer.data(fl + 359);
    const auto *fl_363 = buffer.data(fl + 363);
    const auto *fl_366 = buffer.data(fl + 366);
    const auto *fl_370 = buffer.data(fl + 370);
    const auto *fl_372 = buffer.data(fl + 372);
    const auto *fl_375 = buffer.data(fl + 375);
    const auto *fl_377 = buffer.data(fl + 377);
    const auto *fl_378 = buffer.data(fl + 378);
    const auto *fl_381 = buffer.data(fl + 381);
    const auto *fl_383 = buffer.data(fl + 383);
    const auto *fl_384 = buffer.data(fl + 384);
    const auto *fl_385 = buffer.data(fl + 385);
    const auto *fl_396 = buffer.data(fl + 396);
    const auto *fl_397 = buffer.data(fl + 397);
    const auto *fl_398 = buffer.data(fl + 398);
    const auto *fl_399 = buffer.data(fl + 399);
    const auto *fl_400 = buffer.data(fl + 400);
    const auto *fl_401 = buffer.data(fl + 401);
    const auto *fl_402 = buffer.data(fl + 402);
    const auto *fl_403 = buffer.data(fl + 403);
    const auto *fl_404 = buffer.data(fl + 404);
    const auto *fl_405 = buffer.data(fl + 405);
    const auto *fl_407 = buffer.data(fl + 407);
    const auto *fl_408 = buffer.data(fl + 408);
    const auto *fl_410 = buffer.data(fl + 410);
    const auto *fl_411 = buffer.data(fl + 411);
    const auto *fl_414 = buffer.data(fl + 414);
    const auto *fl_415 = buffer.data(fl + 415);
    const auto *fl_417 = buffer.data(fl + 417);
    const auto *fl_419 = buffer.data(fl + 419);
    const auto *fl_420 = buffer.data(fl + 420);
    const auto *fl_422 = buffer.data(fl + 422);
    const auto *fl_423 = buffer.data(fl + 423);
    const auto *fl_425 = buffer.data(fl + 425);
    const auto *fl_426 = buffer.data(fl + 426);
    const auto *fl_428 = buffer.data(fl + 428);
    const auto *fl_429 = buffer.data(fl + 429);
    const auto *fl_430 = buffer.data(fl + 430);
    const auto *fl_432 = buffer.data(fl + 432);
    const auto *fl_441 = buffer.data(fl + 441);
    const auto *fl_442 = buffer.data(fl + 442);
    const auto *fl_443 = buffer.data(fl + 443);
    const auto *fl_444 = buffer.data(fl + 444);
    const auto *fl_445 = buffer.data(fl + 445);
    const auto *fl_446 = buffer.data(fl + 446);
    const auto *fl_447 = buffer.data(fl + 447);
    const auto *fl_449 = buffer.data(fl + 449);

    const auto *gi0_0 = buffer.data(gi0 + 0);
    const auto *gi0_1 = buffer.data(gi0 + 1);
    const auto *gi0_2 = buffer.data(gi0 + 2);
    const auto *gi0_3 = buffer.data(gi0 + 3);
    const auto *gi0_5 = buffer.data(gi0 + 5);
    const auto *gi0_6 = buffer.data(gi0 + 6);
    const auto *gi0_8 = buffer.data(gi0 + 8);
    const auto *gi0_9 = buffer.data(gi0 + 9);
    const auto *gi0_10 = buffer.data(gi0 + 10);
    const auto *gi0_12 = buffer.data(gi0 + 12);
    const auto *gi0_13 = buffer.data(gi0 + 13);
    const auto *gi0_14 = buffer.data(gi0 + 14);
    const auto *gi0_21 = buffer.data(gi0 + 21);
    const auto *gi0_23 = buffer.data(gi0 + 23);
    const auto *gi0_24 = buffer.data(gi0 + 24);
    const auto *gi0_25 = buffer.data(gi0 + 25);
    const auto *gi0_26 = buffer.data(gi0 + 26);
    const auto *gi0_27 = buffer.data(gi0 + 27);
    const auto *gi0_84 = buffer.data(gi0 + 84);
    const auto *gi0_86 = buffer.data(gi0 + 86);
    const auto *gi0_87 = buffer.data(gi0 + 87);
    const auto *gi0_89 = buffer.data(gi0 + 89);
    const auto *gi0_90 = buffer.data(gi0 + 90);
    const auto *gi0_91 = buffer.data(gi0 + 91);
    const auto *gi0_93 = buffer.data(gi0 + 93);
    const auto *gi0_94 = buffer.data(gi0 + 94);
    const auto *gi0_95 = buffer.data(gi0 + 95);
    const auto *gi0_96 = buffer.data(gi0 + 96);
    const auto *gi0_98 = buffer.data(gi0 + 98);
    const auto *gi0_99 = buffer.data(gi0 + 99);
    const auto *gi0_105 = buffer.data(gi0 + 105);
    const auto *gi0_106 = buffer.data(gi0 + 106);
    const auto *gi0_107 = buffer.data(gi0 + 107);
    const auto *gi0_108 = buffer.data(gi0 + 108);
    const auto *gi0_109 = buffer.data(gi0 + 109);
    const auto *gi0_111 = buffer.data(gi0 + 111);
    const auto *gi0_140 = buffer.data(gi0 + 140);
    const auto *gi0_141 = buffer.data(gi0 + 141);
    const auto *gi0_143 = buffer.data(gi0 + 143);
    const auto *gi0_145 = buffer.data(gi0 + 145);
    const auto *gi0_146 = buffer.data(gi0 + 146);
    const auto *gi0_148 = buffer.data(gi0 + 148);
    const auto *gi0_149 = buffer.data(gi0 + 149);
    const auto *gi0_150 = buffer.data(gi0 + 150);
    const auto *gi0_152 = buffer.data(gi0 + 152);
    const auto *gi0_153 = buffer.data(gi0 + 153);
    const auto *gi0_154 = buffer.data(gi0 + 154);
    const auto *gi0_160 = buffer.data(gi0 + 160);
    const auto *gi0_161 = buffer.data(gi0 + 161);
    const auto *gi0_163 = buffer.data(gi0 + 163);
    const auto *gi0_164 = buffer.data(gi0 + 164);
    const auto *gi0_165 = buffer.data(gi0 + 165);
    const auto *gi0_166 = buffer.data(gi0 + 166);
    const auto *gi0_167 = buffer.data(gi0 + 167);
    const auto *gi0_280 = buffer.data(gi0 + 280);
    const auto *gi0_283 = buffer.data(gi0 + 283);
    const auto *gi0_285 = buffer.data(gi0 + 285);
    const auto *gi0_286 = buffer.data(gi0 + 286);
    const auto *gi0_289 = buffer.data(gi0 + 289);
    const auto *gi0_290 = buffer.data(gi0 + 290);
    const auto *gi0_292 = buffer.data(gi0 + 292);
    const auto *gi0_294 = buffer.data(gi0 + 294);
    const auto *gi0_295 = buffer.data(gi0 + 295);
    const auto *gi0_297 = buffer.data(gi0 + 297);
    const auto *gi0_298 = buffer.data(gi0 + 298);
    const auto *gi0_300 = buffer.data(gi0 + 300);
    const auto *gi0_301 = buffer.data(gi0 + 301);
    const auto *gi0_302 = buffer.data(gi0 + 302);
    const auto *gi0_303 = buffer.data(gi0 + 303);
    const auto *gi0_304 = buffer.data(gi0 + 304);
    const auto *gi0_305 = buffer.data(gi0 + 305);
    const auto *gi0_307 = buffer.data(gi0 + 307);
    const auto *gi0_336 = buffer.data(gi0 + 336);
    const auto *gi0_339 = buffer.data(gi0 + 339);
    const auto *gi0_341 = buffer.data(gi0 + 341);
    const auto *gi0_342 = buffer.data(gi0 + 342);
    const auto *gi0_345 = buffer.data(gi0 + 345);
    const auto *gi0_346 = buffer.data(gi0 + 346);
    const auto *gi0_348 = buffer.data(gi0 + 348);
    const auto *gi0_350 = buffer.data(gi0 + 350);
    const auto *gi0_351 = buffer.data(gi0 + 351);
    const auto *gi0_353 = buffer.data(gi0 + 353);
    const auto *gi0_354 = buffer.data(gi0 + 354);
    const auto *gi0_356 = buffer.data(gi0 + 356);
    const auto *gi0_357 = buffer.data(gi0 + 357);
    const auto *gi0_359 = buffer.data(gi0 + 359);
    const auto *gi0_360 = buffer.data(gi0 + 360);
    const auto *gi0_361 = buffer.data(gi0 + 361);
    const auto *gi0_362 = buffer.data(gi0 + 362);
    const auto *gi0_363 = buffer.data(gi0 + 363);
    const auto *gi0_392 = buffer.data(gi0 + 392);
    const auto *gi0_395 = buffer.data(gi0 + 395);
    const auto *gi0_397 = buffer.data(gi0 + 397);
    const auto *gi0_398 = buffer.data(gi0 + 398);
    const auto *gi0_401 = buffer.data(gi0 + 401);
    const auto *gi0_402 = buffer.data(gi0 + 402);
    const auto *gi0_404 = buffer.data(gi0 + 404);
    const auto *gi0_406 = buffer.data(gi0 + 406);
    const auto *gi0_407 = buffer.data(gi0 + 407);
    const auto *gi0_409 = buffer.data(gi0 + 409);
    const auto *gi0_410 = buffer.data(gi0 + 410);
    const auto *gi0_412 = buffer.data(gi0 + 412);
    const auto *gi0_413 = buffer.data(gi0 + 413);
    const auto *gi0_415 = buffer.data(gi0 + 415);
    const auto *gi0_416 = buffer.data(gi0 + 416);
    const auto *gi0_417 = buffer.data(gi0 + 417);
    const auto *gi0_418 = buffer.data(gi0 + 418);
    const auto *gi0_419 = buffer.data(gi0 + 419);

    const auto *gi1_0 = buffer.data(gi1 + 0);
    const auto *gi1_1 = buffer.data(gi1 + 1);
    const auto *gi1_2 = buffer.data(gi1 + 2);
    const auto *gi1_3 = buffer.data(gi1 + 3);
    const auto *gi1_5 = buffer.data(gi1 + 5);
    const auto *gi1_6 = buffer.data(gi1 + 6);
    const auto *gi1_8 = buffer.data(gi1 + 8);
    const auto *gi1_9 = buffer.data(gi1 + 9);
    const auto *gi1_10 = buffer.data(gi1 + 10);
    const auto *gi1_12 = buffer.data(gi1 + 12);
    const auto *gi1_13 = buffer.data(gi1 + 13);
    const auto *gi1_14 = buffer.data(gi1 + 14);
    const auto *gi1_21 = buffer.data(gi1 + 21);
    const auto *gi1_23 = buffer.data(gi1 + 23);
    const auto *gi1_24 = buffer.data(gi1 + 24);
    const auto *gi1_25 = buffer.data(gi1 + 25);
    const auto *gi1_26 = buffer.data(gi1 + 26);
    const auto *gi1_27 = buffer.data(gi1 + 27);
    const auto *gi1_84 = buffer.data(gi1 + 84);
    const auto *gi1_86 = buffer.data(gi1 + 86);
    const auto *gi1_87 = buffer.data(gi1 + 87);
    const auto *gi1_89 = buffer.data(gi1 + 89);
    const auto *gi1_90 = buffer.data(gi1 + 90);
    const auto *gi1_91 = buffer.data(gi1 + 91);
    const auto *gi1_93 = buffer.data(gi1 + 93);
    const auto *gi1_94 = buffer.data(gi1 + 94);
    const auto *gi1_95 = buffer.data(gi1 + 95);
    const auto *gi1_96 = buffer.data(gi1 + 96);
    const auto *gi1_98 = buffer.data(gi1 + 98);
    const auto *gi1_99 = buffer.data(gi1 + 99);
    const auto *gi1_105 = buffer.data(gi1 + 105);
    const auto *gi1_106 = buffer.data(gi1 + 106);
    const auto *gi1_107 = buffer.data(gi1 + 107);
    const auto *gi1_108 = buffer.data(gi1 + 108);
    const auto *gi1_109 = buffer.data(gi1 + 109);
    const auto *gi1_111 = buffer.data(gi1 + 111);
    const auto *gi1_140 = buffer.data(gi1 + 140);
    const auto *gi1_141 = buffer.data(gi1 + 141);
    const auto *gi1_143 = buffer.data(gi1 + 143);
    const auto *gi1_145 = buffer.data(gi1 + 145);
    const auto *gi1_146 = buffer.data(gi1 + 146);
    const auto *gi1_148 = buffer.data(gi1 + 148);
    const auto *gi1_149 = buffer.data(gi1 + 149);
    const auto *gi1_150 = buffer.data(gi1 + 150);
    const auto *gi1_152 = buffer.data(gi1 + 152);
    const auto *gi1_153 = buffer.data(gi1 + 153);
    const auto *gi1_154 = buffer.data(gi1 + 154);
    const auto *gi1_160 = buffer.data(gi1 + 160);
    const auto *gi1_161 = buffer.data(gi1 + 161);
    const auto *gi1_163 = buffer.data(gi1 + 163);
    const auto *gi1_164 = buffer.data(gi1 + 164);
    const auto *gi1_165 = buffer.data(gi1 + 165);
    const auto *gi1_166 = buffer.data(gi1 + 166);
    const auto *gi1_167 = buffer.data(gi1 + 167);
    const auto *gi1_280 = buffer.data(gi1 + 280);
    const auto *gi1_283 = buffer.data(gi1 + 283);
    const auto *gi1_285 = buffer.data(gi1 + 285);
    const auto *gi1_286 = buffer.data(gi1 + 286);
    const auto *gi1_289 = buffer.data(gi1 + 289);
    const auto *gi1_290 = buffer.data(gi1 + 290);
    const auto *gi1_292 = buffer.data(gi1 + 292);
    const auto *gi1_294 = buffer.data(gi1 + 294);
    const auto *gi1_295 = buffer.data(gi1 + 295);
    const auto *gi1_297 = buffer.data(gi1 + 297);
    const auto *gi1_298 = buffer.data(gi1 + 298);
    const auto *gi1_300 = buffer.data(gi1 + 300);
    const auto *gi1_301 = buffer.data(gi1 + 301);
    const auto *gi1_302 = buffer.data(gi1 + 302);
    const auto *gi1_303 = buffer.data(gi1 + 303);
    const auto *gi1_304 = buffer.data(gi1 + 304);
    const auto *gi1_305 = buffer.data(gi1 + 305);
    const auto *gi1_307 = buffer.data(gi1 + 307);
    const auto *gi1_336 = buffer.data(gi1 + 336);
    const auto *gi1_339 = buffer.data(gi1 + 339);
    const auto *gi1_341 = buffer.data(gi1 + 341);
    const auto *gi1_342 = buffer.data(gi1 + 342);
    const auto *gi1_345 = buffer.data(gi1 + 345);
    const auto *gi1_346 = buffer.data(gi1 + 346);
    const auto *gi1_348 = buffer.data(gi1 + 348);
    const auto *gi1_350 = buffer.data(gi1 + 350);
    const auto *gi1_351 = buffer.data(gi1 + 351);
    const auto *gi1_353 = buffer.data(gi1 + 353);
    const auto *gi1_354 = buffer.data(gi1 + 354);
    const auto *gi1_356 = buffer.data(gi1 + 356);
    const auto *gi1_357 = buffer.data(gi1 + 357);
    const auto *gi1_359 = buffer.data(gi1 + 359);
    const auto *gi1_360 = buffer.data(gi1 + 360);
    const auto *gi1_361 = buffer.data(gi1 + 361);
    const auto *gi1_362 = buffer.data(gi1 + 362);
    const auto *gi1_363 = buffer.data(gi1 + 363);
    const auto *gi1_392 = buffer.data(gi1 + 392);
    const auto *gi1_395 = buffer.data(gi1 + 395);
    const auto *gi1_397 = buffer.data(gi1 + 397);
    const auto *gi1_398 = buffer.data(gi1 + 398);
    const auto *gi1_401 = buffer.data(gi1 + 401);
    const auto *gi1_402 = buffer.data(gi1 + 402);
    const auto *gi1_404 = buffer.data(gi1 + 404);
    const auto *gi1_406 = buffer.data(gi1 + 406);
    const auto *gi1_407 = buffer.data(gi1 + 407);
    const auto *gi1_409 = buffer.data(gi1 + 409);
    const auto *gi1_410 = buffer.data(gi1 + 410);
    const auto *gi1_412 = buffer.data(gi1 + 412);
    const auto *gi1_413 = buffer.data(gi1 + 413);
    const auto *gi1_415 = buffer.data(gi1 + 415);
    const auto *gi1_416 = buffer.data(gi1 + 416);
    const auto *gi1_417 = buffer.data(gi1 + 417);
    const auto *gi1_418 = buffer.data(gi1 + 418);
    const auto *gi1_419 = buffer.data(gi1 + 419);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_1 = buffer.data(gk + 1);
    const auto *gk_2 = buffer.data(gk + 2);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_99 = buffer.data(gk + 99);
    const auto *gk_100 = buffer.data(gk + 100);
    const auto *gk_101 = buffer.data(gk + 101);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_110 = buffer.data(gk + 110);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_115 = buffer.data(gk + 115);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_119 = buffer.data(gk + 119);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_124 = buffer.data(gk + 124);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_126 = buffer.data(gk + 126);
    const auto *gk_128 = buffer.data(gk + 128);
    const auto *gk_129 = buffer.data(gk + 129);
    const auto *gk_136 = buffer.data(gk + 136);
    const auto *gk_137 = buffer.data(gk + 137);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_149 = buffer.data(gk + 149);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_164 = buffer.data(gk + 164);
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
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_188 = buffer.data(gk + 188);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_193 = buffer.data(gk + 193);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_199 = buffer.data(gk + 199);
    const auto *gk_200 = buffer.data(gk + 200);
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
    const auto *gk_219 = buffer.data(gk + 219);
    const auto *gk_221 = buffer.data(gk + 221);
    const auto *gk_222 = buffer.data(gk + 222);
    const auto *gk_225 = buffer.data(gk + 225);
    const auto *gk_226 = buffer.data(gk + 226);
    const auto *gk_230 = buffer.data(gk + 230);
    const auto *gk_231 = buffer.data(gk + 231);
    const auto *gk_236 = buffer.data(gk + 236);
    const auto *gk_237 = buffer.data(gk + 237);
    const auto *gk_244 = buffer.data(gk + 244);
    const auto *gk_246 = buffer.data(gk + 246);
    const auto *gk_247 = buffer.data(gk + 247);
    const auto *gk_248 = buffer.data(gk + 248);
    const auto *gk_249 = buffer.data(gk + 249);
    const auto *gk_250 = buffer.data(gk + 250);
    const auto *gk_251 = buffer.data(gk + 251);
    const auto *gk_252 = buffer.data(gk + 252);
    const auto *gk_254 = buffer.data(gk + 254);
    const auto *gk_255 = buffer.data(gk + 255);
    const auto *gk_257 = buffer.data(gk + 257);
    const auto *gk_258 = buffer.data(gk + 258);
    const auto *gk_261 = buffer.data(gk + 261);
    const auto *gk_262 = buffer.data(gk + 262);
    const auto *gk_266 = buffer.data(gk + 266);
    const auto *gk_267 = buffer.data(gk + 267);
    const auto *gk_272 = buffer.data(gk + 272);
    const auto *gk_281 = buffer.data(gk + 281);
    const auto *gk_282 = buffer.data(gk + 282);
    const auto *gk_283 = buffer.data(gk + 283);
    const auto *gk_284 = buffer.data(gk + 284);
    const auto *gk_285 = buffer.data(gk + 285);
    const auto *gk_286 = buffer.data(gk + 286);
    const auto *gk_287 = buffer.data(gk + 287);
    const auto *gk_288 = buffer.data(gk + 288);
    const auto *gk_290 = buffer.data(gk + 290);
    const auto *gk_291 = buffer.data(gk + 291);
    const auto *gk_293 = buffer.data(gk + 293);
    const auto *gk_294 = buffer.data(gk + 294);
    const auto *gk_297 = buffer.data(gk + 297);
    const auto *gk_298 = buffer.data(gk + 298);
    const auto *gk_302 = buffer.data(gk + 302);
    const auto *gk_303 = buffer.data(gk + 303);
    const auto *gk_308 = buffer.data(gk + 308);
    const auto *gk_316 = buffer.data(gk + 316);
    const auto *gk_317 = buffer.data(gk + 317);
    const auto *gk_318 = buffer.data(gk + 318);
    const auto *gk_319 = buffer.data(gk + 319);
    const auto *gk_320 = buffer.data(gk + 320);
    const auto *gk_321 = buffer.data(gk + 321);
    const auto *gk_322 = buffer.data(gk + 322);
    const auto *gk_324 = buffer.data(gk + 324);
    const auto *gk_326 = buffer.data(gk + 326);
    const auto *gk_327 = buffer.data(gk + 327);
    const auto *gk_329 = buffer.data(gk + 329);
    const auto *gk_330 = buffer.data(gk + 330);
    const auto *gk_333 = buffer.data(gk + 333);
    const auto *gk_334 = buffer.data(gk + 334);
    const auto *gk_338 = buffer.data(gk + 338);
    const auto *gk_339 = buffer.data(gk + 339);
    const auto *gk_344 = buffer.data(gk + 344);
    const auto *gk_351 = buffer.data(gk + 351);
    const auto *gk_352 = buffer.data(gk + 352);
    const auto *gk_353 = buffer.data(gk + 353);
    const auto *gk_354 = buffer.data(gk + 354);
    const auto *gk_355 = buffer.data(gk + 355);
    const auto *gk_356 = buffer.data(gk + 356);
    const auto *gk_357 = buffer.data(gk + 357);
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
    const auto *gk_381 = buffer.data(gk + 381);
    const auto *gk_383 = buffer.data(gk + 383);
    const auto *gk_384 = buffer.data(gk + 384);
    const auto *gk_385 = buffer.data(gk + 385);
    const auto *gk_387 = buffer.data(gk + 387);
    const auto *gk_388 = buffer.data(gk + 388);
    const auto *gk_389 = buffer.data(gk + 389);
    const auto *gk_390 = buffer.data(gk + 390);
    const auto *gk_391 = buffer.data(gk + 391);
    const auto *gk_392 = buffer.data(gk + 392);
    const auto *gk_393 = buffer.data(gk + 393);
    const auto *gk_394 = buffer.data(gk + 394);
    const auto *gk_395 = buffer.data(gk + 395);
    const auto *gk_396 = buffer.data(gk + 396);
    const auto *gk_398 = buffer.data(gk + 398);
    const auto *gk_399 = buffer.data(gk + 399);
    const auto *gk_401 = buffer.data(gk + 401);
    const auto *gk_402 = buffer.data(gk + 402);
    const auto *gk_405 = buffer.data(gk + 405);
    const auto *gk_406 = buffer.data(gk + 406);
    const auto *gk_410 = buffer.data(gk + 410);
    const auto *gk_411 = buffer.data(gk + 411);
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
    const auto *gk_434 = buffer.data(gk + 434);
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
    const auto *gk_453 = buffer.data(gk + 453);
    const auto *gk_455 = buffer.data(gk + 455);
    const auto *gk_456 = buffer.data(gk + 456);
    const auto *gk_457 = buffer.data(gk + 457);
    const auto *gk_459 = buffer.data(gk + 459);
    const auto *gk_460 = buffer.data(gk + 460);
    const auto *gk_461 = buffer.data(gk + 461);
    const auto *gk_462 = buffer.data(gk + 462);
    const auto *gk_463 = buffer.data(gk + 463);
    const auto *gk_464 = buffer.data(gk + 464);
    const auto *gk_465 = buffer.data(gk + 465);
    const auto *gk_466 = buffer.data(gk + 466);
    const auto *gk_467 = buffer.data(gk + 467);
    const auto *gk_468 = buffer.data(gk + 468);
    const auto *gk_470 = buffer.data(gk + 470);
    const auto *gk_471 = buffer.data(gk + 471);
    const auto *gk_473 = buffer.data(gk + 473);
    const auto *gk_474 = buffer.data(gk + 474);
    const auto *gk_477 = buffer.data(gk + 477);
    const auto *gk_478 = buffer.data(gk + 478);
    const auto *gk_482 = buffer.data(gk + 482);
    const auto *gk_483 = buffer.data(gk + 483);
    const auto *gk_488 = buffer.data(gk + 488);
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
    const auto *gk_525 = buffer.data(gk + 525);
    const auto *gk_527 = buffer.data(gk + 527);
    const auto *gk_528 = buffer.data(gk + 528);
    const auto *gk_529 = buffer.data(gk + 529);
    const auto *gk_531 = buffer.data(gk + 531);
    const auto *gk_532 = buffer.data(gk + 532);
    const auto *gk_533 = buffer.data(gk + 533);
    const auto *gk_534 = buffer.data(gk + 534);
    const auto *gk_535 = buffer.data(gk + 535);
    const auto *gk_536 = buffer.data(gk + 536);
    const auto *gk_537 = buffer.data(gk + 537);
    const auto *gk_538 = buffer.data(gk + 538);
    const auto *gk_539 = buffer.data(gk + 539);

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
                         gi1_2, gi1_3, gk_3, gk_5, gk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * gi0_1[k]
                 - f_6 * gi1_1[k]
                 + pb_y[k] * gk_3[k];

        t_7[k] = pb_z[k] * gk_3[k];

        t_8[k] = pb_y[k] * gk_5[k];

        t_9[k] = f_5 * gi0_2[k]
                 - f_6 * gi1_2[k]
                 + pb_z[k] * gk_5[k];

        t_10[k] = f_7 * gi0_3[k]
                  - f_8 * gi1_3[k]
                  + pb_y[k] * gk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, gi0_5, gi0_6, gi1_5, \
                         gi1_6, gk_6, gk_8, gk_9, gk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gk_6[k];

        t_12[k] = f_3 * gi0_5[k]
                  - f_4 * gi1_5[k]
                  + pb_y[k] * gk_8[k];

        t_13[k] = pb_y[k] * gk_9[k];

        t_14[k] = f_7 * gi0_5[k]
                  - f_8 * gi1_5[k]
                  + pb_z[k] * gk_9[k];

        t_15[k] = f_9 * gi0_6[k]
                  - f_10 * gi1_6[k]
                  + pb_y[k] * gk_10[k];

        t_16[k] = pb_z[k] * gk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, gi0_8, gi0_9, gi1_8, gi1_9, \
                         gk_12, gk_13, gk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * gi0_8[k]
                  - f_6 * gi1_8[k]
                  + pb_y[k] * gk_12[k];

        t_18[k] = f_3 * gi0_9[k]
                  - f_4 * gi1_9[k]
                  + pb_y[k] * gk_13[k];

        t_19[k] = pb_y[k] * gk_14[k];

        t_20[k] = f_9 * gi0_9[k]
                  - f_10 * gi1_9[k]
                  + pb_z[k] * gk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, gi0_10, gi0_12, gi0_13, gi1_10, \
                         gi1_12, gi1_13, gk_15, gk_17, gk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * gi0_10[k]
                  - f_12 * gi1_10[k]
                  + pb_y[k] * gk_15[k];

        t_22[k] = pb_z[k] * gk_15[k];

        t_23[k] = f_7 * gi0_12[k]
                  - f_8 * gi1_12[k]
                  + pb_y[k] * gk_17[k];

        t_24[k] = f_5 * gi0_13[k]
                  - f_6 * gi1_13[k]
                  + pb_y[k] * gk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, fk_28, gi0_14, \
                         gi1_14, gk_19, gk_20, gk_21, gk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gi0_14[k]
                  - f_4 * gi1_14[k]
                  + pb_y[k] * gk_19[k];

        t_26[k] = pb_y[k] * gk_20[k];

        t_27[k] = f_11 * gi0_14[k]
                  - f_12 * gi1_14[k]
                  + pb_z[k] * gk_20[k];

        t_28[k] = f_0 * fk_28[k]
                  + pb_x[k] * gk_28[k];

        t_29[k] = pb_z[k] * gk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, fk_30, fk_31, fk_32, fk_33, \
                         gk_27, gk_30, gk_31, gk_32, gk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * fk_30[k]
                  + pb_x[k] * gk_30[k];

        t_31[k] = f_0 * fk_31[k]
                  + pb_x[k] * gk_31[k];

        t_32[k] = f_0 * fk_32[k]
                  + pb_x[k] * gk_32[k];

        t_33[k] = f_0 * fk_33[k]
                  + pb_x[k] * gk_33[k];

        t_34[k] = pb_y[k] * gk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, fk_35, gi0_21, gi0_23, \
                         gi1_21, gi1_23, gk_28, gk_30, gk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * fk_35[k]
                  + pb_x[k] * gk_35[k];

        t_36[k] = f_1 * gi0_21[k]
                  - f_2 * gi1_21[k]
                  + pb_y[k] * gk_28[k];

        t_37[k] = pb_z[k] * gk_28[k];

        t_38[k] = f_11 * gi0_23[k]
                  - f_12 * gi1_23[k]
                  + pb_y[k] * gk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, gi0_24, gi0_25, gi0_26, gi1_24, gi1_25, \
                         gi1_26, gk_31, gk_32, gk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * gi0_24[k]
                  - f_10 * gi1_24[k]
                  + pb_y[k] * gk_31[k];

        t_40[k] = f_7 * gi0_25[k]
                  - f_8 * gi1_25[k]
                  + pb_y[k] * gk_32[k];

        t_41[k] = f_5 * gi0_26[k]
                  - f_6 * gi1_26[k]
                  + pb_y[k] * gk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, fk_0, fl_0, \
                         gi0_27, gi1_27, gk_34, gk_35, gk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * gi0_27[k]
                  - f_4 * gi1_27[k]
                  + pb_y[k] * gk_34[k];

        t_43[k] = pb_y[k] * gk_35[k];

        t_44[k] = f_1 * gi0_27[k]
                  - f_2 * gi1_27[k]
                  + pb_z[k] * gk_35[k];

        t_45[k] = pa_y[k] * fl_0[k];

        t_46[k] = f_13 * fk_0[k]
                  + pb_y[k] * gk_36[k];

        t_47[k] = pb_z[k] * gk_36[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_z, fk_1, fk_3, fl_3, fl_5, \
                         fl_6, gk_37, gk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * fk_1[k]
                  + pa_y[k] * fl_3[k];

        t_49[k] = pb_z[k] * gk_37[k];

        t_50[k] = pa_y[k] * fl_5[k];

        t_51[k] = f_15 * fk_3[k]
                  + pa_y[k] * fl_6[k];

        t_52[k] = pb_z[k] * gk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_y, pb_z, fk_5, fk_6, fk_8, \
                         fl_9, fl_10, fl_12, gk_41, gk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * fk_5[k]
                  + pb_y[k] * gk_41[k];

        t_54[k] = pa_y[k] * fl_9[k];

        t_55[k] = f_0 * fk_6[k]
                  + pa_y[k] * fl_10[k];

        t_56[k] = pb_z[k] * gk_42[k];

        t_57[k] = f_14 * fk_8[k]
                  + pa_y[k] * fl_12[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, fk_9, fk_10, fk_12, \
                         fl_14, fl_15, fl_17, gk_45, gk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_13 * fk_9[k]
                  + pb_y[k] * gk_45[k];

        t_59[k] = pa_y[k] * fl_14[k];

        t_60[k] = f_16 * fk_10[k]
                  + pa_y[k] * fl_15[k];

        t_61[k] = pb_z[k] * gk_46[k];

        t_62[k] = f_15 * fk_12[k]
                  + pa_y[k] * fl_17[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, fk_13, fk_14, fk_15, \
                         fl_18, fl_20, fl_21, gk_50, gk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * fk_13[k]
                  + pa_y[k] * fl_18[k];

        t_64[k] = f_13 * fk_14[k]
                  + pb_y[k] * gk_50[k];

        t_65[k] = pa_y[k] * fl_20[k];

        t_66[k] = f_17 * fk_15[k]
                  + pa_y[k] * fl_21[k];

        t_67[k] = pb_z[k] * gk_51[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pb_y, fk_17, fk_18, fk_19, fk_20, \
                         fl_23, fl_24, fl_25, fl_27, gk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * fk_17[k]
                  + pa_y[k] * fl_23[k];

        t_69[k] = f_15 * fk_18[k]
                  + pa_y[k] * fl_24[k];

        t_70[k] = f_14 * fk_19[k]
                  + pa_y[k] * fl_25[k];

        t_71[k] = f_13 * fk_20[k]
                  + pb_y[k] * gk_56[k];

        t_72[k] = pa_y[k] * fl_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, fk_64, fk_66, fk_67, fk_68, \
                         gk_57, gk_64, gk_66, gk_67, gk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * fk_64[k]
                  + pb_x[k] * gk_64[k];

        t_74[k] = pb_z[k] * gk_57[k];

        t_75[k] = f_15 * fk_66[k]
                  + pb_x[k] * gk_66[k];

        t_76[k] = f_15 * fk_67[k]
                  + pb_x[k] * gk_67[k];

        t_77[k] = f_15 * fk_68[k]
                  + pb_x[k] * gk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, fk_28, fk_69, fk_70, \
                         fl_35, fl_36, gk_64, gk_69, gk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_15 * fk_69[k]
                  + pb_x[k] * gk_69[k];

        t_79[k] = f_15 * fk_70[k]
                  + pb_x[k] * gk_70[k];

        t_80[k] = pa_y[k] * fl_35[k];

        t_81[k] = f_18 * fk_28[k]
                  + pa_y[k] * fl_36[k];

        t_82[k] = pb_z[k] * gk_64[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, fk_30, fk_31, fk_32, fk_33, \
                         fk_34, fl_38, fl_39, fl_40, fl_41, fl_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_17 * fk_30[k]
                  + pa_y[k] * fl_38[k];

        t_84[k] = f_16 * fk_31[k]
                  + pa_y[k] * fl_39[k];

        t_85[k] = f_0 * fk_32[k]
                  + pa_y[k] * fl_40[k];

        t_86[k] = f_15 * fk_33[k]
                  + pa_y[k] * fl_41[k];

        t_87[k] = f_14 * fk_34[k]
                  + pa_y[k] * fl_42[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, fk_0, fk_35, \
                         fl_0, fl_44, gk_71, gk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_13 * fk_35[k]
                  + pb_y[k] * gk_71[k];

        t_89[k] = pa_y[k] * fl_44[k];

        t_90[k] = pa_z[k] * fl_0[k];

        t_91[k] = pb_y[k] * gk_72[k];

        t_92[k] = f_13 * fk_0[k]
                  + pb_z[k] * gk_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, fk_2, fk_3, fl_3, \
                         fl_5, fl_6, gk_74, gk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * fl_3[k];

        t_94[k] = pb_y[k] * gk_74[k];

        t_95[k] = f_14 * fk_2[k]
                  + pa_z[k] * fl_5[k];

        t_96[k] = pa_z[k] * fl_6[k];

        t_97[k] = f_13 * fk_3[k]
                  + pb_z[k] * gk_75[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, fk_5, fk_6, fk_7, \
                         fl_9, fl_10, fl_12, gk_77, gk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * gk_77[k];

        t_99[k] = f_15 * fk_5[k]
                  + pa_z[k] * fl_9[k];

        t_100[k] = pa_z[k] * fl_10[k];

        t_101[k] = f_13 * fk_6[k]
                   + pb_z[k] * gk_78[k];

        t_102[k] = f_14 * fk_7[k]
                   + pa_z[k] * fl_12[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, fk_9, fk_10, \
                         fk_11, fl_14, fl_15, fl_17, gk_81, gk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * gk_81[k];

        t_104[k] = f_0 * fk_9[k]
                   + pa_z[k] * fl_14[k];

        t_105[k] = pa_z[k] * fl_15[k];

        t_106[k] = f_13 * fk_10[k]
                   + pb_z[k] * gk_82[k];

        t_107[k] = f_14 * fk_11[k]
                   + pa_z[k] * fl_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, fk_12, fk_14, \
                         fk_15, fl_18, fl_20, fl_21, gk_86, gk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * fk_12[k]
                   + pa_z[k] * fl_18[k];

        t_109[k] = pb_y[k] * gk_86[k];

        t_110[k] = f_16 * fk_14[k]
                   + pa_z[k] * fl_20[k];

        t_111[k] = pa_z[k] * fl_21[k];

        t_112[k] = f_13 * fk_15[k]
                   + pb_z[k] * gk_87[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, fk_16, fk_17, fk_18, \
                         fk_20, fl_23, fl_24, fl_25, fl_27, gk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * fk_16[k]
                   + pa_z[k] * fl_23[k];

        t_114[k] = f_15 * fk_17[k]
                   + pa_z[k] * fl_24[k];

        t_115[k] = f_0 * fk_18[k]
                   + pa_z[k] * fl_25[k];

        t_116[k] = pb_y[k] * gk_92[k];

        t_117[k] = f_17 * fk_20[k]
                   + pa_z[k] * fl_27[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, fk_101, fk_102, \
                         fk_103, fk_104, fl_28, gk_101, gk_102, gk_103, \
                         gk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * fl_28[k];

        t_119[k] = f_15 * fk_101[k]
                   + pb_x[k] * gk_101[k];

        t_120[k] = f_15 * fk_102[k]
                   + pb_x[k] * gk_102[k];

        t_121[k] = f_15 * fk_103[k]
                   + pb_x[k] * gk_103[k];

        t_122[k] = f_15 * fk_104[k]
                   + pb_x[k] * gk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, pb_y, fk_105, fk_107, fl_36, \
                         gk_99, gk_105, gk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_15 * fk_105[k]
                   + pb_x[k] * gk_105[k];

        t_124[k] = pb_y[k] * gk_99[k];

        t_125[k] = f_15 * fk_107[k]
                   + pb_x[k] * gk_107[k];

        t_126[k] = pa_z[k] * fl_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_z, pb_z, fk_28, fk_29, fk_30, fk_31, \
                         fl_38, fl_39, fl_40, gk_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * fk_28[k]
                   + pb_z[k] * gk_100[k];

        t_128[k] = f_14 * fk_29[k]
                   + pa_z[k] * fl_38[k];

        t_129[k] = f_15 * fk_30[k]
                   + pa_z[k] * fl_39[k];

        t_130[k] = f_0 * fk_31[k]
                   + pa_z[k] * fl_40[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_z, pb_y, fk_32, fk_33, fk_35, fl_41, \
                         fl_42, fl_44, gk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_16 * fk_32[k]
                   + pa_z[k] * fl_41[k];

        t_132[k] = f_17 * fk_33[k]
                   + pa_z[k] * fl_42[k];

        t_133[k] = pb_y[k] * gk_107[k];

        t_134[k] = f_18 * fk_35[k]
                   + pa_z[k] * fl_44[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pb_y, pb_z, dl0_0, dl1_0, fk_36, fl_45, \
                         gk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_19 * dl0_0[k]
                   - f_20 * dl1_0[k]
                   + pa_y[k] * fl_45[k];

        t_136[k] = f_14 * fk_36[k]
                   + pb_y[k] * gk_108[k];

        t_137[k] = pb_z[k] * gk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_x, pb_z, fk_111, gi0_84, gi0_87, gi1_84, \
                         gi1_87, gk_109, gk_110, gk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_14 * fk_111[k]
                   + f_11 * gi0_87[k]
                   - f_12 * gi1_87[k]
                   + pb_x[k] * gk_111[k];

        t_139[k] = pb_z[k] * gk_109[k];

        t_140[k] = f_3 * gi0_84[k]
                   - f_4 * gi1_84[k]
                   + pb_z[k] * gk_110[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_x, pb_y, pb_z, fk_41, fk_114, gi0_86, \
                         gi0_90, gi1_86, gi1_90, gk_111, gk_113, \
                         gk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_14 * fk_114[k]
                   + f_9 * gi0_90[k]
                   - f_10 * gi1_90[k]
                   + pb_x[k] * gk_114[k];

        t_142[k] = pb_z[k] * gk_111[k];

        t_143[k] = f_14 * fk_41[k]
                   + pb_y[k] * gk_113[k];

        t_144[k] = f_5 * gi0_86[k]
                   - f_6 * gi1_86[k]
                   + pb_z[k] * gk_113[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_z, fk_118, gi0_87, gi0_94, gi1_87, \
                         gi1_94, gk_114, gk_115, gk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_14 * fk_118[k]
                   + f_7 * gi0_94[k]
                   - f_8 * gi1_94[k]
                   + pb_x[k] * gk_118[k];

        t_146[k] = pb_z[k] * gk_114[k];

        t_147[k] = f_3 * gi0_87[k]
                   - f_4 * gi1_87[k]
                   + pb_z[k] * gk_115[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pb_y, pb_z, fk_45, fk_123, gi0_89, \
                         gi0_99, gi1_89, gi1_99, gk_117, gk_118, \
                         gk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * fk_45[k]
                   + pb_y[k] * gk_117[k];

        t_149[k] = f_7 * gi0_89[k]
                   - f_8 * gi1_89[k]
                   + pb_z[k] * gk_117[k];

        t_150[k] = f_14 * fk_123[k]
                   + f_5 * gi0_99[k]
                   - f_6 * gi1_99[k]
                   + pb_x[k] * gk_123[k];

        t_151[k] = pb_z[k] * gk_118[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, fk_50, gi0_90, gi0_91, \
                         gi0_93, gi1_90, gi1_91, gi1_93, gk_119, gk_120, \
                         gk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * gi0_90[k]
                   - f_4 * gi1_90[k]
                   + pb_z[k] * gk_119[k];

        t_153[k] = f_5 * gi0_91[k]
                   - f_6 * gi1_91[k]
                   + pb_z[k] * gk_120[k];

        t_154[k] = f_14 * fk_50[k]
                   + pb_y[k] * gk_122[k];

        t_155[k] = f_9 * gi0_93[k]
                   - f_10 * gi1_93[k]
                   + pb_z[k] * gk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_z, fk_129, gi0_94, gi0_105, gi1_94, \
                         gi1_105, gk_123, gk_124, gk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_14 * fk_129[k]
                   + f_3 * gi0_105[k]
                   - f_4 * gi1_105[k]
                   + pb_x[k] * gk_129[k];

        t_157[k] = pb_z[k] * gk_123[k];

        t_158[k] = f_3 * gi0_94[k]
                   - f_4 * gi1_94[k]
                   + pb_z[k] * gk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_y, pb_z, fk_56, gi0_95, gi0_96, \
                         gi0_98, gi1_95, gi1_96, gi1_98, gk_125, gk_126, \
                         gk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * gi0_95[k]
                   - f_6 * gi1_95[k]
                   + pb_z[k] * gk_125[k];

        t_160[k] = f_7 * gi0_96[k]
                   - f_8 * gi1_96[k]
                   + pb_z[k] * gk_126[k];

        t_161[k] = f_14 * fk_56[k]
                   + pb_y[k] * gk_128[k];

        t_162[k] = f_11 * gi0_98[k]
                   - f_12 * gi1_98[k]
                   + pb_z[k] * gk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_x, pb_z, fk_136, fk_138, \
                         fk_139, fk_140, gk_129, gk_136, gk_138, gk_139, \
                         gk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_14 * fk_136[k]
                   + pb_x[k] * gk_136[k];

        t_164[k] = pb_z[k] * gk_129[k];

        t_165[k] = f_14 * fk_138[k]
                   + pb_x[k] * gk_138[k];

        t_166[k] = f_14 * fk_139[k]
                   + pb_x[k] * gk_139[k];

        t_167[k] = f_14 * fk_140[k]
                   + pb_x[k] * gk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pb_x, dl0_171, dl1_171, fk_141, \
                         fk_142, fk_143, fl_171, gk_141, gk_142, \
                         gk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_14 * fk_141[k]
                   + pb_x[k] * gk_141[k];

        t_169[k] = f_14 * fk_142[k]
                   + pb_x[k] * gk_142[k];

        t_170[k] = f_14 * fk_143[k]
                   + pb_x[k] * gk_143[k];

        t_171[k] = f_19 * dl0_171[k]
                   - f_20 * dl1_171[k]
                   + pa_x[k] * fl_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, gi0_105, gi0_106, gi0_107, gi1_105, \
                         gi1_106, gi1_107, gk_136, gk_137, gk_138, \
                         gk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * gk_136[k];

        t_173[k] = f_3 * gi0_105[k]
                   - f_4 * gi1_105[k]
                   + pb_z[k] * gk_137[k];

        t_174[k] = f_5 * gi0_106[k]
                   - f_6 * gi1_106[k]
                   + pb_z[k] * gk_138[k];

        t_175[k] = f_7 * gi0_107[k]
                   - f_8 * gi1_107[k]
                   + pb_z[k] * gk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pb_z, fk_71, gi0_108, gi0_109, \
                         gi0_111, gi1_108, gi1_109, gi1_111, gk_140, gk_141, \
                         gk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * gi0_108[k]
                   - f_10 * gi1_108[k]
                   + pb_z[k] * gk_140[k];

        t_177[k] = f_11 * gi0_109[k]
                   - f_12 * gi1_109[k]
                   + pb_z[k] * gk_141[k];

        t_178[k] = f_14 * fk_71[k]
                   + pb_y[k] * gk_143[k];

        t_179[k] = f_1 * gi0_111[k]
                   - f_2 * gi1_111[k]
                   + pb_z[k] * gk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pa_y, pa_z, pb_y, fk_74, \
                         fl_46, fl_48, fl_90, fl_92, fl_95, gk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * fl_90[k];

        t_181[k] = pa_z[k] * fl_46[k];

        t_182[k] = pa_y[k] * fl_92[k];

        t_183[k] = pa_z[k] * fl_48[k];

        t_184[k] = f_13 * fk_74[k]
                   + pb_y[k] * gk_146[k];

        t_185[k] = pa_y[k] * fl_95[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, fk_39, \
                         fk_77, fl_51, fl_55, fl_99, gk_147, gk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * fl_51[k];

        t_187[k] = f_13 * fk_39[k]
                   + pb_z[k] * gk_147[k];

        t_188[k] = f_13 * fk_77[k]
                   + pb_y[k] * gk_149[k];

        t_189[k] = pa_y[k] * fl_99[k];

        t_190[k] = pa_z[k] * fl_55[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, fk_42, fk_80, fk_81, \
                         fl_102, fl_104, gk_150, gk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * fk_42[k]
                   + pb_z[k] * gk_150[k];

        t_192[k] = f_14 * fk_80[k]
                   + pa_y[k] * fl_102[k];

        t_193[k] = f_13 * fk_81[k]
                   + pb_y[k] * gk_153[k];

        t_194[k] = pa_y[k] * fl_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pa_z, pb_z, fk_46, fk_84, fk_85, \
                         fl_60, fl_107, fl_108, gk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * fl_60[k];

        t_196[k] = f_13 * fk_46[k]
                   + pb_z[k] * gk_154[k];

        t_197[k] = f_15 * fk_84[k]
                   + pa_y[k] * fl_107[k];

        t_198[k] = f_14 * fk_85[k]
                   + pa_y[k] * fl_108[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, fk_51, fk_86, \
                         fl_66, fl_110, gk_158, gk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * fk_86[k]
                   + pb_y[k] * gk_158[k];

        t_200[k] = pa_y[k] * fl_110[k];

        t_201[k] = pa_z[k] * fl_66[k];

        t_202[k] = f_13 * fk_51[k]
                   + pb_z[k] * gk_159[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pa_y, pb_y, fk_89, fk_90, fk_91, \
                         fk_92, fl_113, fl_114, fl_115, fl_117, \
                         gk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * fk_89[k]
                   + pa_y[k] * fl_113[k];

        t_204[k] = f_15 * fk_90[k]
                   + pa_y[k] * fl_114[k];

        t_205[k] = f_14 * fk_91[k]
                   + pa_y[k] * fl_115[k];

        t_206[k] = f_13 * fk_92[k]
                   + pb_y[k] * gk_164[k];

        t_207[k] = pa_y[k] * fl_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pa_z, pb_x, fk_173, fk_174, \
                         fk_175, fk_176, fl_73, gk_173, gk_174, gk_175, \
                         gk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_z[k] * fl_73[k];

        t_209[k] = f_14 * fk_173[k]
                   + pb_x[k] * gk_173[k];

        t_210[k] = f_14 * fk_174[k]
                   + pb_x[k] * gk_174[k];

        t_211[k] = f_14 * fk_175[k]
                   + pb_x[k] * gk_175[k];

        t_212[k] = f_14 * fk_176[k]
                   + pb_x[k] * gk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pa_z, pb_x, fk_177, fk_178, fl_81, \
                         fl_125, gk_177, gk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_14 * fk_177[k]
                   + pb_x[k] * gk_177[k];

        t_214[k] = f_14 * fk_178[k]
                   + pb_x[k] * gk_178[k];

        t_215[k] = pa_y[k] * fl_125[k];

        t_216[k] = pa_z[k] * fl_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_y, pb_z, fk_64, fk_102, fk_103, \
                         fk_104, fl_128, fl_129, fl_130, gk_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_13 * fk_64[k]
                   + pb_z[k] * gk_172[k];

        t_218[k] = f_17 * fk_102[k]
                   + pa_y[k] * fl_128[k];

        t_219[k] = f_16 * fk_103[k]
                   + pa_y[k] * fl_129[k];

        t_220[k] = f_0 * fk_104[k]
                   + pa_y[k] * fl_130[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pb_y, fk_105, fk_106, fk_107, \
                         fl_131, fl_132, fl_134, gk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_15 * fk_105[k]
                   + pa_y[k] * fl_131[k];

        t_222[k] = f_14 * fk_106[k]
                   + pa_y[k] * fl_132[k];

        t_223[k] = f_13 * fk_107[k]
                   + pb_y[k] * gk_179[k];

        t_224[k] = pa_y[k] * fl_134[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_z, pb_y, pb_z, dl0_0, dl1_0, fk_72, \
                         fl_90, gi0_140, gi1_140, gk_180, gk_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_19 * dl0_0[k]
                   - f_20 * dl1_0[k]
                   + pa_z[k] * fl_90[k];

        t_226[k] = pb_y[k] * gk_180[k];

        t_227[k] = f_14 * fk_72[k]
                   + pb_z[k] * gk_180[k];

        t_228[k] = f_3 * gi0_140[k]
                   - f_4 * gi1_140[k]
                   + pb_y[k] * gk_181[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, fk_75, fk_185, gi0_141, \
                         gi0_145, gi1_141, gi1_145, gk_182, gk_183, \
                         gk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * gk_182[k];

        t_230[k] = f_14 * fk_185[k]
                   + f_11 * gi0_145[k]
                   - f_12 * gi1_145[k]
                   + pb_x[k] * gk_185[k];

        t_231[k] = f_5 * gi0_141[k]
                   - f_6 * gi1_141[k]
                   + pb_y[k] * gk_183[k];

        t_232[k] = f_14 * fk_75[k]
                   + pb_z[k] * gk_183[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pb_y, pb_z, fk_78, fk_189, gi0_143, \
                         gi0_149, gi1_143, gi1_149, gk_185, gk_186, \
                         gk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * gk_185[k];

        t_234[k] = f_14 * fk_189[k]
                   + f_9 * gi0_149[k]
                   - f_10 * gi1_149[k]
                   + pb_x[k] * gk_189[k];

        t_235[k] = f_7 * gi0_143[k]
                   - f_8 * gi1_143[k]
                   + pb_y[k] * gk_186[k];

        t_236[k] = f_14 * fk_78[k]
                   + pb_z[k] * gk_186[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_x, pb_y, fk_194, gi0_145, gi0_154, gi1_145, \
                         gi1_154, gk_188, gk_189, gk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_3 * gi0_145[k]
                   - f_4 * gi1_145[k]
                   + pb_y[k] * gk_188[k];

        t_238[k] = pb_y[k] * gk_189[k];

        t_239[k] = f_14 * fk_194[k]
                   + f_7 * gi0_154[k]
                   - f_8 * gi1_154[k]
                   + pb_x[k] * gk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_y, pb_z, fk_82, gi0_146, gi0_148, \
                         gi0_149, gi1_146, gi1_148, gi1_149, gk_190, gk_192, \
                         gk_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_9 * gi0_146[k]
                   - f_10 * gi1_146[k]
                   + pb_y[k] * gk_190[k];

        t_241[k] = f_14 * fk_82[k]
                   + pb_z[k] * gk_190[k];

        t_242[k] = f_5 * gi0_148[k]
                   - f_6 * gi1_148[k]
                   + pb_y[k] * gk_192[k];

        t_243[k] = f_3 * gi0_149[k]
                   - f_4 * gi1_149[k]
                   + pb_y[k] * gk_193[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pb_y, pb_z, fk_87, fk_200, gi0_150, \
                         gi0_160, gi1_150, gi1_160, gk_194, gk_195, \
                         gk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_y[k] * gk_194[k];

        t_245[k] = f_14 * fk_200[k]
                   + f_5 * gi0_160[k]
                   - f_6 * gi1_160[k]
                   + pb_x[k] * gk_200[k];

        t_246[k] = f_11 * gi0_150[k]
                   - f_12 * gi1_150[k]
                   + pb_y[k] * gk_195[k];

        t_247[k] = f_14 * fk_87[k]
                   + pb_z[k] * gk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_y, gi0_152, gi0_153, gi0_154, gi1_152, \
                         gi1_153, gi1_154, gk_197, gk_198, gk_199, \
                         gk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * gi0_152[k]
                   - f_8 * gi1_152[k]
                   + pb_y[k] * gk_197[k];

        t_249[k] = f_5 * gi0_153[k]
                   - f_6 * gi1_153[k]
                   + pb_y[k] * gk_198[k];

        t_250[k] = f_3 * gi0_154[k]
                   - f_4 * gi1_154[k]
                   + pb_y[k] * gk_199[k];

        t_251[k] = pb_y[k] * gk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, fk_207, fk_208, fk_209, fk_210, \
                         gi0_167, gi1_167, gk_207, gk_208, gk_209, \
                         gk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * fk_207[k]
                   + f_3 * gi0_167[k]
                   - f_4 * gi1_167[k]
                   + pb_x[k] * gk_207[k];

        t_253[k] = f_14 * fk_208[k]
                   + pb_x[k] * gk_208[k];

        t_254[k] = f_14 * fk_209[k]
                   + pb_x[k] * gk_209[k];

        t_255[k] = f_14 * fk_210[k]
                   + pb_x[k] * gk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, fk_211, fk_212, \
                         fk_213, fk_215, gk_207, gk_211, gk_212, gk_213, \
                         gk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * fk_211[k]
                   + pb_x[k] * gk_211[k];

        t_257[k] = f_14 * fk_212[k]
                   + pb_x[k] * gk_212[k];

        t_258[k] = f_14 * fk_213[k]
                   + pb_x[k] * gk_213[k];

        t_259[k] = pb_y[k] * gk_207[k];

        t_260[k] = f_14 * fk_215[k]
                   + pb_x[k] * gk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_y, pb_z, fk_100, gi0_161, gi0_163, \
                         gi0_164, gi1_161, gi1_163, gi1_164, gk_208, gk_210, \
                         gk_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * gi0_161[k]
                   - f_2 * gi1_161[k]
                   + pb_y[k] * gk_208[k];

        t_262[k] = f_14 * fk_100[k]
                   + pb_z[k] * gk_208[k];

        t_263[k] = f_11 * gi0_163[k]
                   - f_12 * gi1_163[k]
                   + pb_y[k] * gk_210[k];

        t_264[k] = f_9 * gi0_164[k]
                   - f_10 * gi1_164[k]
                   + pb_y[k] * gk_211[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_y, gi0_165, gi0_166, gi0_167, gi1_165, \
                         gi1_166, gi1_167, gk_212, gk_213, gk_214, \
                         gk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_7 * gi0_165[k]
                   - f_8 * gi1_165[k]
                   + pb_y[k] * gk_212[k];

        t_266[k] = f_5 * gi0_166[k]
                   - f_6 * gi1_166[k]
                   + pb_y[k] * gk_213[k];

        t_267[k] = f_3 * gi0_167[k]
                   - f_4 * gi1_167[k]
                   + pb_y[k] * gk_214[k];

        t_268[k] = pb_y[k] * gk_215[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pb_y, pb_z, dl0_269, dl1_269, \
                         fk_108, fk_216, fl_269, fl_270, gk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_19 * dl0_269[k]
                   - f_20 * dl1_269[k]
                   + pa_x[k] * fl_269[k];

        t_270[k] = f_18 * fk_216[k]
                   + pa_x[k] * fl_270[k];

        t_271[k] = f_15 * fk_108[k]
                   + pb_y[k] * gk_216[k];

        t_272[k] = pb_z[k] * gk_216[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, pa_x, pb_z, fk_219, fk_221, \
                         fk_222, fl_273, fl_275, fl_276, gk_217, \
                         gk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_17 * fk_219[k]
                   + pa_x[k] * fl_273[k];

        t_274[k] = pb_z[k] * gk_217[k];

        t_275[k] = f_17 * fk_221[k]
                   + pa_x[k] * fl_275[k];

        t_276[k] = f_16 * fk_222[k]
                   + pa_x[k] * fl_276[k];

        t_277[k] = pb_z[k] * gk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_x, pb_y, pb_z, fk_113, fk_225, fk_226, \
                         fl_279, fl_280, gk_221, gk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_15 * fk_113[k]
                   + pb_y[k] * gk_221[k];

        t_279[k] = f_16 * fk_225[k]
                   + pa_x[k] * fl_279[k];

        t_280[k] = f_0 * fk_226[k]
                   + pa_x[k] * fl_280[k];

        t_281[k] = pb_z[k] * gk_222[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_x, pb_y, fk_117, fk_228, fk_230, \
                         fk_231, fl_282, fl_284, fl_285, gk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * fk_228[k]
                   + pa_x[k] * fl_282[k];

        t_283[k] = f_15 * fk_117[k]
                   + pb_y[k] * gk_225[k];

        t_284[k] = f_0 * fk_230[k]
                   + pa_x[k] * fl_284[k];

        t_285[k] = f_15 * fk_231[k]
                   + pa_x[k] * fl_285[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_x, pb_y, pb_z, fk_122, fk_233, fk_234, \
                         fl_287, fl_288, gk_226, gk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = pb_z[k] * gk_226[k];

        t_287[k] = f_15 * fk_233[k]
                   + pa_x[k] * fl_287[k];

        t_288[k] = f_15 * fk_234[k]
                   + pa_x[k] * fl_288[k];

        t_289[k] = f_15 * fk_122[k]
                   + pb_y[k] * gk_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pa_x, pb_z, fk_236, fk_237, \
                         fk_239, fk_240, fl_290, fl_291, fl_293, fl_294, \
                         gk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_15 * fk_236[k]
                   + pa_x[k] * fl_290[k];

        t_291[k] = f_14 * fk_237[k]
                   + pa_x[k] * fl_291[k];

        t_292[k] = pb_z[k] * gk_231[k];

        t_293[k] = f_14 * fk_239[k]
                   + pa_x[k] * fl_293[k];

        t_294[k] = f_14 * fk_240[k]
                   + pa_x[k] * fl_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pa_x, pb_x, pb_y, fk_128, fk_241, fk_243, \
                         fk_244, fl_295, fl_297, gk_236, gk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_14 * fk_241[k]
                   + pa_x[k] * fl_295[k];

        t_296[k] = f_15 * fk_128[k]
                   + pb_y[k] * gk_236[k];

        t_297[k] = f_14 * fk_243[k]
                   + pa_x[k] * fl_297[k];

        t_298[k] = f_13 * fk_244[k]
                   + pb_x[k] * gk_244[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pb_x, pb_z, fk_246, fk_247, \
                         fk_248, fk_249, gk_237, gk_246, gk_247, gk_248, \
                         gk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pb_z[k] * gk_237[k];

        t_300[k] = f_13 * fk_246[k]
                   + pb_x[k] * gk_246[k];

        t_301[k] = f_13 * fk_247[k]
                   + pb_x[k] * gk_247[k];

        t_302[k] = f_13 * fk_248[k]
                   + pb_x[k] * gk_248[k];

        t_303[k] = f_13 * fk_249[k]
                   + pb_x[k] * gk_249[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pa_x, pb_x, pb_z, fk_250, fk_251, \
                         fl_306, fl_308, gk_244, gk_250, gk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_13 * fk_250[k]
                   + pb_x[k] * gk_250[k];

        t_305[k] = f_13 * fk_251[k]
                   + pb_x[k] * gk_251[k];

        t_306[k] = pa_x[k] * fl_306[k];

        t_307[k] = pb_z[k] * gk_244[k];

        t_308[k] = pa_x[k] * fl_308[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, t_314, t_315, pa_x, pa_z, fl_135, \
                         fl_309, fl_310, fl_311, fl_312, fl_313, \
                         fl_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_x[k] * fl_309[k];

        t_310[k] = pa_x[k] * fl_310[k];

        t_311[k] = pa_x[k] * fl_311[k];

        t_312[k] = pa_x[k] * fl_312[k];

        t_313[k] = pa_x[k] * fl_313[k];

        t_314[k] = pa_x[k] * fl_314[k];

        t_315[k] = pa_z[k] * fl_135[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_z, pb_y, pb_z, fk_108, fk_146, fl_136, \
                         fl_138, gk_252, gk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_z[k] * fl_136[k];

        t_317[k] = f_13 * fk_108[k]
                   + pb_z[k] * gk_252[k];

        t_318[k] = pa_z[k] * fl_138[k];

        t_319[k] = f_14 * fk_146[k]
                   + pb_y[k] * gk_254[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_x, pa_z, pb_y, pb_z, fk_111, fk_149, \
                         fk_257, fl_141, fl_320, gk_255, gk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_17 * fk_257[k]
                   + pa_x[k] * fl_320[k];

        t_321[k] = pa_z[k] * fl_141[k];

        t_322[k] = f_13 * fk_111[k]
                   + pb_z[k] * gk_255[k];

        t_323[k] = f_14 * fk_149[k]
                   + pb_y[k] * gk_257[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_x, pa_z, pb_z, fk_114, fk_261, fk_264, \
                         fl_145, fl_324, fl_327, gk_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_16 * fk_261[k]
                   + pa_x[k] * fl_324[k];

        t_325[k] = pa_z[k] * fl_145[k];

        t_326[k] = f_13 * fk_114[k]
                   + pb_z[k] * gk_258[k];

        t_327[k] = f_0 * fk_264[k]
                   + pa_x[k] * fl_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_x, pa_z, pb_y, pb_z, fk_118, fk_153, \
                         fk_266, fl_150, fl_329, gk_261, gk_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_14 * fk_153[k]
                   + pb_y[k] * gk_261[k];

        t_329[k] = f_0 * fk_266[k]
                   + pa_x[k] * fl_329[k];

        t_330[k] = pa_z[k] * fl_150[k];

        t_331[k] = f_13 * fk_118[k]
                   + pb_z[k] * gk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pa_x, pb_y, fk_158, fk_269, fk_270, \
                         fk_272, fl_332, fl_333, fl_335, gk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_15 * fk_269[k]
                   + pa_x[k] * fl_332[k];

        t_333[k] = f_15 * fk_270[k]
                   + pa_x[k] * fl_333[k];

        t_334[k] = f_14 * fk_158[k]
                   + pb_y[k] * gk_266[k];

        t_335[k] = f_15 * fk_272[k]
                   + pa_x[k] * fl_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_x, pa_z, pb_z, fk_123, fk_275, fk_276, \
                         fl_156, fl_338, fl_339, gk_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_z[k] * fl_156[k];

        t_337[k] = f_13 * fk_123[k]
                   + pb_z[k] * gk_267[k];

        t_338[k] = f_14 * fk_275[k]
                   + pa_x[k] * fl_338[k];

        t_339[k] = f_14 * fk_276[k]
                   + pa_x[k] * fl_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_x, pa_z, pb_y, fk_164, fk_277, fk_279, \
                         fl_163, fl_340, fl_342, gk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_14 * fk_277[k]
                   + pa_x[k] * fl_340[k];

        t_341[k] = f_14 * fk_164[k]
                   + pb_y[k] * gk_272[k];

        t_342[k] = f_14 * fk_279[k]
                   + pa_x[k] * fl_342[k];

        t_343[k] = pa_z[k] * fl_163[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, pb_x, fk_281, fk_282, fk_283, \
                         fk_284, fk_285, gk_281, gk_282, gk_283, gk_284, \
                         gk_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_13 * fk_281[k]
                   + pb_x[k] * gk_281[k];

        t_345[k] = f_13 * fk_282[k]
                   + pb_x[k] * gk_282[k];

        t_346[k] = f_13 * fk_283[k]
                   + pb_x[k] * gk_283[k];

        t_347[k] = f_13 * fk_284[k]
                   + pb_x[k] * gk_284[k];

        t_348[k] = f_13 * fk_285[k]
                   + pb_x[k] * gk_285[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, t_354, pa_x, pb_x, fk_286, fk_287, \
                         fl_351, fl_352, fl_353, fl_354, gk_286, \
                         gk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_13 * fk_286[k]
                   + pb_x[k] * gk_286[k];

        t_350[k] = f_13 * fk_287[k]
                   + pb_x[k] * gk_287[k];

        t_351[k] = pa_x[k] * fl_351[k];

        t_352[k] = pa_x[k] * fl_352[k];

        t_353[k] = pa_x[k] * fl_353[k];

        t_354[k] = pa_x[k] * fl_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, t_360, pa_x, pa_y, fl_225, fl_355, \
                         fl_356, fl_357, fl_358, fl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = pa_x[k] * fl_355[k];

        t_356[k] = pa_x[k] * fl_356[k];

        t_357[k] = pa_x[k] * fl_357[k];

        t_358[k] = pa_x[k] * fl_358[k];

        t_359[k] = pa_x[k] * fl_359[k];

        t_360[k] = pa_y[k] * fl_225[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, t_365, pa_x, pa_y, pb_y, fk_180, fk_182, \
                         fk_291, fl_227, fl_230, fl_363, gk_288, \
                         gk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_13 * fk_180[k]
                   + pb_y[k] * gk_288[k];

        t_362[k] = pa_y[k] * fl_227[k];

        t_363[k] = f_17 * fk_291[k]
                   + pa_x[k] * fl_363[k];

        t_364[k] = f_13 * fk_182[k]
                   + pb_y[k] * gk_290[k];

        t_365[k] = pa_y[k] * fl_230[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_x, pa_y, pb_y, pb_z, fk_147, fk_185, \
                         fk_294, fl_234, fl_366, gk_291, gk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_16 * fk_294[k]
                   + pa_x[k] * fl_366[k];

        t_367[k] = f_14 * fk_147[k]
                   + pb_z[k] * gk_291[k];

        t_368[k] = f_13 * fk_185[k]
                   + pb_y[k] * gk_293[k];

        t_369[k] = pa_y[k] * fl_234[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pb_y, pb_z, fk_150, fk_189, fk_298, \
                         fk_300, fl_370, fl_372, gk_294, gk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_0 * fk_298[k]
                   + pa_x[k] * fl_370[k];

        t_371[k] = f_14 * fk_150[k]
                   + pb_z[k] * gk_294[k];

        t_372[k] = f_0 * fk_300[k]
                   + pa_x[k] * fl_372[k];

        t_373[k] = f_13 * fk_189[k]
                   + pb_y[k] * gk_297[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pa_x, pa_y, pb_z, fk_154, fk_303, fk_305, \
                         fl_239, fl_375, fl_377, gk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * fl_239[k];

        t_375[k] = f_15 * fk_303[k]
                   + pa_x[k] * fl_375[k];

        t_376[k] = f_14 * fk_154[k]
                   + pb_z[k] * gk_298[k];

        t_377[k] = f_15 * fk_305[k]
                   + pa_x[k] * fl_377[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_x, pa_y, pb_y, fk_194, fk_306, fk_309, \
                         fl_245, fl_378, fl_381, gk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_15 * fk_306[k]
                   + pa_x[k] * fl_378[k];

        t_379[k] = f_13 * fk_194[k]
                   + pb_y[k] * gk_302[k];

        t_380[k] = pa_y[k] * fl_245[k];

        t_381[k] = f_14 * fk_309[k]
                   + pa_x[k] * fl_381[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_x, pb_z, fk_159, fk_311, fk_312, \
                         fk_313, fl_383, fl_384, fl_385, gk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_14 * fk_159[k]
                   + pb_z[k] * gk_303[k];

        t_383[k] = f_14 * fk_311[k]
                   + pa_x[k] * fl_383[k];

        t_384[k] = f_14 * fk_312[k]
                   + pa_x[k] * fl_384[k];

        t_385[k] = f_14 * fk_313[k]
                   + pa_x[k] * fl_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_y, pb_x, pb_y, fk_200, fk_316, fk_317, \
                         fl_252, gk_308, gk_316, gk_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_13 * fk_200[k]
                   + pb_y[k] * gk_308[k];

        t_387[k] = pa_y[k] * fl_252[k];

        t_388[k] = f_13 * fk_316[k]
                   + pb_x[k] * gk_316[k];

        t_389[k] = f_13 * fk_317[k]
                   + pb_x[k] * gk_317[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, pb_x, fk_318, fk_319, fk_320, \
                         fk_321, fk_322, gk_318, gk_319, gk_320, gk_321, \
                         gk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_13 * fk_318[k]
                   + pb_x[k] * gk_318[k];

        t_391[k] = f_13 * fk_319[k]
                   + pb_x[k] * gk_319[k];

        t_392[k] = f_13 * fk_320[k]
                   + pb_x[k] * gk_320[k];

        t_393[k] = f_13 * fk_321[k]
                   + pb_x[k] * gk_321[k];

        t_394[k] = f_13 * fk_322[k]
                   + pb_x[k] * gk_322[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, t_400, t_401, pa_x, pa_y, fl_260, \
                         fl_396, fl_397, fl_398, fl_399, fl_400, \
                         fl_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * fl_260[k];

        t_396[k] = pa_x[k] * fl_396[k];

        t_397[k] = pa_x[k] * fl_397[k];

        t_398[k] = pa_x[k] * fl_398[k];

        t_399[k] = pa_x[k] * fl_399[k];

        t_400[k] = pa_x[k] * fl_400[k];

        t_401[k] = pa_x[k] * fl_401[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, t_407, pa_x, pb_y, pb_z, fk_180, \
                         fk_324, fl_402, fl_403, fl_404, fl_405, \
                         gk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_x[k] * fl_402[k];

        t_403[k] = pa_x[k] * fl_403[k];

        t_404[k] = pa_x[k] * fl_404[k];

        t_405[k] = f_18 * fk_324[k]
                   + pa_x[k] * fl_405[k];

        t_406[k] = pb_y[k] * gk_324[k];

        t_407[k] = f_15 * fk_180[k]
                   + pb_z[k] * gk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pa_x, pb_y, fk_327, fk_329, fk_330, \
                         fl_408, fl_410, fl_411, gk_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_17 * fk_327[k]
                   + pa_x[k] * fl_408[k];

        t_409[k] = pb_y[k] * gk_326[k];

        t_410[k] = f_17 * fk_329[k]
                   + pa_x[k] * fl_410[k];

        t_411[k] = f_16 * fk_330[k]
                   + pa_x[k] * fl_411[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_x, pb_y, pb_z, fk_183, fk_333, fk_334, \
                         fl_414, fl_415, gk_327, gk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_15 * fk_183[k]
                   + pb_z[k] * gk_327[k];

        t_413[k] = pb_y[k] * gk_329[k];

        t_414[k] = f_16 * fk_333[k]
                   + pa_x[k] * fl_414[k];

        t_415[k] = f_0 * fk_334[k]
                   + pa_x[k] * fl_415[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pa_x, pb_y, pb_z, fk_186, fk_336, fk_338, \
                         fl_417, fl_419, gk_330, gk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_15 * fk_186[k]
                   + pb_z[k] * gk_330[k];

        t_417[k] = f_0 * fk_336[k]
                   + pa_x[k] * fl_417[k];

        t_418[k] = pb_y[k] * gk_333[k];

        t_419[k] = f_0 * fk_338[k]
                   + pa_x[k] * fl_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pa_x, pb_z, fk_190, fk_339, fk_341, \
                         fk_342, fl_420, fl_422, fl_423, gk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_15 * fk_339[k]
                   + pa_x[k] * fl_420[k];

        t_421[k] = f_15 * fk_190[k]
                   + pb_z[k] * gk_334[k];

        t_422[k] = f_15 * fk_341[k]
                   + pa_x[k] * fl_422[k];

        t_423[k] = f_15 * fk_342[k]
                   + pa_x[k] * fl_423[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pa_x, pb_y, pb_z, fk_195, fk_344, fk_345, \
                         fl_425, fl_426, gk_338, gk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_y[k] * gk_338[k];

        t_425[k] = f_15 * fk_344[k]
                   + pa_x[k] * fl_425[k];

        t_426[k] = f_14 * fk_345[k]
                   + pa_x[k] * fl_426[k];

        t_427[k] = f_15 * fk_195[k]
                   + pb_z[k] * gk_339[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pa_x, pb_y, fk_347, fk_348, \
                         fk_349, fk_351, fl_428, fl_429, fl_430, fl_432, \
                         gk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_14 * fk_347[k]
                   + pa_x[k] * fl_428[k];

        t_429[k] = f_14 * fk_348[k]
                   + pa_x[k] * fl_429[k];

        t_430[k] = f_14 * fk_349[k]
                   + pa_x[k] * fl_430[k];

        t_431[k] = pb_y[k] * gk_344[k];

        t_432[k] = f_14 * fk_351[k]
                   + pa_x[k] * fl_432[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, pb_x, fk_352, fk_353, fk_354, \
                         fk_355, fk_356, gk_352, gk_353, gk_354, gk_355, \
                         gk_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_13 * fk_352[k]
                   + pb_x[k] * gk_352[k];

        t_434[k] = f_13 * fk_353[k]
                   + pb_x[k] * gk_353[k];

        t_435[k] = f_13 * fk_354[k]
                   + pb_x[k] * gk_354[k];

        t_436[k] = f_13 * fk_355[k]
                   + pb_x[k] * gk_355[k];

        t_437[k] = f_13 * fk_356[k]
                   + pb_x[k] * gk_356[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, pa_x, pb_x, pb_y, fk_357, fk_359, \
                         fl_441, fl_442, gk_351, gk_357, gk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_13 * fk_357[k]
                   + pb_x[k] * gk_357[k];

        t_439[k] = pb_y[k] * gk_351[k];

        t_440[k] = f_13 * fk_359[k]
                   + pb_x[k] * gk_359[k];

        t_441[k] = pa_x[k] * fl_441[k];

        t_442[k] = pa_x[k] * fl_442[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, t_449, pa_x, pb_y, fl_443, \
                         fl_444, fl_445, fl_446, fl_447, fl_449, \
                         gk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pa_x[k] * fl_443[k];

        t_444[k] = pa_x[k] * fl_444[k];

        t_445[k] = pa_x[k] * fl_445[k];

        t_446[k] = pa_x[k] * fl_446[k];

        t_447[k] = pa_x[k] * fl_447[k];

        t_448[k] = pb_y[k] * gk_359[k];

        t_449[k] = pa_x[k] * fl_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, pb_x, pb_y, pb_z, fk_216, gi0_280, \
                         gi0_283, gi1_280, gi1_283, gk_360, gk_361, \
                         gk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * gi0_280[k]
                   - f_2 * gi1_280[k]
                   + pb_x[k] * gk_360[k];

        t_451[k] = f_0 * fk_216[k]
                   + pb_y[k] * gk_360[k];

        t_452[k] = pb_z[k] * gk_360[k];

        t_453[k] = f_11 * gi0_283[k]
                   - f_12 * gi1_283[k]
                   + pb_x[k] * gk_363[k];

        t_454[k] = pb_z[k] * gk_361[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pb_x, pb_y, pb_z, fk_221, gi0_285, \
                         gi0_286, gi1_285, gi1_286, gk_363, gk_365, \
                         gk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * gi0_285[k]
                   - f_12 * gi1_285[k]
                   + pb_x[k] * gk_365[k];

        t_456[k] = f_9 * gi0_286[k]
                   - f_10 * gi1_286[k]
                   + pb_x[k] * gk_366[k];

        t_457[k] = pb_z[k] * gk_363[k];

        t_458[k] = f_0 * fk_221[k]
                   + pb_y[k] * gk_365[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, pb_x, pb_z, gi0_289, gi0_290, gi0_292, \
                         gi1_289, gi1_290, gi1_292, gk_366, gk_369, gk_370, \
                         gk_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_9 * gi0_289[k]
                   - f_10 * gi1_289[k]
                   + pb_x[k] * gk_369[k];

        t_460[k] = f_7 * gi0_290[k]
                   - f_8 * gi1_290[k]
                   + pb_x[k] * gk_370[k];

        t_461[k] = pb_z[k] * gk_366[k];

        t_462[k] = f_7 * gi0_292[k]
                   - f_8 * gi1_292[k]
                   + pb_x[k] * gk_372[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pb_x, pb_y, pb_z, fk_225, gi0_294, \
                         gi0_295, gi1_294, gi1_295, gk_369, gk_370, gk_374, \
                         gk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_0 * fk_225[k]
                   + pb_y[k] * gk_369[k];

        t_464[k] = f_7 * gi0_294[k]
                   - f_8 * gi1_294[k]
                   + pb_x[k] * gk_374[k];

        t_465[k] = f_5 * gi0_295[k]
                   - f_6 * gi1_295[k]
                   + pb_x[k] * gk_375[k];

        t_466[k] = pb_z[k] * gk_370[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_x, pb_y, fk_230, gi0_297, gi0_298, gi1_297, \
                         gi1_298, gk_374, gk_377, gk_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_5 * gi0_297[k]
                   - f_6 * gi1_297[k]
                   + pb_x[k] * gk_377[k];

        t_468[k] = f_5 * gi0_298[k]
                   - f_6 * gi1_298[k]
                   + pb_x[k] * gk_378[k];

        t_469[k] = f_0 * fk_230[k]
                   + pb_y[k] * gk_374[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_x, pb_z, gi0_300, gi0_301, gi0_303, \
                         gi1_300, gi1_301, gi1_303, gk_375, gk_380, gk_381, \
                         gk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_5 * gi0_300[k]
                   - f_6 * gi1_300[k]
                   + pb_x[k] * gk_380[k];

        t_471[k] = f_3 * gi0_301[k]
                   - f_4 * gi1_301[k]
                   + pb_x[k] * gk_381[k];

        t_472[k] = pb_z[k] * gk_375[k];

        t_473[k] = f_3 * gi0_303[k]
                   - f_4 * gi1_303[k]
                   + pb_x[k] * gk_383[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pb_x, pb_y, fk_236, gi0_304, gi0_305, gi1_304, \
                         gi1_305, gk_380, gk_384, gk_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_3 * gi0_304[k]
                   - f_4 * gi1_304[k]
                   + pb_x[k] * gk_384[k];

        t_475[k] = f_3 * gi0_305[k]
                   - f_4 * gi1_305[k]
                   + pb_x[k] * gk_385[k];

        t_476[k] = f_0 * fk_236[k]
                   + pb_y[k] * gk_380[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, t_482, pb_x, gi0_307, gi1_307, \
                         gk_387, gk_388, gk_389, gk_390, gk_391, \
                         gk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_3 * gi0_307[k]
                   - f_4 * gi1_307[k]
                   + pb_x[k] * gk_387[k];

        t_478[k] = pb_x[k] * gk_388[k];

        t_479[k] = pb_x[k] * gk_389[k];

        t_480[k] = pb_x[k] * gk_390[k];

        t_481[k] = pb_x[k] * gk_391[k];

        t_482[k] = pb_x[k] * gk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, pb_x, pb_y, pb_z, fk_244, gi0_301, \
                         gi1_301, gk_388, gk_393, gk_394, gk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pb_x[k] * gk_393[k];

        t_484[k] = pb_x[k] * gk_394[k];

        t_485[k] = pb_x[k] * gk_395[k];

        t_486[k] = f_0 * fk_244[k]
                   + f_1 * gi0_301[k]
                   - f_2 * gi1_301[k]
                   + pb_y[k] * gk_388[k];

        t_487[k] = pb_z[k] * gk_388[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pb_z, gi0_301, gi0_302, gi0_303, gi1_301, \
                         gi1_302, gi1_303, gk_389, gk_390, gk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_3 * gi0_301[k]
                   - f_4 * gi1_301[k]
                   + pb_z[k] * gk_389[k];

        t_489[k] = f_5 * gi0_302[k]
                   - f_6 * gi1_302[k]
                   + pb_z[k] * gk_390[k];

        t_490[k] = f_7 * gi0_303[k]
                   - f_8 * gi1_303[k]
                   + pb_z[k] * gk_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pb_y, pb_z, fk_251, gi0_304, gi0_305, \
                         gi0_307, gi1_304, gi1_305, gi1_307, gk_392, gk_393, \
                         gk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_9 * gi0_304[k]
                   - f_10 * gi1_304[k]
                   + pb_z[k] * gk_392[k];

        t_492[k] = f_11 * gi0_305[k]
                   - f_12 * gi1_305[k]
                   + pb_z[k] * gk_393[k];

        t_493[k] = f_0 * fk_251[k]
                   + pb_y[k] * gk_395[k];

        t_494[k] = f_1 * gi0_307[k]
                   - f_2 * gi1_307[k]
                   + pb_z[k] * gk_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, pa_z, pb_y, pb_z, fk_216, fk_254, \
                         fl_270, fl_271, fl_273, gk_396, gk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * fl_270[k];

        t_496[k] = pa_z[k] * fl_271[k];

        t_497[k] = f_13 * fk_216[k]
                   + pb_z[k] * gk_396[k];

        t_498[k] = pa_z[k] * fl_273[k];

        t_499[k] = f_15 * fk_254[k]
                   + pb_y[k] * gk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_z, pb_y, pb_z, fk_218, fk_219, fk_257, \
                         fl_275, fl_276, gk_399, gk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_14 * fk_218[k]
                   + pa_z[k] * fl_275[k];

        t_501[k] = pa_z[k] * fl_276[k];

        t_502[k] = f_13 * fk_219[k]
                   + pb_z[k] * gk_399[k];

        t_503[k] = f_15 * fk_257[k]
                   + pb_y[k] * gk_401[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_z, pb_z, fk_221, fk_222, fk_223, \
                         fl_279, fl_280, fl_282, gk_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_15 * fk_221[k]
                   + pa_z[k] * fl_279[k];

        t_505[k] = pa_z[k] * fl_280[k];

        t_506[k] = f_13 * fk_222[k]
                   + pb_z[k] * gk_402[k];

        t_507[k] = f_14 * fk_223[k]
                   + pa_z[k] * fl_282[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_z, pb_y, pb_z, fk_225, fk_226, fk_261, \
                         fl_284, fl_285, gk_405, gk_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_15 * fk_261[k]
                   + pb_y[k] * gk_405[k];

        t_509[k] = f_0 * fk_225[k]
                   + pa_z[k] * fl_284[k];

        t_510[k] = pa_z[k] * fl_285[k];

        t_511[k] = f_13 * fk_226[k]
                   + pb_z[k] * gk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pa_z, pb_y, fk_227, fk_228, \
                         fk_230, fk_266, fl_287, fl_288, fl_290, fl_291, \
                         gk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_14 * fk_227[k]
                   + pa_z[k] * fl_287[k];

        t_513[k] = f_15 * fk_228[k]
                   + pa_z[k] * fl_288[k];

        t_514[k] = f_15 * fk_266[k]
                   + pb_y[k] * gk_410[k];

        t_515[k] = f_16 * fk_230[k]
                   + pa_z[k] * fl_290[k];

        t_516[k] = pa_z[k] * fl_291[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_z, pb_z, fk_231, fk_232, fk_233, \
                         fk_234, fl_293, fl_294, fl_295, gk_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_13 * fk_231[k]
                   + pb_z[k] * gk_411[k];

        t_518[k] = f_14 * fk_232[k]
                   + pa_z[k] * fl_293[k];

        t_519[k] = f_15 * fk_233[k]
                   + pa_z[k] * fl_294[k];

        t_520[k] = f_0 * fk_234[k]
                   + pa_z[k] * fl_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, pa_z, pb_x, pb_y, fk_236, fk_272, \
                         fl_297, gk_416, gk_424, gk_425, gk_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_15 * fk_272[k]
                   + pb_y[k] * gk_416[k];

        t_522[k] = f_17 * fk_236[k]
                   + pa_z[k] * fl_297[k];

        t_523[k] = pb_x[k] * gk_424[k];

        t_524[k] = pb_x[k] * gk_425[k];

        t_525[k] = pb_x[k] * gk_426[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, t_531, pa_z, pb_x, fl_306, gk_427, \
                         gk_428, gk_429, gk_430, gk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = pb_x[k] * gk_427[k];

        t_527[k] = pb_x[k] * gk_428[k];

        t_528[k] = pb_x[k] * gk_429[k];

        t_529[k] = pb_x[k] * gk_430[k];

        t_530[k] = pb_x[k] * gk_431[k];

        t_531[k] = pa_z[k] * fl_306[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pa_z, pb_z, fk_244, fk_245, fk_246, \
                         fk_247, fl_308, fl_309, fl_310, gk_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_13 * fk_244[k]
                   + pb_z[k] * gk_424[k];

        t_533[k] = f_14 * fk_245[k]
                   + pa_z[k] * fl_308[k];

        t_534[k] = f_15 * fk_246[k]
                   + pa_z[k] * fl_309[k];

        t_535[k] = f_0 * fk_247[k]
                   + pa_z[k] * fl_310[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_z, pb_y, fk_248, fk_249, fk_251, \
                         fk_287, fl_311, fl_312, fl_314, gk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_16 * fk_248[k]
                   + pa_z[k] * fl_311[k];

        t_537[k] = f_17 * fk_249[k]
                   + pa_z[k] * fl_312[k];

        t_538[k] = f_15 * fk_287[k]
                   + pb_y[k] * gk_431[k];

        t_539[k] = f_18 * fk_251[k]
                   + pa_z[k] * fl_314[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pb_x, pb_y, pb_z, fk_252, fk_288, \
                         gi0_336, gi0_339, gi1_336, gi1_339, gk_432, \
                         gk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * gi0_336[k]
                   - f_2 * gi1_336[k]
                   + pb_x[k] * gk_432[k];

        t_541[k] = f_14 * fk_288[k]
                   + pb_y[k] * gk_432[k];

        t_542[k] = f_14 * fk_252[k]
                   + pb_z[k] * gk_432[k];

        t_543[k] = f_11 * gi0_339[k]
                   - f_12 * gi1_339[k]
                   + pb_x[k] * gk_435[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pb_x, pb_y, fk_290, gi0_341, gi0_342, gi1_341, \
                         gi1_342, gk_434, gk_437, gk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_14 * fk_290[k]
                   + pb_y[k] * gk_434[k];

        t_545[k] = f_11 * gi0_341[k]
                   - f_12 * gi1_341[k]
                   + pb_x[k] * gk_437[k];

        t_546[k] = f_9 * gi0_342[k]
                   - f_10 * gi1_342[k]
                   + pb_x[k] * gk_438[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pb_x, pb_y, pb_z, fk_255, fk_293, gi0_345, \
                         gi1_345, gk_435, gk_437, gk_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_14 * fk_255[k]
                   + pb_z[k] * gk_435[k];

        t_548[k] = f_14 * fk_293[k]
                   + pb_y[k] * gk_437[k];

        t_549[k] = f_9 * gi0_345[k]
                   - f_10 * gi1_345[k]
                   + pb_x[k] * gk_441[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pb_x, pb_z, fk_258, gi0_346, gi0_348, gi1_346, \
                         gi1_348, gk_438, gk_442, gk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_7 * gi0_346[k]
                   - f_8 * gi1_346[k]
                   + pb_x[k] * gk_442[k];

        t_551[k] = f_14 * fk_258[k]
                   + pb_z[k] * gk_438[k];

        t_552[k] = f_7 * gi0_348[k]
                   - f_8 * gi1_348[k]
                   + pb_x[k] * gk_444[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pb_x, pb_y, fk_297, gi0_350, gi0_351, gi1_350, \
                         gi1_351, gk_441, gk_446, gk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_14 * fk_297[k]
                   + pb_y[k] * gk_441[k];

        t_554[k] = f_7 * gi0_350[k]
                   - f_8 * gi1_350[k]
                   + pb_x[k] * gk_446[k];

        t_555[k] = f_5 * gi0_351[k]
                   - f_6 * gi1_351[k]
                   + pb_x[k] * gk_447[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pb_x, pb_z, fk_262, gi0_353, gi0_354, gi1_353, \
                         gi1_354, gk_442, gk_449, gk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_14 * fk_262[k]
                   + pb_z[k] * gk_442[k];

        t_557[k] = f_5 * gi0_353[k]
                   - f_6 * gi1_353[k]
                   + pb_x[k] * gk_449[k];

        t_558[k] = f_5 * gi0_354[k]
                   - f_6 * gi1_354[k]
                   + pb_x[k] * gk_450[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, pb_x, pb_y, fk_302, gi0_356, gi0_357, gi1_356, \
                         gi1_357, gk_446, gk_452, gk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_14 * fk_302[k]
                   + pb_y[k] * gk_446[k];

        t_560[k] = f_5 * gi0_356[k]
                   - f_6 * gi1_356[k]
                   + pb_x[k] * gk_452[k];

        t_561[k] = f_3 * gi0_357[k]
                   - f_4 * gi1_357[k]
                   + pb_x[k] * gk_453[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, pb_x, pb_z, fk_267, gi0_359, gi0_360, gi1_359, \
                         gi1_360, gk_447, gk_455, gk_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_14 * fk_267[k]
                   + pb_z[k] * gk_447[k];

        t_563[k] = f_3 * gi0_359[k]
                   - f_4 * gi1_359[k]
                   + pb_x[k] * gk_455[k];

        t_564[k] = f_3 * gi0_360[k]
                   - f_4 * gi1_360[k]
                   + pb_x[k] * gk_456[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pb_x, pb_y, fk_308, gi0_361, gi0_363, \
                         gi1_361, gi1_363, gk_452, gk_457, gk_459, \
                         gk_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_3 * gi0_361[k]
                   - f_4 * gi1_361[k]
                   + pb_x[k] * gk_457[k];

        t_566[k] = f_14 * fk_308[k]
                   + pb_y[k] * gk_452[k];

        t_567[k] = f_3 * gi0_363[k]
                   - f_4 * gi1_363[k]
                   + pb_x[k] * gk_459[k];

        t_568[k] = pb_x[k] * gk_460[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, t_574, t_575, pb_x, gk_461, \
                         gk_462, gk_463, gk_464, gk_465, gk_466, \
                         gk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = pb_x[k] * gk_461[k];

        t_570[k] = pb_x[k] * gk_462[k];

        t_571[k] = pb_x[k] * gk_463[k];

        t_572[k] = pb_x[k] * gk_464[k];

        t_573[k] = pb_x[k] * gk_465[k];

        t_574[k] = pb_x[k] * gk_466[k];

        t_575[k] = pb_x[k] * gk_467[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, pa_z, pb_y, pb_z, dl0_171, dl1_171, fk_280, \
                         fk_318, fl_351, gi0_359, gi1_359, gk_460, \
                         gk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_19 * dl0_171[k]
                   - f_20 * dl1_171[k]
                   + pa_z[k] * fl_351[k];

        t_577[k] = f_14 * fk_280[k]
                   + pb_z[k] * gk_460[k];

        t_578[k] = f_14 * fk_318[k]
                   + f_11 * gi0_359[k]
                   - f_12 * gi1_359[k]
                   + pb_y[k] * gk_462[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pb_y, fk_319, fk_320, fk_321, gi0_360, gi0_361, \
                         gi0_362, gi1_360, gi1_361, gi1_362, gk_463, gk_464, \
                         gk_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_14 * fk_319[k]
                   + f_9 * gi0_360[k]
                   - f_10 * gi1_360[k]
                   + pb_y[k] * gk_463[k];

        t_580[k] = f_14 * fk_320[k]
                   + f_7 * gi0_361[k]
                   - f_8 * gi1_361[k]
                   + pb_y[k] * gk_464[k];

        t_581[k] = f_14 * fk_321[k]
                   + f_5 * gi0_362[k]
                   - f_6 * gi1_362[k]
                   + pb_y[k] * gk_465[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pa_y, pb_y, dl0_269, dl1_269, fk_322, \
                         fk_323, fl_404, fl_405, gi0_363, gi1_363, gk_466, \
                         gk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_14 * fk_322[k]
                   + f_3 * gi0_363[k]
                   - f_4 * gi1_363[k]
                   + pb_y[k] * gk_466[k];

        t_583[k] = f_14 * fk_323[k]
                   + pb_y[k] * gk_467[k];

        t_584[k] = f_19 * dl0_269[k]
                   - f_20 * dl1_269[k]
                   + pa_y[k] * fl_404[k];

        t_585[k] = pa_y[k] * fl_405[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, pa_y, pb_y, fk_324, fk_325, \
                         fk_326, fl_407, fl_408, fl_410, gk_468, \
                         gk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_13 * fk_324[k]
                   + pb_y[k] * gk_468[k];

        t_587[k] = pa_y[k] * fl_407[k];

        t_588[k] = f_14 * fk_325[k]
                   + pa_y[k] * fl_408[k];

        t_589[k] = f_13 * fk_326[k]
                   + pb_y[k] * gk_470[k];

        t_590[k] = pa_y[k] * fl_410[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pa_y, pb_y, pb_z, fk_291, fk_327, fk_329, \
                         fl_411, fl_414, gk_471, gk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_15 * fk_327[k]
                   + pa_y[k] * fl_411[k];

        t_592[k] = f_15 * fk_291[k]
                   + pb_z[k] * gk_471[k];

        t_593[k] = f_13 * fk_329[k]
                   + pb_y[k] * gk_473[k];

        t_594[k] = pa_y[k] * fl_414[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pa_y, pb_y, pb_z, fk_294, fk_330, fk_332, \
                         fk_333, fl_415, fl_417, gk_474, gk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_0 * fk_330[k]
                   + pa_y[k] * fl_415[k];

        t_596[k] = f_15 * fk_294[k]
                   + pb_z[k] * gk_474[k];

        t_597[k] = f_14 * fk_332[k]
                   + pa_y[k] * fl_417[k];

        t_598[k] = f_13 * fk_333[k]
                   + pb_y[k] * gk_477[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pa_y, pb_z, fk_298, fk_334, \
                         fk_336, fk_337, fl_419, fl_420, fl_422, fl_423, \
                         gk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pa_y[k] * fl_419[k];

        t_600[k] = f_16 * fk_334[k]
                   + pa_y[k] * fl_420[k];

        t_601[k] = f_15 * fk_298[k]
                   + pb_z[k] * gk_478[k];

        t_602[k] = f_15 * fk_336[k]
                   + pa_y[k] * fl_422[k];

        t_603[k] = f_14 * fk_337[k]
                   + pa_y[k] * fl_423[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_y, pb_y, pb_z, fk_303, fk_338, fk_339, \
                         fl_425, fl_426, gk_482, gk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_13 * fk_338[k]
                   + pb_y[k] * gk_482[k];

        t_605[k] = pa_y[k] * fl_425[k];

        t_606[k] = f_17 * fk_339[k]
                   + pa_y[k] * fl_426[k];

        t_607[k] = f_15 * fk_303[k]
                   + pb_z[k] * gk_483[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, t_612, pa_y, pb_y, fk_341, fk_342, \
                         fk_343, fk_344, fl_428, fl_429, fl_430, fl_432, \
                         gk_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_0 * fk_341[k]
                   + pa_y[k] * fl_428[k];

        t_609[k] = f_15 * fk_342[k]
                   + pa_y[k] * fl_429[k];

        t_610[k] = f_14 * fk_343[k]
                   + pa_y[k] * fl_430[k];

        t_611[k] = f_13 * fk_344[k]
                   + pb_y[k] * gk_488[k];

        t_612[k] = pa_y[k] * fl_432[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, t_617, t_618, t_619, pb_x, gk_496, \
                         gk_497, gk_498, gk_499, gk_500, gk_501, \
                         gk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = pb_x[k] * gk_496[k];

        t_614[k] = pb_x[k] * gk_497[k];

        t_615[k] = pb_x[k] * gk_498[k];

        t_616[k] = pb_x[k] * gk_499[k];

        t_617[k] = pb_x[k] * gk_500[k];

        t_618[k] = pb_x[k] * gk_501[k];

        t_619[k] = pb_x[k] * gk_502[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pa_y, pb_x, pb_z, fk_316, fk_352, fk_354, \
                         fl_441, fl_443, gk_496, gk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = pb_x[k] * gk_503[k];

        t_621[k] = f_18 * fk_352[k]
                   + pa_y[k] * fl_441[k];

        t_622[k] = f_15 * fk_316[k]
                   + pb_z[k] * gk_496[k];

        t_623[k] = f_17 * fk_354[k]
                   + pa_y[k] * fl_443[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pa_y, fk_355, fk_356, fk_357, fk_358, \
                         fl_444, fl_445, fl_446, fl_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_16 * fk_355[k]
                   + pa_y[k] * fl_444[k];

        t_625[k] = f_0 * fk_356[k]
                   + pa_y[k] * fl_445[k];

        t_626[k] = f_15 * fk_357[k]
                   + pa_y[k] * fl_446[k];

        t_627[k] = f_14 * fk_358[k]
                   + pa_y[k] * fl_447[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, t_632, pa_y, pb_x, pb_y, pb_z, fk_324, \
                         fk_359, fl_449, gi0_392, gi1_392, gk_503, \
                         gk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_13 * fk_359[k]
                   + pb_y[k] * gk_503[k];

        t_629[k] = pa_y[k] * fl_449[k];

        t_630[k] = f_1 * gi0_392[k]
                   - f_2 * gi1_392[k]
                   + pb_x[k] * gk_504[k];

        t_631[k] = pb_y[k] * gk_504[k];

        t_632[k] = f_0 * fk_324[k]
                   + pb_z[k] * gk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pb_x, pb_y, gi0_395, gi0_397, gi0_398, \
                         gi1_395, gi1_397, gi1_398, gk_506, gk_507, gk_509, \
                         gk_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_11 * gi0_395[k]
                   - f_12 * gi1_395[k]
                   + pb_x[k] * gk_507[k];

        t_634[k] = pb_y[k] * gk_506[k];

        t_635[k] = f_11 * gi0_397[k]
                   - f_12 * gi1_397[k]
                   + pb_x[k] * gk_509[k];

        t_636[k] = f_9 * gi0_398[k]
                   - f_10 * gi1_398[k]
                   + pb_x[k] * gk_510[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pb_x, pb_y, pb_z, fk_327, gi0_401, \
                         gi0_402, gi1_401, gi1_402, gk_507, gk_509, gk_513, \
                         gk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = f_0 * fk_327[k]
                   + pb_z[k] * gk_507[k];

        t_638[k] = pb_y[k] * gk_509[k];

        t_639[k] = f_9 * gi0_401[k]
                   - f_10 * gi1_401[k]
                   + pb_x[k] * gk_513[k];

        t_640[k] = f_7 * gi0_402[k]
                   - f_8 * gi1_402[k]
                   + pb_x[k] * gk_514[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, pb_x, pb_y, pb_z, fk_330, gi0_404, \
                         gi0_406, gi1_404, gi1_406, gk_510, gk_513, gk_516, \
                         gk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_0 * fk_330[k]
                   + pb_z[k] * gk_510[k];

        t_642[k] = f_7 * gi0_404[k]
                   - f_8 * gi1_404[k]
                   + pb_x[k] * gk_516[k];

        t_643[k] = pb_y[k] * gk_513[k];

        t_644[k] = f_7 * gi0_406[k]
                   - f_8 * gi1_406[k]
                   + pb_x[k] * gk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pb_x, pb_z, fk_334, gi0_407, gi0_409, gi1_407, \
                         gi1_409, gk_514, gk_519, gk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_5 * gi0_407[k]
                   - f_6 * gi1_407[k]
                   + pb_x[k] * gk_519[k];

        t_646[k] = f_0 * fk_334[k]
                   + pb_z[k] * gk_514[k];

        t_647[k] = f_5 * gi0_409[k]
                   - f_6 * gi1_409[k]
                   + pb_x[k] * gk_521[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, t_651, pb_x, pb_y, gi0_410, gi0_412, gi0_413, \
                         gi1_410, gi1_412, gi1_413, gk_518, gk_522, gk_524, \
                         gk_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_5 * gi0_410[k]
                   - f_6 * gi1_410[k]
                   + pb_x[k] * gk_522[k];

        t_649[k] = pb_y[k] * gk_518[k];

        t_650[k] = f_5 * gi0_412[k]
                   - f_6 * gi1_412[k]
                   + pb_x[k] * gk_524[k];

        t_651[k] = f_3 * gi0_413[k]
                   - f_4 * gi1_413[k]
                   + pb_x[k] * gk_525[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pb_x, pb_z, fk_339, gi0_415, gi0_416, gi1_415, \
                         gi1_416, gk_519, gk_527, gk_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_0 * fk_339[k]
                   + pb_z[k] * gk_519[k];

        t_653[k] = f_3 * gi0_415[k]
                   - f_4 * gi1_415[k]
                   + pb_x[k] * gk_527[k];

        t_654[k] = f_3 * gi0_416[k]
                   - f_4 * gi1_416[k]
                   + pb_x[k] * gk_528[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, pb_x, pb_y, gi0_417, gi0_419, \
                         gi1_417, gi1_419, gk_524, gk_529, gk_531, gk_532, \
                         gk_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_3 * gi0_417[k]
                   - f_4 * gi1_417[k]
                   + pb_x[k] * gk_529[k];

        t_656[k] = pb_y[k] * gk_524[k];

        t_657[k] = f_3 * gi0_419[k]
                   - f_4 * gi1_419[k]
                   + pb_x[k] * gk_531[k];

        t_658[k] = pb_x[k] * gk_532[k];

        t_659[k] = pb_x[k] * gk_533[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, pb_x, gk_534, gk_535, \
                         gk_536, gk_537, gk_538, gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = pb_x[k] * gk_534[k];

        t_661[k] = pb_x[k] * gk_535[k];

        t_662[k] = pb_x[k] * gk_536[k];

        t_663[k] = pb_x[k] * gk_537[k];

        t_664[k] = pb_x[k] * gk_538[k];

        t_665[k] = pb_x[k] * gk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_y, pb_z, fk_352, gi0_413, gi0_415, \
                         gi0_416, gi1_413, gi1_415, gi1_416, gk_532, gk_534, \
                         gk_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * gi0_413[k]
                   - f_2 * gi1_413[k]
                   + pb_y[k] * gk_532[k];

        t_667[k] = f_0 * fk_352[k]
                   + pb_z[k] * gk_532[k];

        t_668[k] = f_11 * gi0_415[k]
                   - f_12 * gi1_415[k]
                   + pb_y[k] * gk_534[k];

        t_669[k] = f_9 * gi0_416[k]
                   - f_10 * gi1_416[k]
                   + pb_y[k] * gk_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_y, gi0_417, gi0_418, gi0_419, gi1_417, \
                         gi1_418, gi1_419, gk_536, gk_537, gk_538, \
                         gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_7 * gi0_417[k]
                   - f_8 * gi1_417[k]
                   + pb_y[k] * gk_536[k];

        t_671[k] = f_5 * gi0_418[k]
                   - f_6 * gi1_418[k]
                   + pb_y[k] * gk_537[k];

        t_672[k] = f_3 * gi0_419[k]
                   - f_4 * gi1_419[k]
                   + pb_y[k] * gk_538[k];

        t_673[k] = pb_y[k] * gk_539[k];
    }

#pragma omp simd aligned(t_674, pb_z, fk_359, gi0_419, gi1_419, \
                         gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_0 * fk_359[k]
                   + f_1 * gi0_419[k]
                   - f_2 * gi1_419[k]
                   + pb_z[k] * gk_539[k];
    }
}

}  // namespace simdt2ceri
