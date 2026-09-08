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


#include "SimdElectronRepulsionVrrRecHL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_hl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
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
    const auto f_16 = 2.0 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);

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
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_315 = buffer.data(gk + 315);
    const auto *gk_316 = buffer.data(gk + 316);
    const auto *gk_317 = buffer.data(gk + 317);
    const auto *gk_318 = buffer.data(gk + 318);
    const auto *gk_319 = buffer.data(gk + 319);
    const auto *gk_320 = buffer.data(gk + 320);
    const auto *gk_321 = buffer.data(gk + 321);
    const auto *gk_322 = buffer.data(gk + 322);
    const auto *gk_323 = buffer.data(gk + 323);
    const auto *gk_324 = buffer.data(gk + 324);
    const auto *gk_325 = buffer.data(gk + 325);
    const auto *gk_326 = buffer.data(gk + 326);
    const auto *gk_327 = buffer.data(gk + 327);
    const auto *gk_328 = buffer.data(gk + 328);
    const auto *gk_329 = buffer.data(gk + 329);
    const auto *gk_330 = buffer.data(gk + 330);
    const auto *gk_331 = buffer.data(gk + 331);
    const auto *gk_332 = buffer.data(gk + 332);
    const auto *gk_333 = buffer.data(gk + 333);
    const auto *gk_334 = buffer.data(gk + 334);
    const auto *gk_335 = buffer.data(gk + 335);
    const auto *gk_336 = buffer.data(gk + 336);
    const auto *gk_337 = buffer.data(gk + 337);
    const auto *gk_338 = buffer.data(gk + 338);
    const auto *gk_339 = buffer.data(gk + 339);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);
    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);
    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_1 = buffer.data(hi0 + 1);
    const auto *hi0_2 = buffer.data(hi0 + 2);
    const auto *hi0_3 = buffer.data(hi0 + 3);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_6 = buffer.data(hi0 + 6);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_9 = buffer.data(hi0 + 9);
    const auto *hi0_10 = buffer.data(hi0 + 10);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_15 = buffer.data(hi0 + 15);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_18 = buffer.data(hi0 + 18);
    const auto *hi0_19 = buffer.data(hi0 + 19);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_21 = buffer.data(hi0 + 21);
    const auto *hi0_22 = buffer.data(hi0 + 22);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_25 = buffer.data(hi0 + 25);
    const auto *hi0_26 = buffer.data(hi0 + 26);
    const auto *hi0_27 = buffer.data(hi0 + 27);
    const auto *hi0_28 = buffer.data(hi0 + 28);
    const auto *hi0_29 = buffer.data(hi0 + 29);
    const auto *hi0_30 = buffer.data(hi0 + 30);
    const auto *hi0_31 = buffer.data(hi0 + 31);
    const auto *hi0_32 = buffer.data(hi0 + 32);
    const auto *hi0_33 = buffer.data(hi0 + 33);
    const auto *hi0_34 = buffer.data(hi0 + 34);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_36 = buffer.data(hi0 + 36);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_38 = buffer.data(hi0 + 38);
    const auto *hi0_39 = buffer.data(hi0 + 39);
    const auto *hi0_40 = buffer.data(hi0 + 40);
    const auto *hi0_41 = buffer.data(hi0 + 41);
    const auto *hi0_42 = buffer.data(hi0 + 42);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_44 = buffer.data(hi0 + 44);
    const auto *hi0_45 = buffer.data(hi0 + 45);
    const auto *hi0_46 = buffer.data(hi0 + 46);
    const auto *hi0_47 = buffer.data(hi0 + 47);
    const auto *hi0_48 = buffer.data(hi0 + 48);
    const auto *hi0_49 = buffer.data(hi0 + 49);
    const auto *hi0_50 = buffer.data(hi0 + 50);
    const auto *hi0_51 = buffer.data(hi0 + 51);
    const auto *hi0_52 = buffer.data(hi0 + 52);
    const auto *hi0_53 = buffer.data(hi0 + 53);
    const auto *hi0_54 = buffer.data(hi0 + 54);
    const auto *hi0_55 = buffer.data(hi0 + 55);
    const auto *hi0_56 = buffer.data(hi0 + 56);
    const auto *hi0_57 = buffer.data(hi0 + 57);
    const auto *hi0_58 = buffer.data(hi0 + 58);
    const auto *hi0_59 = buffer.data(hi0 + 59);
    const auto *hi0_60 = buffer.data(hi0 + 60);
    const auto *hi0_61 = buffer.data(hi0 + 61);
    const auto *hi0_62 = buffer.data(hi0 + 62);
    const auto *hi0_63 = buffer.data(hi0 + 63);
    const auto *hi0_64 = buffer.data(hi0 + 64);
    const auto *hi0_65 = buffer.data(hi0 + 65);
    const auto *hi0_66 = buffer.data(hi0 + 66);
    const auto *hi0_67 = buffer.data(hi0 + 67);
    const auto *hi0_68 = buffer.data(hi0 + 68);
    const auto *hi0_69 = buffer.data(hi0 + 69);
    const auto *hi0_70 = buffer.data(hi0 + 70);
    const auto *hi0_71 = buffer.data(hi0 + 71);
    const auto *hi0_72 = buffer.data(hi0 + 72);
    const auto *hi0_73 = buffer.data(hi0 + 73);
    const auto *hi0_74 = buffer.data(hi0 + 74);
    const auto *hi0_75 = buffer.data(hi0 + 75);
    const auto *hi0_76 = buffer.data(hi0 + 76);
    const auto *hi0_77 = buffer.data(hi0 + 77);
    const auto *hi0_78 = buffer.data(hi0 + 78);
    const auto *hi0_79 = buffer.data(hi0 + 79);
    const auto *hi0_80 = buffer.data(hi0 + 80);
    const auto *hi0_81 = buffer.data(hi0 + 81);
    const auto *hi0_82 = buffer.data(hi0 + 82);
    const auto *hi0_83 = buffer.data(hi0 + 83);
    const auto *hi0_84 = buffer.data(hi0 + 84);
    const auto *hi0_85 = buffer.data(hi0 + 85);
    const auto *hi0_86 = buffer.data(hi0 + 86);
    const auto *hi0_87 = buffer.data(hi0 + 87);
    const auto *hi0_88 = buffer.data(hi0 + 88);
    const auto *hi0_89 = buffer.data(hi0 + 89);
    const auto *hi0_90 = buffer.data(hi0 + 90);
    const auto *hi0_91 = buffer.data(hi0 + 91);
    const auto *hi0_92 = buffer.data(hi0 + 92);
    const auto *hi0_93 = buffer.data(hi0 + 93);
    const auto *hi0_94 = buffer.data(hi0 + 94);
    const auto *hi0_95 = buffer.data(hi0 + 95);
    const auto *hi0_96 = buffer.data(hi0 + 96);
    const auto *hi0_97 = buffer.data(hi0 + 97);
    const auto *hi0_98 = buffer.data(hi0 + 98);
    const auto *hi0_99 = buffer.data(hi0 + 99);
    const auto *hi0_100 = buffer.data(hi0 + 100);
    const auto *hi0_101 = buffer.data(hi0 + 101);
    const auto *hi0_102 = buffer.data(hi0 + 102);
    const auto *hi0_103 = buffer.data(hi0 + 103);
    const auto *hi0_104 = buffer.data(hi0 + 104);
    const auto *hi0_105 = buffer.data(hi0 + 105);
    const auto *hi0_106 = buffer.data(hi0 + 106);
    const auto *hi0_107 = buffer.data(hi0 + 107);
    const auto *hi0_108 = buffer.data(hi0 + 108);
    const auto *hi0_109 = buffer.data(hi0 + 109);
    const auto *hi0_110 = buffer.data(hi0 + 110);
    const auto *hi0_111 = buffer.data(hi0 + 111);
    const auto *hi0_112 = buffer.data(hi0 + 112);
    const auto *hi0_113 = buffer.data(hi0 + 113);
    const auto *hi0_114 = buffer.data(hi0 + 114);
    const auto *hi0_115 = buffer.data(hi0 + 115);
    const auto *hi0_116 = buffer.data(hi0 + 116);
    const auto *hi0_117 = buffer.data(hi0 + 117);
    const auto *hi0_118 = buffer.data(hi0 + 118);
    const auto *hi0_119 = buffer.data(hi0 + 119);
    const auto *hi0_120 = buffer.data(hi0 + 120);
    const auto *hi0_121 = buffer.data(hi0 + 121);
    const auto *hi0_122 = buffer.data(hi0 + 122);
    const auto *hi0_123 = buffer.data(hi0 + 123);
    const auto *hi0_124 = buffer.data(hi0 + 124);
    const auto *hi0_125 = buffer.data(hi0 + 125);
    const auto *hi0_126 = buffer.data(hi0 + 126);
    const auto *hi0_127 = buffer.data(hi0 + 127);
    const auto *hi0_128 = buffer.data(hi0 + 128);
    const auto *hi0_129 = buffer.data(hi0 + 129);
    const auto *hi0_130 = buffer.data(hi0 + 130);
    const auto *hi0_131 = buffer.data(hi0 + 131);
    const auto *hi0_132 = buffer.data(hi0 + 132);
    const auto *hi0_133 = buffer.data(hi0 + 133);
    const auto *hi0_134 = buffer.data(hi0 + 134);
    const auto *hi0_135 = buffer.data(hi0 + 135);
    const auto *hi0_136 = buffer.data(hi0 + 136);
    const auto *hi0_137 = buffer.data(hi0 + 137);
    const auto *hi0_138 = buffer.data(hi0 + 138);
    const auto *hi0_139 = buffer.data(hi0 + 139);
    const auto *hi0_140 = buffer.data(hi0 + 140);
    const auto *hi0_141 = buffer.data(hi0 + 141);
    const auto *hi0_142 = buffer.data(hi0 + 142);
    const auto *hi0_143 = buffer.data(hi0 + 143);
    const auto *hi0_144 = buffer.data(hi0 + 144);
    const auto *hi0_145 = buffer.data(hi0 + 145);
    const auto *hi0_146 = buffer.data(hi0 + 146);
    const auto *hi0_147 = buffer.data(hi0 + 147);
    const auto *hi0_148 = buffer.data(hi0 + 148);
    const auto *hi0_149 = buffer.data(hi0 + 149);
    const auto *hi0_150 = buffer.data(hi0 + 150);
    const auto *hi0_151 = buffer.data(hi0 + 151);
    const auto *hi0_152 = buffer.data(hi0 + 152);
    const auto *hi0_153 = buffer.data(hi0 + 153);
    const auto *hi0_154 = buffer.data(hi0 + 154);
    const auto *hi0_155 = buffer.data(hi0 + 155);
    const auto *hi0_156 = buffer.data(hi0 + 156);
    const auto *hi0_157 = buffer.data(hi0 + 157);
    const auto *hi0_158 = buffer.data(hi0 + 158);
    const auto *hi0_159 = buffer.data(hi0 + 159);
    const auto *hi0_160 = buffer.data(hi0 + 160);
    const auto *hi0_161 = buffer.data(hi0 + 161);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_1 = buffer.data(hi1 + 1);
    const auto *hi1_2 = buffer.data(hi1 + 2);
    const auto *hi1_3 = buffer.data(hi1 + 3);
    const auto *hi1_4 = buffer.data(hi1 + 4);
    const auto *hi1_5 = buffer.data(hi1 + 5);
    const auto *hi1_6 = buffer.data(hi1 + 6);
    const auto *hi1_7 = buffer.data(hi1 + 7);
    const auto *hi1_8 = buffer.data(hi1 + 8);
    const auto *hi1_9 = buffer.data(hi1 + 9);
    const auto *hi1_10 = buffer.data(hi1 + 10);
    const auto *hi1_11 = buffer.data(hi1 + 11);
    const auto *hi1_12 = buffer.data(hi1 + 12);
    const auto *hi1_13 = buffer.data(hi1 + 13);
    const auto *hi1_14 = buffer.data(hi1 + 14);
    const auto *hi1_15 = buffer.data(hi1 + 15);
    const auto *hi1_16 = buffer.data(hi1 + 16);
    const auto *hi1_17 = buffer.data(hi1 + 17);
    const auto *hi1_18 = buffer.data(hi1 + 18);
    const auto *hi1_19 = buffer.data(hi1 + 19);
    const auto *hi1_20 = buffer.data(hi1 + 20);
    const auto *hi1_21 = buffer.data(hi1 + 21);
    const auto *hi1_22 = buffer.data(hi1 + 22);
    const auto *hi1_23 = buffer.data(hi1 + 23);
    const auto *hi1_24 = buffer.data(hi1 + 24);
    const auto *hi1_25 = buffer.data(hi1 + 25);
    const auto *hi1_26 = buffer.data(hi1 + 26);
    const auto *hi1_27 = buffer.data(hi1 + 27);
    const auto *hi1_28 = buffer.data(hi1 + 28);
    const auto *hi1_29 = buffer.data(hi1 + 29);
    const auto *hi1_30 = buffer.data(hi1 + 30);
    const auto *hi1_31 = buffer.data(hi1 + 31);
    const auto *hi1_32 = buffer.data(hi1 + 32);
    const auto *hi1_33 = buffer.data(hi1 + 33);
    const auto *hi1_34 = buffer.data(hi1 + 34);
    const auto *hi1_35 = buffer.data(hi1 + 35);
    const auto *hi1_36 = buffer.data(hi1 + 36);
    const auto *hi1_37 = buffer.data(hi1 + 37);
    const auto *hi1_38 = buffer.data(hi1 + 38);
    const auto *hi1_39 = buffer.data(hi1 + 39);
    const auto *hi1_40 = buffer.data(hi1 + 40);
    const auto *hi1_41 = buffer.data(hi1 + 41);
    const auto *hi1_42 = buffer.data(hi1 + 42);
    const auto *hi1_43 = buffer.data(hi1 + 43);
    const auto *hi1_44 = buffer.data(hi1 + 44);
    const auto *hi1_45 = buffer.data(hi1 + 45);
    const auto *hi1_46 = buffer.data(hi1 + 46);
    const auto *hi1_47 = buffer.data(hi1 + 47);
    const auto *hi1_48 = buffer.data(hi1 + 48);
    const auto *hi1_49 = buffer.data(hi1 + 49);
    const auto *hi1_50 = buffer.data(hi1 + 50);
    const auto *hi1_51 = buffer.data(hi1 + 51);
    const auto *hi1_52 = buffer.data(hi1 + 52);
    const auto *hi1_53 = buffer.data(hi1 + 53);
    const auto *hi1_54 = buffer.data(hi1 + 54);
    const auto *hi1_55 = buffer.data(hi1 + 55);
    const auto *hi1_56 = buffer.data(hi1 + 56);
    const auto *hi1_57 = buffer.data(hi1 + 57);
    const auto *hi1_58 = buffer.data(hi1 + 58);
    const auto *hi1_59 = buffer.data(hi1 + 59);
    const auto *hi1_60 = buffer.data(hi1 + 60);
    const auto *hi1_61 = buffer.data(hi1 + 61);
    const auto *hi1_62 = buffer.data(hi1 + 62);
    const auto *hi1_63 = buffer.data(hi1 + 63);
    const auto *hi1_64 = buffer.data(hi1 + 64);
    const auto *hi1_65 = buffer.data(hi1 + 65);
    const auto *hi1_66 = buffer.data(hi1 + 66);
    const auto *hi1_67 = buffer.data(hi1 + 67);
    const auto *hi1_68 = buffer.data(hi1 + 68);
    const auto *hi1_69 = buffer.data(hi1 + 69);
    const auto *hi1_70 = buffer.data(hi1 + 70);
    const auto *hi1_71 = buffer.data(hi1 + 71);
    const auto *hi1_72 = buffer.data(hi1 + 72);
    const auto *hi1_73 = buffer.data(hi1 + 73);
    const auto *hi1_74 = buffer.data(hi1 + 74);
    const auto *hi1_75 = buffer.data(hi1 + 75);
    const auto *hi1_76 = buffer.data(hi1 + 76);
    const auto *hi1_77 = buffer.data(hi1 + 77);
    const auto *hi1_78 = buffer.data(hi1 + 78);
    const auto *hi1_79 = buffer.data(hi1 + 79);
    const auto *hi1_80 = buffer.data(hi1 + 80);
    const auto *hi1_81 = buffer.data(hi1 + 81);
    const auto *hi1_82 = buffer.data(hi1 + 82);
    const auto *hi1_83 = buffer.data(hi1 + 83);
    const auto *hi1_84 = buffer.data(hi1 + 84);
    const auto *hi1_85 = buffer.data(hi1 + 85);
    const auto *hi1_86 = buffer.data(hi1 + 86);
    const auto *hi1_87 = buffer.data(hi1 + 87);
    const auto *hi1_88 = buffer.data(hi1 + 88);
    const auto *hi1_89 = buffer.data(hi1 + 89);
    const auto *hi1_90 = buffer.data(hi1 + 90);
    const auto *hi1_91 = buffer.data(hi1 + 91);
    const auto *hi1_92 = buffer.data(hi1 + 92);
    const auto *hi1_93 = buffer.data(hi1 + 93);
    const auto *hi1_94 = buffer.data(hi1 + 94);
    const auto *hi1_95 = buffer.data(hi1 + 95);
    const auto *hi1_96 = buffer.data(hi1 + 96);
    const auto *hi1_97 = buffer.data(hi1 + 97);
    const auto *hi1_98 = buffer.data(hi1 + 98);
    const auto *hi1_99 = buffer.data(hi1 + 99);
    const auto *hi1_100 = buffer.data(hi1 + 100);
    const auto *hi1_101 = buffer.data(hi1 + 101);
    const auto *hi1_102 = buffer.data(hi1 + 102);
    const auto *hi1_103 = buffer.data(hi1 + 103);
    const auto *hi1_104 = buffer.data(hi1 + 104);
    const auto *hi1_105 = buffer.data(hi1 + 105);
    const auto *hi1_106 = buffer.data(hi1 + 106);
    const auto *hi1_107 = buffer.data(hi1 + 107);
    const auto *hi1_108 = buffer.data(hi1 + 108);
    const auto *hi1_109 = buffer.data(hi1 + 109);
    const auto *hi1_110 = buffer.data(hi1 + 110);
    const auto *hi1_111 = buffer.data(hi1 + 111);
    const auto *hi1_112 = buffer.data(hi1 + 112);
    const auto *hi1_113 = buffer.data(hi1 + 113);
    const auto *hi1_114 = buffer.data(hi1 + 114);
    const auto *hi1_115 = buffer.data(hi1 + 115);
    const auto *hi1_116 = buffer.data(hi1 + 116);
    const auto *hi1_117 = buffer.data(hi1 + 117);
    const auto *hi1_118 = buffer.data(hi1 + 118);
    const auto *hi1_119 = buffer.data(hi1 + 119);
    const auto *hi1_120 = buffer.data(hi1 + 120);
    const auto *hi1_121 = buffer.data(hi1 + 121);
    const auto *hi1_122 = buffer.data(hi1 + 122);
    const auto *hi1_123 = buffer.data(hi1 + 123);
    const auto *hi1_124 = buffer.data(hi1 + 124);
    const auto *hi1_125 = buffer.data(hi1 + 125);
    const auto *hi1_126 = buffer.data(hi1 + 126);
    const auto *hi1_127 = buffer.data(hi1 + 127);
    const auto *hi1_128 = buffer.data(hi1 + 128);
    const auto *hi1_129 = buffer.data(hi1 + 129);
    const auto *hi1_130 = buffer.data(hi1 + 130);
    const auto *hi1_131 = buffer.data(hi1 + 131);
    const auto *hi1_132 = buffer.data(hi1 + 132);
    const auto *hi1_133 = buffer.data(hi1 + 133);
    const auto *hi1_134 = buffer.data(hi1 + 134);
    const auto *hi1_135 = buffer.data(hi1 + 135);
    const auto *hi1_136 = buffer.data(hi1 + 136);
    const auto *hi1_137 = buffer.data(hi1 + 137);
    const auto *hi1_138 = buffer.data(hi1 + 138);
    const auto *hi1_139 = buffer.data(hi1 + 139);
    const auto *hi1_140 = buffer.data(hi1 + 140);
    const auto *hi1_141 = buffer.data(hi1 + 141);
    const auto *hi1_142 = buffer.data(hi1 + 142);
    const auto *hi1_143 = buffer.data(hi1 + 143);
    const auto *hi1_144 = buffer.data(hi1 + 144);
    const auto *hi1_145 = buffer.data(hi1 + 145);
    const auto *hi1_146 = buffer.data(hi1 + 146);
    const auto *hi1_147 = buffer.data(hi1 + 147);
    const auto *hi1_148 = buffer.data(hi1 + 148);
    const auto *hi1_149 = buffer.data(hi1 + 149);
    const auto *hi1_150 = buffer.data(hi1 + 150);
    const auto *hi1_151 = buffer.data(hi1 + 151);
    const auto *hi1_152 = buffer.data(hi1 + 152);
    const auto *hi1_153 = buffer.data(hi1 + 153);
    const auto *hi1_154 = buffer.data(hi1 + 154);
    const auto *hi1_155 = buffer.data(hi1 + 155);
    const auto *hi1_156 = buffer.data(hi1 + 156);
    const auto *hi1_157 = buffer.data(hi1 + 157);
    const auto *hi1_158 = buffer.data(hi1 + 158);
    const auto *hi1_159 = buffer.data(hi1 + 159);
    const auto *hi1_160 = buffer.data(hi1 + 160);
    const auto *hi1_161 = buffer.data(hi1 + 161);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_24 = buffer.data(hk + 24);
    const auto *hk_25 = buffer.data(hk + 25);
    const auto *hk_26 = buffer.data(hk + 26);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_39 = buffer.data(hk + 39);
    const auto *hk_40 = buffer.data(hk + 40);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_47 = buffer.data(hk + 47);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_95 = buffer.data(hk + 95);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_116 = buffer.data(hk + 116);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_121 = buffer.data(hk + 121);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);
    const auto *hk_135 = buffer.data(hk + 135);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_144 = buffer.data(hk + 144);
    const auto *hk_145 = buffer.data(hk + 145);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_148 = buffer.data(hk + 148);
    const auto *hk_149 = buffer.data(hk + 149);
    const auto *hk_150 = buffer.data(hk + 150);
    const auto *hk_151 = buffer.data(hk + 151);
    const auto *hk_152 = buffer.data(hk + 152);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_155 = buffer.data(hk + 155);
    const auto *hk_156 = buffer.data(hk + 156);
    const auto *hk_157 = buffer.data(hk + 157);
    const auto *hk_158 = buffer.data(hk + 158);
    const auto *hk_159 = buffer.data(hk + 159);
    const auto *hk_160 = buffer.data(hk + 160);
    const auto *hk_161 = buffer.data(hk + 161);
    const auto *hk_162 = buffer.data(hk + 162);
    const auto *hk_163 = buffer.data(hk + 163);
    const auto *hk_164 = buffer.data(hk + 164);
    const auto *hk_165 = buffer.data(hk + 165);
    const auto *hk_166 = buffer.data(hk + 166);
    const auto *hk_167 = buffer.data(hk + 167);
    const auto *hk_168 = buffer.data(hk + 168);
    const auto *hk_169 = buffer.data(hk + 169);
    const auto *hk_170 = buffer.data(hk + 170);
    const auto *hk_171 = buffer.data(hk + 171);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_184 = buffer.data(hk + 184);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_187 = buffer.data(hk + 187);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_193 = buffer.data(hk + 193);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_196 = buffer.data(hk + 196);
    const auto *hk_197 = buffer.data(hk + 197);
    const auto *hk_198 = buffer.data(hk + 198);
    const auto *hk_199 = buffer.data(hk + 199);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_201 = buffer.data(hk + 201);
    const auto *hk_202 = buffer.data(hk + 202);
    const auto *hk_203 = buffer.data(hk + 203);
    const auto *hk_204 = buffer.data(hk + 204);
    const auto *hk_205 = buffer.data(hk + 205);
    const auto *hk_206 = buffer.data(hk + 206);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_217 = buffer.data(hk + 217);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_220 = buffer.data(hk + 220);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_224 = buffer.data(hk + 224);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_229 = buffer.data(hk + 229);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_235 = buffer.data(hk + 235);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_238 = buffer.data(hk + 238);
    const auto *hk_239 = buffer.data(hk + 239);
    const auto *hk_240 = buffer.data(hk + 240);
    const auto *hk_241 = buffer.data(hk + 241);
    const auto *hk_242 = buffer.data(hk + 242);
    const auto *hk_243 = buffer.data(hk + 243);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_253 = buffer.data(hk + 253);
    const auto *hk_254 = buffer.data(hk + 254);
    const auto *hk_255 = buffer.data(hk + 255);
    const auto *hk_256 = buffer.data(hk + 256);
    const auto *hk_257 = buffer.data(hk + 257);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_259 = buffer.data(hk + 259);
    const auto *hk_260 = buffer.data(hk + 260);
    const auto *hk_261 = buffer.data(hk + 261);
    const auto *hk_262 = buffer.data(hk + 262);
    const auto *hk_263 = buffer.data(hk + 263);
    const auto *hk_264 = buffer.data(hk + 264);
    const auto *hk_265 = buffer.data(hk + 265);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_267 = buffer.data(hk + 267);
    const auto *hk_268 = buffer.data(hk + 268);
    const auto *hk_269 = buffer.data(hk + 269);
    const auto *hk_270 = buffer.data(hk + 270);
    const auto *hk_271 = buffer.data(hk + 271);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_273 = buffer.data(hk + 273);
    const auto *hk_274 = buffer.data(hk + 274);
    const auto *hk_275 = buffer.data(hk + 275);
    const auto *hk_276 = buffer.data(hk + 276);
    const auto *hk_277 = buffer.data(hk + 277);
    const auto *hk_278 = buffer.data(hk + 278);
    const auto *hk_279 = buffer.data(hk + 279);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_285 = buffer.data(hk + 285);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_289 = buffer.data(hk + 289);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);
    const auto *hk_292 = buffer.data(hk + 292);
    const auto *hk_293 = buffer.data(hk + 293);
    const auto *hk_294 = buffer.data(hk + 294);
    const auto *hk_295 = buffer.data(hk + 295);
    const auto *hk_296 = buffer.data(hk + 296);
    const auto *hk_297 = buffer.data(hk + 297);
    const auto *hk_298 = buffer.data(hk + 298);
    const auto *hk_299 = buffer.data(hk + 299);
    const auto *hk_300 = buffer.data(hk + 300);
    const auto *hk_301 = buffer.data(hk + 301);
    const auto *hk_302 = buffer.data(hk + 302);
    const auto *hk_303 = buffer.data(hk + 303);
    const auto *hk_304 = buffer.data(hk + 304);
    const auto *hk_305 = buffer.data(hk + 305);
    const auto *hk_306 = buffer.data(hk + 306);
    const auto *hk_307 = buffer.data(hk + 307);
    const auto *hk_308 = buffer.data(hk + 308);
    const auto *hk_309 = buffer.data(hk + 309);
    const auto *hk_310 = buffer.data(hk + 310);
    const auto *hk_311 = buffer.data(hk + 311);
    const auto *hk_312 = buffer.data(hk + 312);
    const auto *hk_313 = buffer.data(hk + 313);
    const auto *hk_314 = buffer.data(hk + 314);
    const auto *hk_315 = buffer.data(hk + 315);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_317 = buffer.data(hk + 317);
    const auto *hk_318 = buffer.data(hk + 318);
    const auto *hk_319 = buffer.data(hk + 319);
    const auto *hk_320 = buffer.data(hk + 320);
    const auto *hk_321 = buffer.data(hk + 321);
    const auto *hk_322 = buffer.data(hk + 322);
    const auto *hk_323 = buffer.data(hk + 323);
    const auto *hk_324 = buffer.data(hk + 324);
    const auto *hk_325 = buffer.data(hk + 325);
    const auto *hk_326 = buffer.data(hk + 326);
    const auto *hk_327 = buffer.data(hk + 327);
    const auto *hk_328 = buffer.data(hk + 328);
    const auto *hk_329 = buffer.data(hk + 329);
    const auto *hk_330 = buffer.data(hk + 330);
    const auto *hk_331 = buffer.data(hk + 331);
    const auto *hk_332 = buffer.data(hk + 332);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_334 = buffer.data(hk + 334);
    const auto *hk_335 = buffer.data(hk + 335);
    const auto *hk_336 = buffer.data(hk + 336);
    const auto *hk_337 = buffer.data(hk + 337);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_339 = buffer.data(hk + 339);
    const auto *hk_340 = buffer.data(hk + 340);
    const auto *hk_341 = buffer.data(hk + 341);
    const auto *hk_342 = buffer.data(hk + 342);
    const auto *hk_343 = buffer.data(hk + 343);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_345 = buffer.data(hk + 345);
    const auto *hk_346 = buffer.data(hk + 346);
    const auto *hk_347 = buffer.data(hk + 347);
    const auto *hk_348 = buffer.data(hk + 348);
    const auto *hk_349 = buffer.data(hk + 349);
    const auto *hk_350 = buffer.data(hk + 350);
    const auto *hk_351 = buffer.data(hk + 351);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_353 = buffer.data(hk + 353);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_358 = buffer.data(hk + 358);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_360 = buffer.data(hk + 360);
    const auto *hk_361 = buffer.data(hk + 361);
    const auto *hk_362 = buffer.data(hk + 362);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_364 = buffer.data(hk + 364);
    const auto *hk_365 = buffer.data(hk + 365);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_367 = buffer.data(hk + 367);
    const auto *hk_368 = buffer.data(hk + 368);
    const auto *hk_369 = buffer.data(hk + 369);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_371 = buffer.data(hk + 371);
    const auto *hk_372 = buffer.data(hk + 372);
    const auto *hk_373 = buffer.data(hk + 373);
    const auto *hk_374 = buffer.data(hk + 374);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_376 = buffer.data(hk + 376);
    const auto *hk_377 = buffer.data(hk + 377);
    const auto *hk_378 = buffer.data(hk + 378);
    const auto *hk_379 = buffer.data(hk + 379);
    const auto *hk_380 = buffer.data(hk + 380);
    const auto *hk_381 = buffer.data(hk + 381);
    const auto *hk_382 = buffer.data(hk + 382);
    const auto *hk_383 = buffer.data(hk + 383);
    const auto *hk_384 = buffer.data(hk + 384);
    const auto *hk_385 = buffer.data(hk + 385);
    const auto *hk_386 = buffer.data(hk + 386);
    const auto *hk_387 = buffer.data(hk + 387);
    const auto *hk_388 = buffer.data(hk + 388);
    const auto *hk_389 = buffer.data(hk + 389);
    const auto *hk_390 = buffer.data(hk + 390);
    const auto *hk_391 = buffer.data(hk + 391);
    const auto *hk_392 = buffer.data(hk + 392);
    const auto *hk_393 = buffer.data(hk + 393);
    const auto *hk_394 = buffer.data(hk + 394);
    const auto *hk_395 = buffer.data(hk + 395);
    const auto *hk_396 = buffer.data(hk + 396);
    const auto *hk_397 = buffer.data(hk + 397);
    const auto *hk_398 = buffer.data(hk + 398);
    const auto *hk_399 = buffer.data(hk + 399);
    const auto *hk_400 = buffer.data(hk + 400);
    const auto *hk_401 = buffer.data(hk + 401);
    const auto *hk_402 = buffer.data(hk + 402);
    const auto *hk_403 = buffer.data(hk + 403);
    const auto *hk_404 = buffer.data(hk + 404);
    const auto *hk_405 = buffer.data(hk + 405);
    const auto *hk_406 = buffer.data(hk + 406);
    const auto *hk_407 = buffer.data(hk + 407);
    const auto *hk_408 = buffer.data(hk + 408);
    const auto *hk_409 = buffer.data(hk + 409);
    const auto *hk_410 = buffer.data(hk + 410);
    const auto *hk_411 = buffer.data(hk + 411);
    const auto *hk_412 = buffer.data(hk + 412);
    const auto *hk_413 = buffer.data(hk + 413);
    const auto *hk_414 = buffer.data(hk + 414);
    const auto *hk_415 = buffer.data(hk + 415);
    const auto *hk_416 = buffer.data(hk + 416);
    const auto *hk_417 = buffer.data(hk + 417);
    const auto *hk_418 = buffer.data(hk + 418);
    const auto *hk_419 = buffer.data(hk + 419);
    const auto *hk_420 = buffer.data(hk + 420);
    const auto *hk_421 = buffer.data(hk + 421);
    const auto *hk_422 = buffer.data(hk + 422);
    const auto *hk_423 = buffer.data(hk + 423);
    const auto *hk_424 = buffer.data(hk + 424);
    const auto *hk_425 = buffer.data(hk + 425);
    const auto *hk_426 = buffer.data(hk + 426);
    const auto *hk_427 = buffer.data(hk + 427);
    const auto *hk_428 = buffer.data(hk + 428);
    const auto *hk_429 = buffer.data(hk + 429);
    const auto *hk_430 = buffer.data(hk + 430);
    const auto *hk_431 = buffer.data(hk + 431);
    const auto *hk_432 = buffer.data(hk + 432);
    const auto *hk_433 = buffer.data(hk + 433);
    const auto *hk_434 = buffer.data(hk + 434);
    const auto *hk_435 = buffer.data(hk + 435);
    const auto *hk_436 = buffer.data(hk + 436);
    const auto *hk_437 = buffer.data(hk + 437);
    const auto *hk_438 = buffer.data(hk + 438);
    const auto *hk_439 = buffer.data(hk + 439);
    const auto *hk_440 = buffer.data(hk + 440);
    const auto *hk_441 = buffer.data(hk + 441);
    const auto *hk_442 = buffer.data(hk + 442);
    const auto *hk_443 = buffer.data(hk + 443);
    const auto *hk_444 = buffer.data(hk + 444);
    const auto *hk_445 = buffer.data(hk + 445);
    const auto *hk_446 = buffer.data(hk + 446);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gk_0, hi0_0, hi1_0, \
                         hk_0, hk_1, hk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = pb_y[k] * hk_0[k];

        t_2[k] = pb_z[k] * hk_0[k];

        t_3[k] = f_3 * hi0_0[k]
                 - f_4 * hi1_0[k]
                 + pb_y[k] * hk_1[k];

        t_4[k] = pb_y[k] * hk_2[k];

        t_5[k] = f_3 * hi0_0[k]
                 - f_4 * hi1_0[k]
                 + pb_z[k] * hk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, hi0_1, hi0_2, hi0_3, hi1_1, \
                         hi1_2, hi1_3, hk_3, hk_4, hk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * hi0_1[k]
                 - f_6 * hi1_1[k]
                 + pb_y[k] * hk_3[k];

        t_7[k] = pb_z[k] * hk_3[k];

        t_8[k] = pb_y[k] * hk_4[k];

        t_9[k] = f_5 * hi0_2[k]
                 - f_6 * hi1_2[k]
                 + pb_z[k] * hk_4[k];

        t_10[k] = f_7 * hi0_3[k]
                  - f_8 * hi1_3[k]
                  + pb_y[k] * hk_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, hi0_4, hi0_5, hi1_4, \
                         hi1_5, hk_5, hk_6, hk_7, hk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * hk_5[k];

        t_12[k] = f_3 * hi0_4[k]
                  - f_4 * hi1_4[k]
                  + pb_y[k] * hk_6[k];

        t_13[k] = pb_y[k] * hk_7[k];

        t_14[k] = f_7 * hi0_4[k]
                  - f_8 * hi1_4[k]
                  + pb_z[k] * hk_7[k];

        t_15[k] = f_9 * hi0_5[k]
                  - f_10 * hi1_5[k]
                  + pb_y[k] * hk_8[k];

        t_16[k] = pb_z[k] * hk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, hi0_6, hi0_7, hi1_6, hi1_7, hk_9, \
                         hk_10, hk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * hi0_6[k]
                  - f_6 * hi1_6[k]
                  + pb_y[k] * hk_9[k];

        t_18[k] = f_3 * hi0_7[k]
                  - f_4 * hi1_7[k]
                  + pb_y[k] * hk_10[k];

        t_19[k] = pb_y[k] * hk_11[k];

        t_20[k] = f_9 * hi0_7[k]
                  - f_10 * hi1_7[k]
                  + pb_z[k] * hk_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, hi0_8, hi0_9, hi0_10, hi1_8, \
                         hi1_9, hi1_10, hk_12, hk_13, hk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * hi0_8[k]
                  - f_12 * hi1_8[k]
                  + pb_y[k] * hk_12[k];

        t_22[k] = pb_z[k] * hk_12[k];

        t_23[k] = f_7 * hi0_9[k]
                  - f_8 * hi1_9[k]
                  + pb_y[k] * hk_13[k];

        t_24[k] = f_5 * hi0_10[k]
                  - f_6 * hi1_10[k]
                  + pb_y[k] * hk_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, gk_20, hi0_11, \
                         hi1_11, hk_15, hk_16, hk_17, hk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * hi0_11[k]
                  - f_4 * hi1_11[k]
                  + pb_y[k] * hk_15[k];

        t_26[k] = pb_y[k] * hk_16[k];

        t_27[k] = f_11 * hi0_11[k]
                  - f_12 * hi1_11[k]
                  + pb_z[k] * hk_16[k];

        t_28[k] = f_0 * gk_20[k]
                  + pb_x[k] * hk_19[k];

        t_29[k] = pb_z[k] * hk_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, gk_22, gk_23, gk_24, gk_25, \
                         hk_18, hk_20, hk_21, hk_22, hk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * gk_22[k]
                  + pb_x[k] * hk_20[k];

        t_31[k] = f_0 * gk_23[k]
                  + pb_x[k] * hk_21[k];

        t_32[k] = f_0 * gk_24[k]
                  + pb_x[k] * hk_22[k];

        t_33[k] = f_0 * gk_25[k]
                  + pb_x[k] * hk_23[k];

        t_34[k] = pb_y[k] * hk_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, gk_27, hi0_12, hi0_13, \
                         hi1_12, hi1_13, hk_19, hk_20, hk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * gk_27[k]
                  + pb_x[k] * hk_25[k];

        t_36[k] = f_1 * hi0_12[k]
                  - f_2 * hi1_12[k]
                  + pb_y[k] * hk_19[k];

        t_37[k] = pb_z[k] * hk_19[k];

        t_38[k] = f_11 * hi0_13[k]
                  - f_12 * hi1_13[k]
                  + pb_y[k] * hk_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, hi0_14, hi0_15, hi0_16, hi1_14, hi1_15, \
                         hi1_16, hk_21, hk_22, hk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * hi0_14[k]
                  - f_10 * hi1_14[k]
                  + pb_y[k] * hk_21[k];

        t_40[k] = f_7 * hi0_15[k]
                  - f_8 * hi1_15[k]
                  + pb_y[k] * hk_22[k];

        t_41[k] = f_5 * hi0_16[k]
                  - f_6 * hi1_16[k]
                  + pb_y[k] * hk_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, gk_0, gl_0, \
                         hi0_17, hi1_17, hk_24, hk_25, hk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hi0_17[k]
                  - f_4 * hi1_17[k]
                  + pb_y[k] * hk_24[k];

        t_43[k] = pb_y[k] * hk_25[k];

        t_44[k] = f_1 * hi0_17[k]
                  - f_2 * hi1_17[k]
                  + pb_z[k] * hk_25[k];

        t_45[k] = pa_y[k] * gl_0[k];

        t_46[k] = f_13 * gk_0[k]
                  + pb_y[k] * hk_26[k];

        t_47[k] = pb_z[k] * hk_26[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_z, gk_1, gk_3, gl_1, gl_2, \
                         gl_3, hk_27, hk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * gk_1[k]
                  + pa_y[k] * gl_1[k];

        t_49[k] = pb_z[k] * hk_27[k];

        t_50[k] = pa_y[k] * gl_2[k];

        t_51[k] = f_15 * gk_3[k]
                  + pa_y[k] * gl_3[k];

        t_52[k] = pb_z[k] * hk_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_y, pb_z, gk_4, gk_5, gk_7, \
                         gl_4, gl_5, gl_6, hk_29, hk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * gk_4[k]
                  + pb_y[k] * hk_29[k];

        t_54[k] = pa_y[k] * gl_4[k];

        t_55[k] = f_16 * gk_5[k]
                  + pa_y[k] * gl_5[k];

        t_56[k] = pb_z[k] * hk_30[k];

        t_57[k] = f_14 * gk_7[k]
                  + pa_y[k] * gl_6[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, gk_8, gk_9, gk_11, \
                         gl_7, gl_8, gl_9, hk_31, hk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_13 * gk_8[k]
                  + pb_y[k] * hk_31[k];

        t_59[k] = pa_y[k] * gl_7[k];

        t_60[k] = f_0 * gk_9[k]
                  + pa_y[k] * gl_8[k];

        t_61[k] = pb_z[k] * hk_32[k];

        t_62[k] = f_15 * gk_11[k]
                  + pa_y[k] * gl_9[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, gk_12, gk_13, gk_14, \
                         gl_10, gl_11, gl_12, hk_33, hk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * gk_12[k]
                  + pa_y[k] * gl_10[k];

        t_64[k] = f_13 * gk_13[k]
                  + pb_y[k] * hk_33[k];

        t_65[k] = pa_y[k] * gl_11[k];

        t_66[k] = f_17 * gk_14[k]
                  + pa_y[k] * gl_12[k];

        t_67[k] = pb_z[k] * hk_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pb_y, gk_16, gk_17, gk_18, gk_19, \
                         gl_13, gl_14, gl_15, gl_16, hk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_16 * gk_16[k]
                  + pa_y[k] * gl_13[k];

        t_69[k] = f_15 * gk_17[k]
                  + pa_y[k] * gl_14[k];

        t_70[k] = f_14 * gk_18[k]
                  + pa_y[k] * gl_15[k];

        t_71[k] = f_13 * gk_19[k]
                  + pb_y[k] * hk_35[k];

        t_72[k] = pa_y[k] * gl_16[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, gk_37, gk_38, gk_39, gk_40, \
                         hk_36, hk_37, hk_38, hk_39, hk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_16 * gk_37[k]
                  + pb_x[k] * hk_37[k];

        t_74[k] = pb_z[k] * hk_36[k];

        t_75[k] = f_16 * gk_38[k]
                  + pb_x[k] * hk_38[k];

        t_76[k] = f_16 * gk_39[k]
                  + pb_x[k] * hk_39[k];

        t_77[k] = f_16 * gk_40[k]
                  + pb_x[k] * hk_40[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, gk_20, gk_41, gk_42, \
                         gl_18, gl_19, hk_37, hk_41, hk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * gk_41[k]
                  + pb_x[k] * hk_41[k];

        t_79[k] = f_16 * gk_42[k]
                  + pb_x[k] * hk_42[k];

        t_80[k] = pa_y[k] * gl_18[k];

        t_81[k] = f_18 * gk_20[k]
                  + pa_y[k] * gl_19[k];

        t_82[k] = pb_z[k] * hk_37[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, gk_22, gk_23, gk_24, gk_25, \
                         gk_26, gl_20, gl_21, gl_22, gl_23, gl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_17 * gk_22[k]
                  + pa_y[k] * gl_20[k];

        t_84[k] = f_0 * gk_23[k]
                  + pa_y[k] * gl_21[k];

        t_85[k] = f_16 * gk_24[k]
                  + pa_y[k] * gl_22[k];

        t_86[k] = f_15 * gk_25[k]
                  + pa_y[k] * gl_23[k];

        t_87[k] = f_14 * gk_26[k]
                  + pa_y[k] * gl_24[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, gk_0, gk_27, \
                         gl_0, gl_25, hk_43, hk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_13 * gk_27[k]
                  + pb_y[k] * hk_43[k];

        t_89[k] = pa_y[k] * gl_25[k];

        t_90[k] = pa_z[k] * gl_0[k];

        t_91[k] = pb_y[k] * hk_44[k];

        t_92[k] = f_13 * gk_0[k]
                  + pb_z[k] * hk_44[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, gk_2, gk_3, gl_1, \
                         gl_2, gl_3, hk_45, hk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * gl_1[k];

        t_94[k] = pb_y[k] * hk_45[k];

        t_95[k] = f_14 * gk_2[k]
                  + pa_z[k] * gl_2[k];

        t_96[k] = pa_z[k] * gl_3[k];

        t_97[k] = f_13 * gk_3[k]
                  + pb_z[k] * hk_46[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, gk_4, gk_5, gk_6, \
                         gl_4, gl_5, gl_6, hk_47, hk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * hk_47[k];

        t_99[k] = f_15 * gk_4[k]
                  + pa_z[k] * gl_4[k];

        t_100[k] = pa_z[k] * gl_5[k];

        t_101[k] = f_13 * gk_5[k]
                   + pb_z[k] * hk_48[k];

        t_102[k] = f_14 * gk_6[k]
                   + pa_z[k] * gl_6[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, gk_8, gk_9, \
                         gk_10, gl_7, gl_8, gl_9, hk_49, hk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * hk_49[k];

        t_104[k] = f_16 * gk_8[k]
                   + pa_z[k] * gl_7[k];

        t_105[k] = pa_z[k] * gl_8[k];

        t_106[k] = f_13 * gk_9[k]
                   + pb_z[k] * hk_50[k];

        t_107[k] = f_14 * gk_10[k]
                   + pa_z[k] * gl_9[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, gk_11, gk_13, \
                         gk_14, gl_10, gl_11, gl_12, hk_51, hk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * gk_11[k]
                   + pa_z[k] * gl_10[k];

        t_109[k] = pb_y[k] * hk_51[k];

        t_110[k] = f_0 * gk_13[k]
                   + pa_z[k] * gl_11[k];

        t_111[k] = pa_z[k] * gl_12[k];

        t_112[k] = f_13 * gk_14[k]
                   + pb_z[k] * hk_52[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, gk_15, gk_16, gk_17, \
                         gk_19, gl_13, gl_14, gl_15, gl_16, hk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * gk_15[k]
                   + pa_z[k] * gl_13[k];

        t_114[k] = f_15 * gk_16[k]
                   + pa_z[k] * gl_14[k];

        t_115[k] = f_16 * gk_17[k]
                   + pa_z[k] * gl_15[k];

        t_116[k] = pb_y[k] * hk_53[k];

        t_117[k] = f_17 * gk_19[k]
                   + pa_z[k] * gl_16[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, gk_61, gk_62, gk_63, \
                         gk_64, gl_17, hk_56, hk_57, hk_58, hk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * gl_17[k];

        t_119[k] = f_16 * gk_61[k]
                   + pb_x[k] * hk_56[k];

        t_120[k] = f_16 * gk_62[k]
                   + pb_x[k] * hk_57[k];

        t_121[k] = f_16 * gk_63[k]
                   + pb_x[k] * hk_58[k];

        t_122[k] = f_16 * gk_64[k]
                   + pb_x[k] * hk_59[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, pb_y, gk_65, gk_67, gl_19, \
                         hk_54, hk_60, hk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_16 * gk_65[k]
                   + pb_x[k] * hk_60[k];

        t_124[k] = pb_y[k] * hk_54[k];

        t_125[k] = f_16 * gk_67[k]
                   + pb_x[k] * hk_61[k];

        t_126[k] = pa_z[k] * gl_19[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_z, pb_z, gk_20, gk_21, gk_22, gk_23, \
                         gl_20, gl_21, gl_22, hk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * gk_20[k]
                   + pb_z[k] * hk_55[k];

        t_128[k] = f_14 * gk_21[k]
                   + pa_z[k] * gl_20[k];

        t_129[k] = f_15 * gk_22[k]
                   + pa_z[k] * gl_21[k];

        t_130[k] = f_16 * gk_23[k]
                   + pa_z[k] * gl_22[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_z, pb_y, gk_24, gk_25, gk_27, gl_23, \
                         gl_24, gl_25, hk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * gk_24[k]
                   + pa_z[k] * gl_23[k];

        t_132[k] = f_17 * gk_25[k]
                   + pa_z[k] * gl_24[k];

        t_133[k] = pb_y[k] * hk_61[k];

        t_134[k] = f_18 * gk_27[k]
                   + pa_z[k] * gl_25[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pb_y, pb_z, fl0_0, fl1_0, gk_28, gl_26, \
                         hk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_19 * fl0_0[k]
                   - f_20 * fl1_0[k]
                   + pa_y[k] * gl_26[k];

        t_136[k] = f_14 * gk_28[k]
                   + pb_y[k] * hk_62[k];

        t_137[k] = pb_z[k] * hk_62[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_x, pb_z, gk_70, hi0_18, hi0_20, hi1_18, \
                         hi1_20, hk_63, hk_64, hk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_15 * gk_70[k]
                   + f_11 * hi0_20[k]
                   - f_12 * hi1_20[k]
                   + pb_x[k] * hk_65[k];

        t_139[k] = pb_z[k] * hk_63[k];

        t_140[k] = f_3 * hi0_18[k]
                   - f_4 * hi1_18[k]
                   + pb_z[k] * hk_64[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_x, pb_y, pb_z, gk_30, gk_72, hi0_19, \
                         hi0_22, hi1_19, hi1_22, hk_65, hk_66, hk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_15 * gk_72[k]
                   + f_9 * hi0_22[k]
                   - f_10 * hi1_22[k]
                   + pb_x[k] * hk_67[k];

        t_142[k] = pb_z[k] * hk_65[k];

        t_143[k] = f_14 * gk_30[k]
                   + pb_y[k] * hk_66[k];

        t_144[k] = f_5 * hi0_19[k]
                   - f_6 * hi1_19[k]
                   + pb_z[k] * hk_66[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_z, gk_75, hi0_20, hi0_25, hi1_20, \
                         hi1_25, hk_67, hk_68, hk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_15 * gk_75[k]
                   + f_7 * hi0_25[k]
                   - f_8 * hi1_25[k]
                   + pb_x[k] * hk_70[k];

        t_146[k] = pb_z[k] * hk_67[k];

        t_147[k] = f_3 * hi0_20[k]
                   - f_4 * hi1_20[k]
                   + pb_z[k] * hk_68[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pb_y, pb_z, gk_32, gk_79, hi0_21, \
                         hi0_29, hi1_21, hi1_29, hk_69, hk_70, hk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * gk_32[k]
                   + pb_y[k] * hk_69[k];

        t_149[k] = f_7 * hi0_21[k]
                   - f_8 * hi1_21[k]
                   + pb_z[k] * hk_69[k];

        t_150[k] = f_15 * gk_79[k]
                   + f_5 * hi0_29[k]
                   - f_6 * hi1_29[k]
                   + pb_x[k] * hk_74[k];

        t_151[k] = pb_z[k] * hk_70[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, gk_34, hi0_22, hi0_23, \
                         hi0_24, hi1_22, hi1_23, hi1_24, hk_71, hk_72, \
                         hk_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * hi0_22[k]
                   - f_4 * hi1_22[k]
                   + pb_z[k] * hk_71[k];

        t_153[k] = f_5 * hi0_23[k]
                   - f_6 * hi1_23[k]
                   + pb_z[k] * hk_72[k];

        t_154[k] = f_14 * gk_34[k]
                   + pb_y[k] * hk_73[k];

        t_155[k] = f_9 * hi0_24[k]
                   - f_10 * hi1_24[k]
                   + pb_z[k] * hk_73[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_z, gk_84, hi0_25, hi0_30, hi1_25, \
                         hi1_30, hk_74, hk_75, hk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_15 * gk_84[k]
                   + f_3 * hi0_30[k]
                   - f_4 * hi1_30[k]
                   + pb_x[k] * hk_79[k];

        t_157[k] = pb_z[k] * hk_74[k];

        t_158[k] = f_3 * hi0_25[k]
                   - f_4 * hi1_25[k]
                   + pb_z[k] * hk_75[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_y, pb_z, gk_36, hi0_26, hi0_27, \
                         hi0_28, hi1_26, hi1_27, hi1_28, hk_76, hk_77, \
                         hk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * hi0_26[k]
                   - f_6 * hi1_26[k]
                   + pb_z[k] * hk_76[k];

        t_160[k] = f_7 * hi0_27[k]
                   - f_8 * hi1_27[k]
                   + pb_z[k] * hk_77[k];

        t_161[k] = f_14 * gk_36[k]
                   + pb_y[k] * hk_78[k];

        t_162[k] = f_11 * hi0_28[k]
                   - f_12 * hi1_28[k]
                   + pb_z[k] * hk_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_x, pb_z, gk_85, gk_87, gk_88, \
                         gk_89, hk_79, hk_80, hk_82, hk_83, hk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_15 * gk_85[k]
                   + pb_x[k] * hk_80[k];

        t_164[k] = pb_z[k] * hk_79[k];

        t_165[k] = f_15 * gk_87[k]
                   + pb_x[k] * hk_82[k];

        t_166[k] = f_15 * gk_88[k]
                   + pb_x[k] * hk_83[k];

        t_167[k] = f_15 * gk_89[k]
                   + pb_x[k] * hk_84[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pb_x, fl0_3, fl1_3, gk_90, gk_91, \
                         gk_92, gl_74, hk_85, hk_86, hk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_15 * gk_90[k]
                   + pb_x[k] * hk_85[k];

        t_169[k] = f_15 * gk_91[k]
                   + pb_x[k] * hk_86[k];

        t_170[k] = f_15 * gk_92[k]
                   + pb_x[k] * hk_87[k];

        t_171[k] = f_21 * fl0_3[k]
                   - f_22 * fl1_3[k]
                   + pa_x[k] * gl_74[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, hi0_30, hi0_31, hi0_32, hi1_30, \
                         hi1_31, hi1_32, hk_80, hk_81, hk_82, hk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * hk_80[k];

        t_173[k] = f_3 * hi0_30[k]
                   - f_4 * hi1_30[k]
                   + pb_z[k] * hk_81[k];

        t_174[k] = f_5 * hi0_31[k]
                   - f_6 * hi1_31[k]
                   + pb_z[k] * hk_82[k];

        t_175[k] = f_7 * hi0_32[k]
                   - f_8 * hi1_32[k]
                   + pb_z[k] * hk_83[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pb_z, gk_43, hi0_33, hi0_34, \
                         hi0_35, hi1_33, hi1_34, hi1_35, hk_84, hk_85, \
                         hk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * hi0_33[k]
                   - f_10 * hi1_33[k]
                   + pb_z[k] * hk_84[k];

        t_177[k] = f_11 * hi0_34[k]
                   - f_12 * hi1_34[k]
                   + pb_z[k] * hk_85[k];

        t_178[k] = f_14 * gk_43[k]
                   + pb_y[k] * hk_87[k];

        t_179[k] = f_1 * hi0_35[k]
                   - f_2 * hi1_35[k]
                   + pb_z[k] * hk_87[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pa_y, pa_z, pb_y, gk_45, \
                         gl_27, gl_28, gl_35, gl_36, gl_37, hk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * gl_35[k];

        t_181[k] = pa_z[k] * gl_27[k];

        t_182[k] = pa_y[k] * gl_36[k];

        t_183[k] = pa_z[k] * gl_28[k];

        t_184[k] = f_13 * gk_45[k]
                   + pb_y[k] * hk_88[k];

        t_185[k] = pa_y[k] * gl_37[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, gk_29, \
                         gk_47, gl_29, gl_30, gl_38, hk_89, hk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * gl_29[k];

        t_187[k] = f_13 * gk_29[k]
                   + pb_z[k] * hk_89[k];

        t_188[k] = f_13 * gk_47[k]
                   + pb_y[k] * hk_90[k];

        t_189[k] = pa_y[k] * gl_38[k];

        t_190[k] = pa_z[k] * gl_30[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, gk_31, gk_49, gk_50, \
                         gl_39, gl_40, hk_91, hk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * gk_31[k]
                   + pb_z[k] * hk_91[k];

        t_192[k] = f_14 * gk_49[k]
                   + pa_y[k] * gl_39[k];

        t_193[k] = f_13 * gk_50[k]
                   + pb_y[k] * hk_92[k];

        t_194[k] = pa_y[k] * gl_40[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pa_z, pb_z, gk_33, gk_52, gk_53, \
                         gl_31, gl_41, gl_42, hk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * gl_31[k];

        t_196[k] = f_13 * gk_33[k]
                   + pb_z[k] * hk_93[k];

        t_197[k] = f_15 * gk_52[k]
                   + pa_y[k] * gl_41[k];

        t_198[k] = f_14 * gk_53[k]
                   + pa_y[k] * gl_42[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, gk_35, gk_54, \
                         gl_32, gl_43, hk_94, hk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * gk_54[k]
                   + pb_y[k] * hk_94[k];

        t_200[k] = pa_y[k] * gl_43[k];

        t_201[k] = pa_z[k] * gl_32[k];

        t_202[k] = f_13 * gk_35[k]
                   + pb_z[k] * hk_95[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pa_y, pb_y, gk_56, gk_57, gk_58, \
                         gk_59, gl_44, gl_45, gl_46, gl_47, hk_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_16 * gk_56[k]
                   + pa_y[k] * gl_44[k];

        t_204[k] = f_15 * gk_57[k]
                   + pa_y[k] * gl_45[k];

        t_205[k] = f_14 * gk_58[k]
                   + pa_y[k] * gl_46[k];

        t_206[k] = f_13 * gk_59[k]
                   + pb_y[k] * hk_96[k];

        t_207[k] = pa_y[k] * gl_47[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pa_z, pb_x, gk_103, gk_104, \
                         gk_105, gk_106, gl_33, hk_98, hk_99, hk_100, \
                         hk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_z[k] * gl_33[k];

        t_209[k] = f_15 * gk_103[k]
                   + pb_x[k] * hk_98[k];

        t_210[k] = f_15 * gk_104[k]
                   + pb_x[k] * hk_99[k];

        t_211[k] = f_15 * gk_105[k]
                   + pb_x[k] * hk_100[k];

        t_212[k] = f_15 * gk_106[k]
                   + pb_x[k] * hk_101[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pa_z, pb_x, gk_107, gk_108, gl_34, \
                         gl_48, hk_102, hk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_15 * gk_107[k]
                   + pb_x[k] * hk_102[k];

        t_214[k] = f_15 * gk_108[k]
                   + pb_x[k] * hk_103[k];

        t_215[k] = pa_y[k] * gl_48[k];

        t_216[k] = pa_z[k] * gl_34[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_y, pb_z, gk_37, gk_62, gk_63, gk_64, \
                         gl_49, gl_50, gl_51, hk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_13 * gk_37[k]
                   + pb_z[k] * hk_97[k];

        t_218[k] = f_17 * gk_62[k]
                   + pa_y[k] * gl_49[k];

        t_219[k] = f_0 * gk_63[k]
                   + pa_y[k] * gl_50[k];

        t_220[k] = f_16 * gk_64[k]
                   + pa_y[k] * gl_51[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pb_y, gk_65, gk_66, gk_67, gl_52, \
                         gl_53, gl_54, hk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_15 * gk_65[k]
                   + pa_y[k] * gl_52[k];

        t_222[k] = f_14 * gk_66[k]
                   + pa_y[k] * gl_53[k];

        t_223[k] = f_13 * gk_67[k]
                   + pb_y[k] * hk_104[k];

        t_224[k] = pa_y[k] * gl_54[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_z, pb_y, pb_z, fl0_0, fl1_0, gk_44, \
                         gl_35, hi0_36, hi1_36, hk_105, hk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_19 * fl0_0[k]
                   - f_20 * fl1_0[k]
                   + pa_z[k] * gl_35[k];

        t_226[k] = pb_y[k] * hk_105[k];

        t_227[k] = f_14 * gk_44[k]
                   + pb_z[k] * hk_105[k];

        t_228[k] = f_3 * hi0_36[k]
                   - f_4 * hi1_36[k]
                   + pb_y[k] * hk_106[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, gk_46, gk_114, hi0_37, \
                         hi0_39, hi1_37, hi1_39, hk_107, hk_108, \
                         hk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * hk_107[k];

        t_230[k] = f_15 * gk_114[k]
                   + f_11 * hi0_39[k]
                   - f_12 * hi1_39[k]
                   + pb_x[k] * hk_109[k];

        t_231[k] = f_5 * hi0_37[k]
                   - f_6 * hi1_37[k]
                   + pb_y[k] * hk_108[k];

        t_232[k] = f_14 * gk_46[k]
                   + pb_z[k] * hk_108[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pb_y, pb_z, gk_48, gk_117, hi0_38, \
                         hi0_42, hi1_38, hi1_42, hk_109, hk_110, \
                         hk_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * hk_109[k];

        t_234[k] = f_15 * gk_117[k]
                   + f_9 * hi0_42[k]
                   - f_10 * hi1_42[k]
                   + pb_x[k] * hk_112[k];

        t_235[k] = f_7 * hi0_38[k]
                   - f_8 * hi1_38[k]
                   + pb_y[k] * hk_110[k];

        t_236[k] = f_14 * gk_48[k]
                   + pb_z[k] * hk_110[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_x, pb_y, gk_121, hi0_39, hi0_46, hi1_39, \
                         hi1_46, hk_111, hk_112, hk_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_3 * hi0_39[k]
                   - f_4 * hi1_39[k]
                   + pb_y[k] * hk_111[k];

        t_238[k] = pb_y[k] * hk_112[k];

        t_239[k] = f_15 * gk_121[k]
                   + f_7 * hi0_46[k]
                   - f_8 * hi1_46[k]
                   + pb_x[k] * hk_116[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_y, pb_z, gk_51, hi0_40, hi0_41, \
                         hi0_42, hi1_40, hi1_41, hi1_42, hk_113, hk_114, \
                         hk_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_9 * hi0_40[k]
                   - f_10 * hi1_40[k]
                   + pb_y[k] * hk_113[k];

        t_241[k] = f_14 * gk_51[k]
                   + pb_z[k] * hk_113[k];

        t_242[k] = f_5 * hi0_41[k]
                   - f_6 * hi1_41[k]
                   + pb_y[k] * hk_114[k];

        t_243[k] = f_3 * hi0_42[k]
                   - f_4 * hi1_42[k]
                   + pb_y[k] * hk_115[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pb_y, pb_z, gk_55, gk_126, hi0_43, \
                         hi0_47, hi1_43, hi1_47, hk_116, hk_117, \
                         hk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_y[k] * hk_116[k];

        t_245[k] = f_15 * gk_126[k]
                   + f_5 * hi0_47[k]
                   - f_6 * hi1_47[k]
                   + pb_x[k] * hk_121[k];

        t_246[k] = f_11 * hi0_43[k]
                   - f_12 * hi1_43[k]
                   + pb_y[k] * hk_117[k];

        t_247[k] = f_14 * gk_55[k]
                   + pb_z[k] * hk_117[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_y, hi0_44, hi0_45, hi0_46, hi1_44, \
                         hi1_45, hi1_46, hk_118, hk_119, hk_120, \
                         hk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * hi0_44[k]
                   - f_8 * hi1_44[k]
                   + pb_y[k] * hk_118[k];

        t_249[k] = f_5 * hi0_45[k]
                   - f_6 * hi1_45[k]
                   + pb_y[k] * hk_119[k];

        t_250[k] = f_3 * hi0_46[k]
                   - f_4 * hi1_46[k]
                   + pb_y[k] * hk_120[k];

        t_251[k] = pb_y[k] * hk_121[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, gk_127, gk_128, gk_129, gk_130, \
                         hi0_53, hi1_53, hk_122, hk_123, hk_124, \
                         hk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_15 * gk_127[k]
                   + f_3 * hi0_53[k]
                   - f_4 * hi1_53[k]
                   + pb_x[k] * hk_122[k];

        t_253[k] = f_15 * gk_128[k]
                   + pb_x[k] * hk_123[k];

        t_254[k] = f_15 * gk_129[k]
                   + pb_x[k] * hk_124[k];

        t_255[k] = f_15 * gk_130[k]
                   + pb_x[k] * hk_125[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, gk_131, gk_132, \
                         gk_133, gk_135, hk_122, hk_126, hk_127, hk_128, \
                         hk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_15 * gk_131[k]
                   + pb_x[k] * hk_126[k];

        t_257[k] = f_15 * gk_132[k]
                   + pb_x[k] * hk_127[k];

        t_258[k] = f_15 * gk_133[k]
                   + pb_x[k] * hk_128[k];

        t_259[k] = pb_y[k] * hk_122[k];

        t_260[k] = f_15 * gk_135[k]
                   + pb_x[k] * hk_130[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_y, pb_z, gk_60, hi0_48, hi0_49, \
                         hi0_50, hi1_48, hi1_49, hi1_50, hk_123, hk_125, \
                         hk_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * hi0_48[k]
                   - f_2 * hi1_48[k]
                   + pb_y[k] * hk_123[k];

        t_262[k] = f_14 * gk_60[k]
                   + pb_z[k] * hk_123[k];

        t_263[k] = f_11 * hi0_49[k]
                   - f_12 * hi1_49[k]
                   + pb_y[k] * hk_125[k];

        t_264[k] = f_9 * hi0_50[k]
                   - f_10 * hi1_50[k]
                   + pb_y[k] * hk_126[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_y, hi0_51, hi0_52, hi0_53, hi1_51, \
                         hi1_52, hi1_53, hk_127, hk_128, hk_129, \
                         hk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_7 * hi0_51[k]
                   - f_8 * hi1_51[k]
                   + pb_y[k] * hk_127[k];

        t_266[k] = f_5 * hi0_52[k]
                   - f_6 * hi1_52[k]
                   + pb_y[k] * hk_128[k];

        t_267[k] = f_3 * hi0_53[k]
                   - f_4 * hi1_53[k]
                   + pb_y[k] * hk_129[k];

        t_268[k] = pb_y[k] * hk_130[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pa_y, pb_y, pb_z, fl0_1, fl0_4, \
                         fl1_1, fl1_4, gk_68, gl_55, gl_106, hk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_21 * fl0_4[k]
                   - f_22 * fl1_4[k]
                   + pa_x[k] * gl_106[k];

        t_270[k] = f_21 * fl0_1[k]
                   - f_22 * fl1_1[k]
                   + pa_y[k] * gl_55[k];

        t_271[k] = f_15 * gk_68[k]
                   + pb_y[k] * hk_131[k];

        t_272[k] = pb_z[k] * hk_131[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pb_x, pb_z, gk_137, hi0_54, hi0_56, hi1_54, \
                         hi1_56, hk_132, hk_133, hk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_14 * gk_137[k]
                   + f_11 * hi0_56[k]
                   - f_12 * hi1_56[k]
                   + pb_x[k] * hk_134[k];

        t_274[k] = pb_z[k] * hk_132[k];

        t_275[k] = f_3 * hi0_54[k]
                   - f_4 * hi1_54[k]
                   + pb_z[k] * hk_133[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, pb_y, pb_z, gk_71, gk_139, hi0_55, \
                         hi0_58, hi1_55, hi1_58, hk_134, hk_135, \
                         hk_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_14 * gk_139[k]
                   + f_9 * hi0_58[k]
                   - f_10 * hi1_58[k]
                   + pb_x[k] * hk_136[k];

        t_277[k] = pb_z[k] * hk_134[k];

        t_278[k] = f_15 * gk_71[k]
                   + pb_y[k] * hk_135[k];

        t_279[k] = f_5 * hi0_55[k]
                   - f_6 * hi1_55[k]
                   + pb_z[k] * hk_135[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_x, pb_z, gk_141, hi0_56, hi0_61, hi1_56, \
                         hi1_61, hk_136, hk_137, hk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_14 * gk_141[k]
                   + f_7 * hi0_61[k]
                   - f_8 * hi1_61[k]
                   + pb_x[k] * hk_139[k];

        t_281[k] = pb_z[k] * hk_136[k];

        t_282[k] = f_3 * hi0_56[k]
                   - f_4 * hi1_56[k]
                   + pb_z[k] * hk_137[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_x, pb_y, pb_z, gk_74, gk_143, hi0_57, \
                         hi0_65, hi1_57, hi1_65, hk_138, hk_139, \
                         hk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * gk_74[k]
                   + pb_y[k] * hk_138[k];

        t_284[k] = f_7 * hi0_57[k]
                   - f_8 * hi1_57[k]
                   + pb_z[k] * hk_138[k];

        t_285[k] = f_14 * gk_143[k]
                   + f_5 * hi0_65[k]
                   - f_6 * hi1_65[k]
                   + pb_x[k] * hk_143[k];

        t_286[k] = pb_z[k] * hk_139[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_y, pb_z, gk_78, hi0_58, hi0_59, \
                         hi0_60, hi1_58, hi1_59, hi1_60, hk_140, hk_141, \
                         hk_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * hi0_58[k]
                   - f_4 * hi1_58[k]
                   + pb_z[k] * hk_140[k];

        t_288[k] = f_5 * hi0_59[k]
                   - f_6 * hi1_59[k]
                   + pb_z[k] * hk_141[k];

        t_289[k] = f_15 * gk_78[k]
                   + pb_y[k] * hk_142[k];

        t_290[k] = f_9 * hi0_60[k]
                   - f_10 * hi1_60[k]
                   + pb_z[k] * hk_142[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_x, pb_z, gk_145, hi0_61, hi0_66, hi1_61, \
                         hi1_66, hk_143, hk_144, hk_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_14 * gk_145[k]
                   + f_3 * hi0_66[k]
                   - f_4 * hi1_66[k]
                   + pb_x[k] * hk_148[k];

        t_292[k] = pb_z[k] * hk_143[k];

        t_293[k] = f_3 * hi0_61[k]
                   - f_4 * hi1_61[k]
                   + pb_z[k] * hk_144[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pb_z, gk_83, hi0_62, hi0_63, \
                         hi0_64, hi1_62, hi1_63, hi1_64, hk_145, hk_146, \
                         hk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_5 * hi0_62[k]
                   - f_6 * hi1_62[k]
                   + pb_z[k] * hk_145[k];

        t_295[k] = f_7 * hi0_63[k]
                   - f_8 * hi1_63[k]
                   + pb_z[k] * hk_146[k];

        t_296[k] = f_15 * gk_83[k]
                   + pb_y[k] * hk_147[k];

        t_297[k] = f_11 * hi0_64[k]
                   - f_12 * hi1_64[k]
                   + pb_z[k] * hk_147[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pb_x, pb_z, gk_146, gk_147, \
                         gk_148, gk_149, hk_148, hk_149, hk_151, hk_152, \
                         hk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_14 * gk_146[k]
                   + pb_x[k] * hk_149[k];

        t_299[k] = pb_z[k] * hk_148[k];

        t_300[k] = f_14 * gk_147[k]
                   + pb_x[k] * hk_151[k];

        t_301[k] = f_14 * gk_148[k]
                   + pb_x[k] * hk_152[k];

        t_302[k] = f_14 * gk_149[k]
                   + pb_x[k] * hk_153[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pb_x, fl0_5, fl1_5, gk_150, gk_151, \
                         gk_152, gl_115, hk_154, hk_155, hk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_14 * gk_150[k]
                   + pb_x[k] * hk_154[k];

        t_304[k] = f_14 * gk_151[k]
                   + pb_x[k] * hk_155[k];

        t_305[k] = f_14 * gk_152[k]
                   + pb_x[k] * hk_156[k];

        t_306[k] = f_19 * fl0_5[k]
                   - f_20 * fl1_5[k]
                   + pa_x[k] * gl_115[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pb_z, hi0_66, hi0_67, hi0_68, hi1_66, \
                         hi1_67, hi1_68, hk_149, hk_150, hk_151, \
                         hk_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = pb_z[k] * hk_149[k];

        t_308[k] = f_3 * hi0_66[k]
                   - f_4 * hi1_66[k]
                   + pb_z[k] * hk_150[k];

        t_309[k] = f_5 * hi0_67[k]
                   - f_6 * hi1_67[k]
                   + pb_z[k] * hk_151[k];

        t_310[k] = f_7 * hi0_68[k]
                   - f_8 * hi1_68[k]
                   + pb_z[k] * hk_152[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_y, pb_z, gk_92, hi0_69, hi0_70, \
                         hi0_71, hi1_69, hi1_70, hi1_71, hk_153, hk_154, \
                         hk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_9 * hi0_69[k]
                   - f_10 * hi1_69[k]
                   + pb_z[k] * hk_153[k];

        t_312[k] = f_11 * hi0_70[k]
                   - f_12 * hi1_70[k]
                   + pb_z[k] * hk_154[k];

        t_313[k] = f_15 * gk_92[k]
                   + pb_y[k] * hk_156[k];

        t_314[k] = f_1 * hi0_71[k]
                   - f_2 * hi1_71[k]
                   + pb_z[k] * hk_156[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, pa_z, pb_y, pb_z, gk_68, gk_93, \
                         gl_55, gl_56, gl_57, hk_157, hk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * gl_55[k];

        t_316[k] = pa_z[k] * gl_56[k];

        t_317[k] = f_13 * gk_68[k]
                   + pb_z[k] * hk_157[k];

        t_318[k] = pa_z[k] * gl_57[k];

        t_319[k] = f_14 * gk_93[k]
                   + pb_y[k] * hk_158[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_z, pb_y, pb_z, gk_69, gk_70, gk_95, \
                         gl_58, gl_59, hk_159, hk_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_14 * gk_69[k]
                   + pa_z[k] * gl_58[k];

        t_321[k] = pa_z[k] * gl_59[k];

        t_322[k] = f_13 * gk_70[k]
                   + pb_z[k] * hk_159[k];

        t_323[k] = f_14 * gk_95[k]
                   + pb_y[k] * hk_160[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_z, gk_71, gk_72, gk_73, gl_60, \
                         gl_61, gl_62, hk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_15 * gk_71[k]
                   + pa_z[k] * gl_60[k];

        t_325[k] = pa_z[k] * gl_61[k];

        t_326[k] = f_13 * gk_72[k]
                   + pb_z[k] * hk_161[k];

        t_327[k] = f_14 * gk_73[k]
                   + pa_z[k] * gl_62[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_z, pb_y, pb_z, gk_74, gk_75, gk_97, \
                         gl_63, gl_64, hk_162, hk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_14 * gk_97[k]
                   + pb_y[k] * hk_162[k];

        t_329[k] = f_16 * gk_74[k]
                   + pa_z[k] * gl_63[k];

        t_330[k] = pa_z[k] * gl_64[k];

        t_331[k] = f_13 * gk_75[k]
                   + pb_z[k] * hk_163[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pa_z, pb_y, gk_76, gk_77, gk_78, \
                         gk_99, gl_65, gl_66, gl_67, gl_68, hk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_14 * gk_76[k]
                   + pa_z[k] * gl_65[k];

        t_333[k] = f_15 * gk_77[k]
                   + pa_z[k] * gl_66[k];

        t_334[k] = f_14 * gk_99[k]
                   + pb_y[k] * hk_164[k];

        t_335[k] = f_0 * gk_78[k]
                   + pa_z[k] * gl_67[k];

        t_336[k] = pa_z[k] * gl_68[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pb_z, gk_79, gk_80, gk_81, gk_82, \
                         gl_69, gl_70, gl_71, hk_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_13 * gk_79[k]
                   + pb_z[k] * hk_165[k];

        t_338[k] = f_14 * gk_80[k]
                   + pa_z[k] * gl_69[k];

        t_339[k] = f_15 * gk_81[k]
                   + pa_z[k] * gl_70[k];

        t_340[k] = f_16 * gk_82[k]
                   + pa_z[k] * gl_71[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_z, pb_x, pb_y, gk_83, gk_101, gk_163, \
                         gl_72, gl_73, hk_166, hk_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_14 * gk_101[k]
                   + pb_y[k] * hk_166[k];

        t_342[k] = f_17 * gk_83[k]
                   + pa_z[k] * gl_72[k];

        t_343[k] = pa_z[k] * gl_73[k];

        t_344[k] = f_14 * gk_163[k]
                   + pb_x[k] * hk_168[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pb_x, gk_164, gk_165, gk_166, \
                         gk_167, gk_168, hk_169, hk_170, hk_171, hk_172, \
                         hk_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_14 * gk_164[k]
                   + pb_x[k] * hk_169[k];

        t_346[k] = f_14 * gk_165[k]
                   + pb_x[k] * hk_170[k];

        t_347[k] = f_14 * gk_166[k]
                   + pb_x[k] * hk_171[k];

        t_348[k] = f_14 * gk_167[k]
                   + pb_x[k] * hk_172[k];

        t_349[k] = f_14 * gk_168[k]
                   + pb_x[k] * hk_173[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_z, pb_x, pb_z, gk_85, gk_86, gk_169, \
                         gl_74, gl_75, hk_167, hk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_14 * gk_169[k]
                   + pb_x[k] * hk_174[k];

        t_351[k] = pa_z[k] * gl_74[k];

        t_352[k] = f_13 * gk_85[k]
                   + pb_z[k] * hk_167[k];

        t_353[k] = f_14 * gk_86[k]
                   + pa_z[k] * gl_75[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_z, gk_87, gk_88, gk_89, gk_90, gl_76, \
                         gl_77, gl_78, gl_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_15 * gk_87[k]
                   + pa_z[k] * gl_76[k];

        t_355[k] = f_16 * gk_88[k]
                   + pa_z[k] * gl_77[k];

        t_356[k] = f_0 * gk_89[k]
                   + pa_z[k] * gl_78[k];

        t_357[k] = f_17 * gk_90[k]
                   + pa_z[k] * gl_79[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, pa_y, pa_z, pb_y, gk_92, gk_109, \
                         gk_110, gl_80, gl_81, gl_82, hk_174, hk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_14 * gk_109[k]
                   + pb_y[k] * hk_174[k];

        t_359[k] = f_18 * gk_92[k]
                   + pa_z[k] * gl_80[k];

        t_360[k] = pa_y[k] * gl_81[k];

        t_361[k] = f_13 * gk_110[k]
                   + pb_y[k] * hk_175[k];

        t_362[k] = pa_y[k] * gl_82[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pb_y, gk_111, gk_112, gk_113, \
                         gl_83, gl_84, gl_85, hk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * gk_111[k]
                   + pa_y[k] * gl_83[k];

        t_364[k] = f_13 * gk_112[k]
                   + pb_y[k] * hk_176[k];

        t_365[k] = pa_y[k] * gl_84[k];

        t_366[k] = f_15 * gk_113[k]
                   + pa_y[k] * gl_85[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pb_y, pb_z, gk_94, gk_114, gk_115, \
                         gl_86, gl_87, hk_177, hk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_14 * gk_94[k]
                   + pb_z[k] * hk_177[k];

        t_368[k] = f_13 * gk_114[k]
                   + pb_y[k] * hk_178[k];

        t_369[k] = pa_y[k] * gl_86[k];

        t_370[k] = f_16 * gk_115[k]
                   + pa_y[k] * gl_87[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pb_y, pb_z, gk_96, gk_116, gk_117, \
                         gl_88, gl_89, hk_179, hk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_14 * gk_96[k]
                   + pb_z[k] * hk_179[k];

        t_372[k] = f_14 * gk_116[k]
                   + pa_y[k] * gl_88[k];

        t_373[k] = f_13 * gk_117[k]
                   + pb_y[k] * hk_180[k];

        t_374[k] = pa_y[k] * gl_89[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pa_y, pb_z, gk_98, gk_118, gk_119, \
                         gk_120, gl_90, gl_91, gl_92, hk_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_0 * gk_118[k]
                   + pa_y[k] * gl_90[k];

        t_376[k] = f_14 * gk_98[k]
                   + pb_z[k] * hk_181[k];

        t_377[k] = f_15 * gk_119[k]
                   + pa_y[k] * gl_91[k];

        t_378[k] = f_14 * gk_120[k]
                   + pa_y[k] * gl_92[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pa_y, pb_y, pb_z, gk_100, gk_121, gk_122, \
                         gl_93, gl_94, hk_182, hk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_13 * gk_121[k]
                   + pb_y[k] * hk_182[k];

        t_380[k] = pa_y[k] * gl_93[k];

        t_381[k] = f_17 * gk_122[k]
                   + pa_y[k] * gl_94[k];

        t_382[k] = f_14 * gk_100[k]
                   + pb_z[k] * hk_183[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pa_y, pb_y, gk_123, gk_124, \
                         gk_125, gk_126, gl_95, gl_96, gl_97, gl_98, \
                         hk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_16 * gk_123[k]
                   + pa_y[k] * gl_95[k];

        t_384[k] = f_15 * gk_124[k]
                   + pa_y[k] * gl_96[k];

        t_385[k] = f_14 * gk_125[k]
                   + pa_y[k] * gl_97[k];

        t_386[k] = f_13 * gk_126[k]
                   + pb_y[k] * hk_184[k];

        t_387[k] = pa_y[k] * gl_98[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, pb_x, gk_180, gk_181, gk_182, \
                         gk_183, gk_184, hk_185, hk_186, hk_187, hk_188, \
                         hk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_14 * gk_180[k]
                   + pb_x[k] * hk_185[k];

        t_389[k] = f_14 * gk_181[k]
                   + pb_x[k] * hk_186[k];

        t_390[k] = f_14 * gk_182[k]
                   + pb_x[k] * hk_187[k];

        t_391[k] = f_14 * gk_183[k]
                   + pb_x[k] * hk_188[k];

        t_392[k] = f_14 * gk_184[k]
                   + pb_x[k] * hk_189[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_y, pb_x, gk_128, gk_185, gk_186, \
                         gl_99, gl_100, hk_190, hk_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_14 * gk_185[k]
                   + pb_x[k] * hk_190[k];

        t_394[k] = f_14 * gk_186[k]
                   + pb_x[k] * hk_191[k];

        t_395[k] = pa_y[k] * gl_99[k];

        t_396[k] = f_18 * gk_128[k]
                   + pa_y[k] * gl_100[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pb_z, gk_102, gk_130, gk_131, \
                         gk_132, gl_101, gl_102, gl_103, hk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_14 * gk_102[k]
                   + pb_z[k] * hk_185[k];

        t_398[k] = f_17 * gk_130[k]
                   + pa_y[k] * gl_101[k];

        t_399[k] = f_0 * gk_131[k]
                   + pa_y[k] * gl_102[k];

        t_400[k] = f_16 * gk_132[k]
                   + pa_y[k] * gl_103[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pb_y, gk_133, gk_134, gk_135, \
                         gl_104, gl_105, gl_106, hk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * gk_133[k]
                   + pa_y[k] * gl_104[k];

        t_402[k] = f_14 * gk_134[k]
                   + pa_y[k] * gl_105[k];

        t_403[k] = f_13 * gk_135[k]
                   + pb_y[k] * hk_192[k];

        t_404[k] = pa_y[k] * gl_106[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_z, pb_y, pb_z, fl0_2, fl1_2, gk_110, \
                         gl_81, hi0_72, hi1_72, hk_193, hk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_21 * fl0_2[k]
                   - f_22 * fl1_2[k]
                   + pa_z[k] * gl_81[k];

        t_406[k] = pb_y[k] * hk_193[k];

        t_407[k] = f_15 * gk_110[k]
                   + pb_z[k] * hk_193[k];

        t_408[k] = f_3 * hi0_72[k]
                   - f_4 * hi1_72[k]
                   + pb_y[k] * hk_194[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, pb_x, pb_y, pb_z, gk_113, gk_190, hi0_73, \
                         hi0_75, hi1_73, hi1_75, hk_195, hk_196, \
                         hk_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * hk_195[k];

        t_410[k] = f_14 * gk_190[k]
                   + f_11 * hi0_75[k]
                   - f_12 * hi1_75[k]
                   + pb_x[k] * hk_197[k];

        t_411[k] = f_5 * hi0_73[k]
                   - f_6 * hi1_73[k]
                   + pb_y[k] * hk_196[k];

        t_412[k] = f_15 * gk_113[k]
                   + pb_z[k] * hk_196[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, pb_x, pb_y, pb_z, gk_115, gk_192, hi0_74, \
                         hi0_78, hi1_74, hi1_78, hk_197, hk_198, \
                         hk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_y[k] * hk_197[k];

        t_414[k] = f_14 * gk_192[k]
                   + f_9 * hi0_78[k]
                   - f_10 * hi1_78[k]
                   + pb_x[k] * hk_200[k];

        t_415[k] = f_7 * hi0_74[k]
                   - f_8 * hi1_74[k]
                   + pb_y[k] * hk_198[k];

        t_416[k] = f_15 * gk_115[k]
                   + pb_z[k] * hk_198[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_x, pb_y, gk_194, hi0_75, hi0_82, hi1_75, \
                         hi1_82, hk_199, hk_200, hk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_3 * hi0_75[k]
                   - f_4 * hi1_75[k]
                   + pb_y[k] * hk_199[k];

        t_418[k] = pb_y[k] * hk_200[k];

        t_419[k] = f_14 * gk_194[k]
                   + f_7 * hi0_82[k]
                   - f_8 * hi1_82[k]
                   + pb_x[k] * hk_204[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_y, pb_z, gk_118, hi0_76, hi0_77, \
                         hi0_78, hi1_76, hi1_77, hi1_78, hk_201, hk_202, \
                         hk_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_9 * hi0_76[k]
                   - f_10 * hi1_76[k]
                   + pb_y[k] * hk_201[k];

        t_421[k] = f_15 * gk_118[k]
                   + pb_z[k] * hk_201[k];

        t_422[k] = f_5 * hi0_77[k]
                   - f_6 * hi1_77[k]
                   + pb_y[k] * hk_202[k];

        t_423[k] = f_3 * hi0_78[k]
                   - f_4 * hi1_78[k]
                   + pb_y[k] * hk_203[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pb_x, pb_y, pb_z, gk_122, gk_196, hi0_79, \
                         hi0_83, hi1_79, hi1_83, hk_204, hk_205, \
                         hk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_y[k] * hk_204[k];

        t_425[k] = f_14 * gk_196[k]
                   + f_5 * hi0_83[k]
                   - f_6 * hi1_83[k]
                   + pb_x[k] * hk_209[k];

        t_426[k] = f_11 * hi0_79[k]
                   - f_12 * hi1_79[k]
                   + pb_y[k] * hk_205[k];

        t_427[k] = f_15 * gk_122[k]
                   + pb_z[k] * hk_205[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pb_y, hi0_80, hi0_81, hi0_82, hi1_80, \
                         hi1_81, hi1_82, hk_206, hk_207, hk_208, \
                         hk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_7 * hi0_80[k]
                   - f_8 * hi1_80[k]
                   + pb_y[k] * hk_206[k];

        t_429[k] = f_5 * hi0_81[k]
                   - f_6 * hi1_81[k]
                   + pb_y[k] * hk_207[k];

        t_430[k] = f_3 * hi0_82[k]
                   - f_4 * hi1_82[k]
                   + pb_y[k] * hk_208[k];

        t_431[k] = pb_y[k] * hk_209[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_x, gk_197, gk_198, gk_199, gk_200, \
                         hi0_89, hi1_89, hk_210, hk_211, hk_212, \
                         hk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_14 * gk_197[k]
                   + f_3 * hi0_89[k]
                   - f_4 * hi1_89[k]
                   + pb_x[k] * hk_210[k];

        t_433[k] = f_14 * gk_198[k]
                   + pb_x[k] * hk_211[k];

        t_434[k] = f_14 * gk_199[k]
                   + pb_x[k] * hk_212[k];

        t_435[k] = f_14 * gk_200[k]
                   + pb_x[k] * hk_213[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pb_x, pb_y, gk_201, gk_202, \
                         gk_203, gk_204, hk_210, hk_214, hk_215, hk_216, \
                         hk_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_14 * gk_201[k]
                   + pb_x[k] * hk_214[k];

        t_437[k] = f_14 * gk_202[k]
                   + pb_x[k] * hk_215[k];

        t_438[k] = f_14 * gk_203[k]
                   + pb_x[k] * hk_216[k];

        t_439[k] = pb_y[k] * hk_210[k];

        t_440[k] = f_14 * gk_204[k]
                   + pb_x[k] * hk_218[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_y, pb_z, gk_128, hi0_84, hi0_85, \
                         hi0_86, hi1_84, hi1_85, hi1_86, hk_211, hk_213, \
                         hk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * hi0_84[k]
                   - f_2 * hi1_84[k]
                   + pb_y[k] * hk_211[k];

        t_442[k] = f_15 * gk_128[k]
                   + pb_z[k] * hk_211[k];

        t_443[k] = f_11 * hi0_85[k]
                   - f_12 * hi1_85[k]
                   + pb_y[k] * hk_213[k];

        t_444[k] = f_9 * hi0_86[k]
                   - f_10 * hi1_86[k]
                   + pb_y[k] * hk_214[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pb_y, hi0_87, hi0_88, hi0_89, hi1_87, \
                         hi1_88, hi1_89, hk_215, hk_216, hk_217, \
                         hk_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_7 * hi0_87[k]
                   - f_8 * hi1_87[k]
                   + pb_y[k] * hk_215[k];

        t_446[k] = f_5 * hi0_88[k]
                   - f_6 * hi1_88[k]
                   + pb_y[k] * hk_216[k];

        t_447[k] = f_3 * hi0_89[k]
                   - f_4 * hi1_89[k]
                   + pb_y[k] * hk_217[k];

        t_448[k] = pb_y[k] * hk_218[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_x, pb_y, pb_z, fl0_8, fl1_8, gk_136, \
                         gk_205, gl_124, gl_125, hk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_19 * fl0_8[k]
                   - f_20 * fl1_8[k]
                   + pa_x[k] * gl_124[k];

        t_450[k] = f_18 * gk_205[k]
                   + pa_x[k] * gl_125[k];

        t_451[k] = f_16 * gk_136[k]
                   + pb_y[k] * hk_219[k];

        t_452[k] = pb_z[k] * hk_219[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, pa_x, pb_z, gk_207, gk_208, \
                         gk_209, gl_127, gl_128, gl_129, hk_220, \
                         hk_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_17 * gk_207[k]
                   + pa_x[k] * gl_127[k];

        t_454[k] = pb_z[k] * hk_220[k];

        t_455[k] = f_17 * gk_208[k]
                   + pa_x[k] * gl_128[k];

        t_456[k] = f_0 * gk_209[k]
                   + pa_x[k] * gl_129[k];

        t_457[k] = pb_z[k] * hk_221[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pa_x, pb_y, pb_z, gk_138, gk_211, gk_212, \
                         gl_130, gl_131, hk_222, hk_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_16 * gk_138[k]
                   + pb_y[k] * hk_222[k];

        t_459[k] = f_0 * gk_211[k]
                   + pa_x[k] * gl_130[k];

        t_460[k] = f_16 * gk_212[k]
                   + pa_x[k] * gl_131[k];

        t_461[k] = pb_z[k] * hk_223[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pa_x, pb_y, gk_140, gk_214, gk_215, \
                         gk_216, gl_132, gl_133, gl_134, hk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_16 * gk_214[k]
                   + pa_x[k] * gl_132[k];

        t_463[k] = f_16 * gk_140[k]
                   + pb_y[k] * hk_224[k];

        t_464[k] = f_16 * gk_215[k]
                   + pa_x[k] * gl_133[k];

        t_465[k] = f_15 * gk_216[k]
                   + pa_x[k] * gl_134[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_x, pb_y, pb_z, gk_142, gk_218, gk_219, \
                         gl_135, gl_136, hk_225, hk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = pb_z[k] * hk_225[k];

        t_467[k] = f_15 * gk_218[k]
                   + pa_x[k] * gl_135[k];

        t_468[k] = f_15 * gk_219[k]
                   + pa_x[k] * gl_136[k];

        t_469[k] = f_16 * gk_142[k]
                   + pb_y[k] * hk_226[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pa_x, pb_z, gk_220, gk_221, \
                         gk_222, gk_223, gl_137, gl_138, gl_139, gl_140, \
                         hk_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * gk_220[k]
                   + pa_x[k] * gl_137[k];

        t_471[k] = f_14 * gk_221[k]
                   + pa_x[k] * gl_138[k];

        t_472[k] = pb_z[k] * hk_227[k];

        t_473[k] = f_14 * gk_222[k]
                   + pa_x[k] * gl_139[k];

        t_474[k] = f_14 * gk_223[k]
                   + pa_x[k] * gl_140[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_x, pb_x, pb_y, gk_144, gk_224, gk_225, \
                         gk_226, gl_141, gl_142, hk_228, hk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_14 * gk_224[k]
                   + pa_x[k] * gl_141[k];

        t_476[k] = f_16 * gk_144[k]
                   + pb_y[k] * hk_228[k];

        t_477[k] = f_14 * gk_225[k]
                   + pa_x[k] * gl_142[k];

        t_478[k] = f_13 * gk_226[k]
                   + pb_x[k] * hk_230[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, pb_x, pb_z, gk_228, gk_229, \
                         gk_230, gk_231, hk_229, hk_231, hk_232, hk_233, \
                         hk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pb_z[k] * hk_229[k];

        t_480[k] = f_13 * gk_228[k]
                   + pb_x[k] * hk_231[k];

        t_481[k] = f_13 * gk_229[k]
                   + pb_x[k] * hk_232[k];

        t_482[k] = f_13 * gk_230[k]
                   + pb_x[k] * hk_233[k];

        t_483[k] = f_13 * gk_231[k]
                   + pb_x[k] * hk_234[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, pa_x, pb_x, pb_z, gk_232, gk_233, \
                         gl_143, gl_144, hk_230, hk_235, hk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_13 * gk_232[k]
                   + pb_x[k] * hk_235[k];

        t_485[k] = f_13 * gk_233[k]
                   + pb_x[k] * hk_236[k];

        t_486[k] = pa_x[k] * gl_143[k];

        t_487[k] = pb_z[k] * hk_230[k];

        t_488[k] = pa_x[k] * gl_144[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, t_494, t_495, pa_x, pa_z, gl_107, \
                         gl_145, gl_146, gl_147, gl_148, gl_149, \
                         gl_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = pa_x[k] * gl_145[k];

        t_490[k] = pa_x[k] * gl_146[k];

        t_491[k] = pa_x[k] * gl_147[k];

        t_492[k] = pa_x[k] * gl_148[k];

        t_493[k] = pa_x[k] * gl_149[k];

        t_494[k] = pa_x[k] * gl_150[k];

        t_495[k] = pa_z[k] * gl_107[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_z, pb_y, pb_z, gk_136, gk_154, gl_108, \
                         gl_109, hk_237, hk_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = pa_z[k] * gl_108[k];

        t_497[k] = f_13 * gk_136[k]
                   + pb_z[k] * hk_237[k];

        t_498[k] = pa_z[k] * gl_109[k];

        t_499[k] = f_15 * gk_154[k]
                   + pb_y[k] * hk_238[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_x, pa_z, pb_y, pb_z, gk_137, gk_156, \
                         gk_237, gl_110, gl_151, hk_239, hk_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_17 * gk_237[k]
                   + pa_x[k] * gl_151[k];

        t_501[k] = pa_z[k] * gl_110[k];

        t_502[k] = f_13 * gk_137[k]
                   + pb_z[k] * hk_239[k];

        t_503[k] = f_15 * gk_156[k]
                   + pb_y[k] * hk_240[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_x, pa_z, pb_z, gk_139, gk_239, gk_241, \
                         gl_111, gl_152, gl_153, hk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_0 * gk_239[k]
                   + pa_x[k] * gl_152[k];

        t_505[k] = pa_z[k] * gl_111[k];

        t_506[k] = f_13 * gk_139[k]
                   + pb_z[k] * hk_241[k];

        t_507[k] = f_16 * gk_241[k]
                   + pa_x[k] * gl_153[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_x, pa_z, pb_y, pb_z, gk_141, gk_158, \
                         gk_242, gl_112, gl_154, hk_242, hk_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_15 * gk_158[k]
                   + pb_y[k] * hk_242[k];

        t_509[k] = f_16 * gk_242[k]
                   + pa_x[k] * gl_154[k];

        t_510[k] = pa_z[k] * gl_112[k];

        t_511[k] = f_13 * gk_141[k]
                   + pb_z[k] * hk_243[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_x, pb_y, gk_160, gk_244, gk_245, \
                         gk_246, gl_155, gl_156, gl_157, hk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_15 * gk_244[k]
                   + pa_x[k] * gl_155[k];

        t_513[k] = f_15 * gk_245[k]
                   + pa_x[k] * gl_156[k];

        t_514[k] = f_15 * gk_160[k]
                   + pb_y[k] * hk_244[k];

        t_515[k] = f_15 * gk_246[k]
                   + pa_x[k] * gl_157[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pa_x, pa_z, pb_z, gk_143, gk_247, gk_248, \
                         gl_113, gl_158, gl_159, hk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = pa_z[k] * gl_113[k];

        t_517[k] = f_13 * gk_143[k]
                   + pb_z[k] * hk_245[k];

        t_518[k] = f_14 * gk_247[k]
                   + pa_x[k] * gl_158[k];

        t_519[k] = f_14 * gk_248[k]
                   + pa_x[k] * gl_159[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pa_x, pa_z, pb_y, gk_162, gk_249, gk_250, \
                         gl_114, gl_160, gl_161, hk_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_14 * gk_249[k]
                   + pa_x[k] * gl_160[k];

        t_521[k] = f_15 * gk_162[k]
                   + pb_y[k] * hk_246[k];

        t_522[k] = f_14 * gk_250[k]
                   + pa_x[k] * gl_161[k];

        t_523[k] = pa_z[k] * gl_114[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, pb_x, gk_252, gk_253, gk_254, \
                         gk_255, gk_256, hk_247, hk_248, hk_249, hk_250, \
                         hk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_13 * gk_252[k]
                   + pb_x[k] * hk_247[k];

        t_525[k] = f_13 * gk_253[k]
                   + pb_x[k] * hk_248[k];

        t_526[k] = f_13 * gk_254[k]
                   + pb_x[k] * hk_249[k];

        t_527[k] = f_13 * gk_255[k]
                   + pb_x[k] * hk_250[k];

        t_528[k] = f_13 * gk_256[k]
                   + pb_x[k] * hk_251[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, t_534, pa_x, pb_x, gk_257, gk_258, \
                         gl_162, gl_163, gl_164, gl_165, hk_252, \
                         hk_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_13 * gk_257[k]
                   + pb_x[k] * hk_252[k];

        t_530[k] = f_13 * gk_258[k]
                   + pb_x[k] * hk_253[k];

        t_531[k] = pa_x[k] * gl_162[k];

        t_532[k] = pa_x[k] * gl_163[k];

        t_533[k] = pa_x[k] * gl_164[k];

        t_534[k] = pa_x[k] * gl_165[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, t_540, pa_x, gk_259, gl_166, \
                         gl_167, gl_168, gl_169, gl_170, gl_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = pa_x[k] * gl_166[k];

        t_536[k] = pa_x[k] * gl_167[k];

        t_537[k] = pa_x[k] * gl_168[k];

        t_538[k] = pa_x[k] * gl_169[k];

        t_539[k] = pa_x[k] * gl_170[k];

        t_540[k] = f_18 * gk_259[k]
                   + pa_x[k] * gl_171[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pa_x, pb_y, pb_z, gk_153, gk_170, gk_171, \
                         gk_261, gl_172, hk_254, hk_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_14 * gk_170[k]
                   + pb_y[k] * hk_254[k];

        t_542[k] = f_14 * gk_153[k]
                   + pb_z[k] * hk_254[k];

        t_543[k] = f_17 * gk_261[k]
                   + pa_x[k] * gl_172[k];

        t_544[k] = f_14 * gk_171[k]
                   + pb_y[k] * hk_255[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pa_x, pb_y, pb_z, gk_155, gk_173, gk_262, \
                         gk_263, gl_173, gl_174, hk_256, hk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_17 * gk_262[k]
                   + pa_x[k] * gl_173[k];

        t_546[k] = f_0 * gk_263[k]
                   + pa_x[k] * gl_174[k];

        t_547[k] = f_14 * gk_155[k]
                   + pb_z[k] * hk_256[k];

        t_548[k] = f_14 * gk_173[k]
                   + pb_y[k] * hk_257[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pa_x, pb_z, gk_157, gk_264, gk_265, \
                         gk_266, gl_175, gl_176, gl_177, hk_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_0 * gk_264[k]
                   + pa_x[k] * gl_175[k];

        t_550[k] = f_16 * gk_265[k]
                   + pa_x[k] * gl_176[k];

        t_551[k] = f_14 * gk_157[k]
                   + pb_z[k] * hk_258[k];

        t_552[k] = f_16 * gk_266[k]
                   + pa_x[k] * gl_177[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pa_x, pb_y, pb_z, gk_159, gk_175, gk_267, \
                         gk_268, gl_178, gl_179, hk_259, hk_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_14 * gk_175[k]
                   + pb_y[k] * hk_259[k];

        t_554[k] = f_16 * gk_267[k]
                   + pa_x[k] * gl_178[k];

        t_555[k] = f_15 * gk_268[k]
                   + pa_x[k] * gl_179[k];

        t_556[k] = f_14 * gk_159[k]
                   + pb_z[k] * hk_260[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, pa_x, pb_y, gk_177, gk_269, gk_270, \
                         gk_271, gl_180, gl_181, gl_182, hk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_15 * gk_269[k]
                   + pa_x[k] * gl_180[k];

        t_558[k] = f_15 * gk_270[k]
                   + pa_x[k] * gl_181[k];

        t_559[k] = f_14 * gk_177[k]
                   + pb_y[k] * hk_261[k];

        t_560[k] = f_15 * gk_271[k]
                   + pa_x[k] * gl_182[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, pa_x, pb_z, gk_161, gk_272, gk_273, \
                         gk_274, gl_183, gl_184, gl_185, hk_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_14 * gk_272[k]
                   + pa_x[k] * gl_183[k];

        t_562[k] = f_14 * gk_161[k]
                   + pb_z[k] * hk_262[k];

        t_563[k] = f_14 * gk_273[k]
                   + pa_x[k] * gl_184[k];

        t_564[k] = f_14 * gk_274[k]
                   + pa_x[k] * gl_185[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pa_x, pb_x, pb_y, gk_179, gk_275, gk_276, \
                         gk_277, gl_186, gl_187, hk_263, hk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_14 * gk_275[k]
                   + pa_x[k] * gl_186[k];

        t_566[k] = f_14 * gk_179[k]
                   + pb_y[k] * hk_263[k];

        t_567[k] = f_14 * gk_276[k]
                   + pa_x[k] * gl_187[k];

        t_568[k] = f_13 * gk_277[k]
                   + pb_x[k] * hk_264[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, pb_x, gk_278, gk_279, gk_280, \
                         gk_281, gk_282, hk_265, hk_266, hk_267, hk_268, \
                         hk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_13 * gk_278[k]
                   + pb_x[k] * hk_265[k];

        t_570[k] = f_13 * gk_279[k]
                   + pb_x[k] * hk_266[k];

        t_571[k] = f_13 * gk_280[k]
                   + pb_x[k] * hk_267[k];

        t_572[k] = f_13 * gk_281[k]
                   + pb_x[k] * hk_268[k];

        t_573[k] = f_13 * gk_282[k]
                   + pb_x[k] * hk_269[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, t_579, pa_x, pb_x, gk_283, gk_284, \
                         gl_188, gl_189, gl_190, gl_191, hk_270, \
                         hk_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_13 * gk_283[k]
                   + pb_x[k] * hk_270[k];

        t_575[k] = f_13 * gk_284[k]
                   + pb_x[k] * hk_271[k];

        t_576[k] = pa_x[k] * gl_188[k];

        t_577[k] = pa_x[k] * gl_189[k];

        t_578[k] = pa_x[k] * gl_190[k];

        t_579[k] = pa_x[k] * gl_191[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, t_585, pa_x, pa_y, gl_116, gl_192, \
                         gl_193, gl_194, gl_195, gl_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = pa_x[k] * gl_192[k];

        t_581[k] = pa_x[k] * gl_193[k];

        t_582[k] = pa_x[k] * gl_194[k];

        t_583[k] = pa_x[k] * gl_195[k];

        t_584[k] = pa_x[k] * gl_196[k];

        t_585[k] = pa_y[k] * gl_116[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, pa_x, pa_y, pb_y, gk_187, gk_188, \
                         gk_287, gl_117, gl_118, gl_197, hk_272, \
                         hk_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_13 * gk_187[k]
                   + pb_y[k] * hk_272[k];

        t_587[k] = pa_y[k] * gl_117[k];

        t_588[k] = f_17 * gk_287[k]
                   + pa_x[k] * gl_197[k];

        t_589[k] = f_13 * gk_188[k]
                   + pb_y[k] * hk_273[k];

        t_590[k] = pa_y[k] * gl_118[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pa_x, pa_y, pb_y, pb_z, gk_172, gk_190, \
                         gk_289, gl_119, gl_198, hk_274, hk_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_0 * gk_289[k]
                   + pa_x[k] * gl_198[k];

        t_592[k] = f_15 * gk_172[k]
                   + pb_z[k] * hk_274[k];

        t_593[k] = f_13 * gk_190[k]
                   + pb_y[k] * hk_275[k];

        t_594[k] = pa_y[k] * gl_119[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pa_x, pb_y, pb_z, gk_174, gk_192, gk_291, \
                         gk_292, gl_199, gl_200, hk_276, hk_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_16 * gk_291[k]
                   + pa_x[k] * gl_199[k];

        t_596[k] = f_15 * gk_174[k]
                   + pb_z[k] * hk_276[k];

        t_597[k] = f_16 * gk_292[k]
                   + pa_x[k] * gl_200[k];

        t_598[k] = f_13 * gk_192[k]
                   + pb_y[k] * hk_277[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, pa_x, pa_y, pb_z, gk_176, gk_294, gk_295, \
                         gl_120, gl_201, gl_202, hk_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pa_y[k] * gl_120[k];

        t_600[k] = f_15 * gk_294[k]
                   + pa_x[k] * gl_201[k];

        t_601[k] = f_15 * gk_176[k]
                   + pb_z[k] * hk_278[k];

        t_602[k] = f_15 * gk_295[k]
                   + pa_x[k] * gl_202[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pa_x, pa_y, pb_y, gk_194, gk_296, gk_298, \
                         gl_121, gl_203, gl_204, hk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_15 * gk_296[k]
                   + pa_x[k] * gl_203[k];

        t_604[k] = f_13 * gk_194[k]
                   + pb_y[k] * hk_279[k];

        t_605[k] = pa_y[k] * gl_121[k];

        t_606[k] = f_14 * gk_298[k]
                   + pa_x[k] * gl_204[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pa_x, pb_z, gk_178, gk_299, gk_300, \
                         gk_301, gl_205, gl_206, gl_207, hk_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_15 * gk_178[k]
                   + pb_z[k] * hk_280[k];

        t_608[k] = f_14 * gk_299[k]
                   + pa_x[k] * gl_205[k];

        t_609[k] = f_14 * gk_300[k]
                   + pa_x[k] * gl_206[k];

        t_610[k] = f_14 * gk_301[k]
                   + pa_x[k] * gl_207[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pa_y, pb_x, pb_y, gk_196, gk_302, gk_303, \
                         gl_122, hk_281, hk_282, hk_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_13 * gk_196[k]
                   + pb_y[k] * hk_281[k];

        t_612[k] = pa_y[k] * gl_122[k];

        t_613[k] = f_13 * gk_302[k]
                   + pb_x[k] * hk_282[k];

        t_614[k] = f_13 * gk_303[k]
                   + pb_x[k] * hk_283[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, pb_x, gk_304, gk_305, gk_306, \
                         gk_307, gk_308, hk_284, hk_285, hk_286, hk_287, \
                         hk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_13 * gk_304[k]
                   + pb_x[k] * hk_284[k];

        t_616[k] = f_13 * gk_305[k]
                   + pb_x[k] * hk_285[k];

        t_617[k] = f_13 * gk_306[k]
                   + pb_x[k] * hk_286[k];

        t_618[k] = f_13 * gk_307[k]
                   + pb_x[k] * hk_287[k];

        t_619[k] = f_13 * gk_308[k]
                   + pb_x[k] * hk_288[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, t_625, t_626, pa_x, pa_y, gl_123, \
                         gl_208, gl_209, gl_210, gl_211, gl_212, \
                         gl_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = pa_y[k] * gl_123[k];

        t_621[k] = pa_x[k] * gl_208[k];

        t_622[k] = pa_x[k] * gl_209[k];

        t_623[k] = pa_x[k] * gl_210[k];

        t_624[k] = pa_x[k] * gl_211[k];

        t_625[k] = pa_x[k] * gl_212[k];

        t_626[k] = pa_x[k] * gl_213[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, t_632, pa_x, pb_y, pb_z, gk_187, \
                         gk_310, gl_214, gl_215, gl_216, gl_217, \
                         hk_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = pa_x[k] * gl_214[k];

        t_628[k] = pa_x[k] * gl_215[k];

        t_629[k] = pa_x[k] * gl_216[k];

        t_630[k] = f_18 * gk_310[k]
                   + pa_x[k] * gl_217[k];

        t_631[k] = pb_y[k] * hk_289[k];

        t_632[k] = f_16 * gk_187[k]
                   + pb_z[k] * hk_289[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pa_x, pb_y, gk_313, gk_314, gk_315, \
                         gl_219, gl_220, gl_221, hk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_17 * gk_313[k]
                   + pa_x[k] * gl_219[k];

        t_634[k] = pb_y[k] * hk_290[k];

        t_635[k] = f_17 * gk_314[k]
                   + pa_x[k] * gl_220[k];

        t_636[k] = f_0 * gk_315[k]
                   + pa_x[k] * gl_221[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pa_x, pb_y, pb_z, gk_189, gk_317, gk_318, \
                         gl_222, gl_223, hk_291, hk_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = f_16 * gk_189[k]
                   + pb_z[k] * hk_291[k];

        t_638[k] = pb_y[k] * hk_292[k];

        t_639[k] = f_0 * gk_317[k]
                   + pa_x[k] * gl_222[k];

        t_640[k] = f_16 * gk_318[k]
                   + pa_x[k] * gl_223[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, pa_x, pb_y, pb_z, gk_191, gk_319, gk_321, \
                         gl_224, gl_225, hk_293, hk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_16 * gk_191[k]
                   + pb_z[k] * hk_293[k];

        t_642[k] = f_16 * gk_319[k]
                   + pa_x[k] * gl_224[k];

        t_643[k] = pb_y[k] * hk_294[k];

        t_644[k] = f_16 * gk_321[k]
                   + pa_x[k] * gl_225[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, pa_x, pb_z, gk_193, gk_322, gk_323, \
                         gk_324, gl_226, gl_227, gl_228, hk_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_15 * gk_322[k]
                   + pa_x[k] * gl_226[k];

        t_646[k] = f_16 * gk_193[k]
                   + pb_z[k] * hk_295[k];

        t_647[k] = f_15 * gk_323[k]
                   + pa_x[k] * gl_227[k];

        t_648[k] = f_15 * gk_324[k]
                   + pa_x[k] * gl_228[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pa_x, pb_y, pb_z, gk_195, gk_326, gk_327, \
                         gl_229, gl_230, hk_296, hk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pb_y[k] * hk_296[k];

        t_650[k] = f_15 * gk_326[k]
                   + pa_x[k] * gl_229[k];

        t_651[k] = f_14 * gk_327[k]
                   + pa_x[k] * gl_230[k];

        t_652[k] = f_16 * gk_195[k]
                   + pb_z[k] * hk_297[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, t_657, pa_x, pb_y, gk_328, gk_329, \
                         gk_330, gk_331, gl_231, gl_232, gl_233, gl_234, \
                         hk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_14 * gk_328[k]
                   + pa_x[k] * gl_231[k];

        t_654[k] = f_14 * gk_329[k]
                   + pa_x[k] * gl_232[k];

        t_655[k] = f_14 * gk_330[k]
                   + pa_x[k] * gl_233[k];

        t_656[k] = pb_y[k] * hk_298[k];

        t_657[k] = f_14 * gk_331[k]
                   + pa_x[k] * gl_234[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, t_662, pb_x, gk_332, gk_333, gk_334, \
                         gk_335, gk_336, hk_300, hk_301, hk_302, hk_303, \
                         hk_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_13 * gk_332[k]
                   + pb_x[k] * hk_300[k];

        t_659[k] = f_13 * gk_333[k]
                   + pb_x[k] * hk_301[k];

        t_660[k] = f_13 * gk_334[k]
                   + pb_x[k] * hk_302[k];

        t_661[k] = f_13 * gk_335[k]
                   + pb_x[k] * hk_303[k];

        t_662[k] = f_13 * gk_336[k]
                   + pb_x[k] * hk_304[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, t_667, pa_x, pb_x, pb_y, gk_337, gk_339, \
                         gl_235, gl_236, hk_299, hk_305, hk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_13 * gk_337[k]
                   + pb_x[k] * hk_305[k];

        t_664[k] = pb_y[k] * hk_299[k];

        t_665[k] = f_13 * gk_339[k]
                   + pb_x[k] * hk_306[k];

        t_666[k] = pa_x[k] * gl_235[k];

        t_667[k] = pa_x[k] * gl_236[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, t_672, t_673, t_674, pa_x, pb_y, gl_237, \
                         gl_238, gl_239, gl_240, gl_241, gl_242, \
                         hk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = pa_x[k] * gl_237[k];

        t_669[k] = pa_x[k] * gl_238[k];

        t_670[k] = pa_x[k] * gl_239[k];

        t_671[k] = pa_x[k] * gl_240[k];

        t_672[k] = pa_x[k] * gl_241[k];

        t_673[k] = pb_y[k] * hk_306[k];

        t_674[k] = pa_x[k] * gl_242[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, pb_x, pb_y, pb_z, gk_205, hi0_90, \
                         hi0_91, hi1_90, hi1_91, hk_307, hk_308, \
                         hk_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = f_1 * hi0_90[k]
                   - f_2 * hi1_90[k]
                   + pb_x[k] * hk_307[k];

        t_676[k] = f_0 * gk_205[k]
                   + pb_y[k] * hk_307[k];

        t_677[k] = pb_z[k] * hk_307[k];

        t_678[k] = f_11 * hi0_91[k]
                   - f_12 * hi1_91[k]
                   + pb_x[k] * hk_309[k];

        t_679[k] = pb_z[k] * hk_308[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pb_x, pb_y, pb_z, gk_208, hi0_92, hi0_93, \
                         hi1_92, hi1_93, hk_309, hk_310, hk_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * hi0_92[k]
                   - f_12 * hi1_92[k]
                   + pb_x[k] * hk_310[k];

        t_681[k] = f_9 * hi0_93[k]
                   - f_10 * hi1_93[k]
                   + pb_x[k] * hk_311[k];

        t_682[k] = pb_z[k] * hk_309[k];

        t_683[k] = f_0 * gk_208[k]
                   + pb_y[k] * hk_310[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pb_x, pb_z, hi0_94, hi0_95, hi0_96, \
                         hi1_94, hi1_95, hi1_96, hk_311, hk_312, hk_313, \
                         hk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_9 * hi0_94[k]
                   - f_10 * hi1_94[k]
                   + pb_x[k] * hk_312[k];

        t_685[k] = f_7 * hi0_95[k]
                   - f_8 * hi1_95[k]
                   + pb_x[k] * hk_313[k];

        t_686[k] = pb_z[k] * hk_311[k];

        t_687[k] = f_7 * hi0_96[k]
                   - f_8 * hi1_96[k]
                   + pb_x[k] * hk_314[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, t_691, pb_x, pb_y, pb_z, gk_211, hi0_97, hi0_98, \
                         hi1_97, hi1_98, hk_312, hk_313, hk_315, \
                         hk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_0 * gk_211[k]
                   + pb_y[k] * hk_312[k];

        t_689[k] = f_7 * hi0_97[k]
                   - f_8 * hi1_97[k]
                   + pb_x[k] * hk_315[k];

        t_690[k] = f_5 * hi0_98[k]
                   - f_6 * hi1_98[k]
                   + pb_x[k] * hk_316[k];

        t_691[k] = pb_z[k] * hk_313[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pb_x, pb_y, gk_215, hi0_99, hi0_100, hi1_99, \
                         hi1_100, hk_315, hk_317, hk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_5 * hi0_99[k]
                   - f_6 * hi1_99[k]
                   + pb_x[k] * hk_317[k];

        t_693[k] = f_5 * hi0_100[k]
                   - f_6 * hi1_100[k]
                   + pb_x[k] * hk_318[k];

        t_694[k] = f_0 * gk_215[k]
                   + pb_y[k] * hk_315[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pb_x, pb_z, hi0_101, hi0_102, hi0_104, \
                         hi1_101, hi1_102, hi1_104, hk_316, hk_319, hk_320, \
                         hk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_5 * hi0_101[k]
                   - f_6 * hi1_101[k]
                   + pb_x[k] * hk_319[k];

        t_696[k] = f_3 * hi0_102[k]
                   - f_4 * hi1_102[k]
                   + pb_x[k] * hk_320[k];

        t_697[k] = pb_z[k] * hk_316[k];

        t_698[k] = f_3 * hi0_104[k]
                   - f_4 * hi1_104[k]
                   + pb_x[k] * hk_321[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pb_x, pb_y, gk_220, hi0_105, hi0_106, hi1_105, \
                         hi1_106, hk_319, hk_322, hk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_3 * hi0_105[k]
                   - f_4 * hi1_105[k]
                   + pb_x[k] * hk_322[k];

        t_700[k] = f_3 * hi0_106[k]
                   - f_4 * hi1_106[k]
                   + pb_x[k] * hk_323[k];

        t_701[k] = f_0 * gk_220[k]
                   + pb_y[k] * hk_319[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, t_707, pb_x, hi0_107, hi1_107, \
                         hk_324, hk_325, hk_326, hk_327, hk_328, \
                         hk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_3 * hi0_107[k]
                   - f_4 * hi1_107[k]
                   + pb_x[k] * hk_324[k];

        t_703[k] = pb_x[k] * hk_325[k];

        t_704[k] = pb_x[k] * hk_326[k];

        t_705[k] = pb_x[k] * hk_327[k];

        t_706[k] = pb_x[k] * hk_328[k];

        t_707[k] = pb_x[k] * hk_329[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, t_712, pb_x, pb_y, pb_z, gk_226, hi0_102, \
                         hi1_102, hk_325, hk_330, hk_331, hk_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = pb_x[k] * hk_330[k];

        t_709[k] = pb_x[k] * hk_331[k];

        t_710[k] = pb_x[k] * hk_332[k];

        t_711[k] = f_0 * gk_226[k]
                   + f_1 * hi0_102[k]
                   - f_2 * hi1_102[k]
                   + pb_y[k] * hk_325[k];

        t_712[k] = pb_z[k] * hk_325[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, pb_z, hi0_102, hi0_103, hi0_104, hi1_102, \
                         hi1_103, hi1_104, hk_326, hk_327, hk_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_3 * hi0_102[k]
                   - f_4 * hi1_102[k]
                   + pb_z[k] * hk_326[k];

        t_714[k] = f_5 * hi0_103[k]
                   - f_6 * hi1_103[k]
                   + pb_z[k] * hk_327[k];

        t_715[k] = f_7 * hi0_104[k]
                   - f_8 * hi1_104[k]
                   + pb_z[k] * hk_328[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pb_y, pb_z, gk_233, hi0_105, hi0_106, \
                         hi0_107, hi1_105, hi1_106, hi1_107, hk_329, hk_330, \
                         hk_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_9 * hi0_105[k]
                   - f_10 * hi1_105[k]
                   + pb_z[k] * hk_329[k];

        t_717[k] = f_11 * hi0_106[k]
                   - f_12 * hi1_106[k]
                   + pb_z[k] * hk_330[k];

        t_718[k] = f_0 * gk_233[k]
                   + pb_y[k] * hk_332[k];

        t_719[k] = f_1 * hi0_107[k]
                   - f_2 * hi1_107[k]
                   + pb_z[k] * hk_332[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pa_z, pb_y, pb_z, gk_205, gk_235, \
                         gl_125, gl_126, gl_127, hk_333, hk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_z[k] * gl_125[k];

        t_721[k] = pa_z[k] * gl_126[k];

        t_722[k] = f_13 * gk_205[k]
                   + pb_z[k] * hk_333[k];

        t_723[k] = pa_z[k] * gl_127[k];

        t_724[k] = f_16 * gk_235[k]
                   + pb_y[k] * hk_334[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pa_z, pb_y, pb_z, gk_206, gk_207, gk_237, \
                         gl_128, gl_129, hk_335, hk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_14 * gk_206[k]
                   + pa_z[k] * gl_128[k];

        t_726[k] = pa_z[k] * gl_129[k];

        t_727[k] = f_13 * gk_207[k]
                   + pb_z[k] * hk_335[k];

        t_728[k] = f_16 * gk_237[k]
                   + pb_y[k] * hk_336[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, pa_z, pb_z, gk_208, gk_209, gk_210, \
                         gl_130, gl_131, gl_132, hk_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_15 * gk_208[k]
                   + pa_z[k] * gl_130[k];

        t_730[k] = pa_z[k] * gl_131[k];

        t_731[k] = f_13 * gk_209[k]
                   + pb_z[k] * hk_337[k];

        t_732[k] = f_14 * gk_210[k]
                   + pa_z[k] * gl_132[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, pa_z, pb_y, pb_z, gk_211, gk_212, gk_239, \
                         gl_133, gl_134, hk_338, hk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_16 * gk_239[k]
                   + pb_y[k] * hk_338[k];

        t_734[k] = f_16 * gk_211[k]
                   + pa_z[k] * gl_133[k];

        t_735[k] = pa_z[k] * gl_134[k];

        t_736[k] = f_13 * gk_212[k]
                   + pb_z[k] * hk_339[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, pa_z, pb_y, gk_213, gk_214, \
                         gk_215, gk_242, gl_135, gl_136, gl_137, gl_138, \
                         hk_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_14 * gk_213[k]
                   + pa_z[k] * gl_135[k];

        t_738[k] = f_15 * gk_214[k]
                   + pa_z[k] * gl_136[k];

        t_739[k] = f_16 * gk_242[k]
                   + pb_y[k] * hk_340[k];

        t_740[k] = f_0 * gk_215[k]
                   + pa_z[k] * gl_137[k];

        t_741[k] = pa_z[k] * gl_138[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pa_z, pb_z, gk_216, gk_217, gk_218, \
                         gk_219, gl_139, gl_140, gl_141, hk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_13 * gk_216[k]
                   + pb_z[k] * hk_341[k];

        t_743[k] = f_14 * gk_217[k]
                   + pa_z[k] * gl_139[k];

        t_744[k] = f_15 * gk_218[k]
                   + pa_z[k] * gl_140[k];

        t_745[k] = f_16 * gk_219[k]
                   + pa_z[k] * gl_141[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, t_750, pa_z, pb_x, pb_y, gk_220, gk_246, \
                         gl_142, hk_342, hk_343, hk_344, hk_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_16 * gk_246[k]
                   + pb_y[k] * hk_342[k];

        t_747[k] = f_17 * gk_220[k]
                   + pa_z[k] * gl_142[k];

        t_748[k] = pb_x[k] * hk_343[k];

        t_749[k] = pb_x[k] * hk_344[k];

        t_750[k] = pb_x[k] * hk_345[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, t_754, t_755, t_756, pa_z, pb_x, gl_143, hk_346, \
                         hk_347, hk_348, hk_349, hk_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = pb_x[k] * hk_346[k];

        t_752[k] = pb_x[k] * hk_347[k];

        t_753[k] = pb_x[k] * hk_348[k];

        t_754[k] = pb_x[k] * hk_349[k];

        t_755[k] = pb_x[k] * hk_350[k];

        t_756[k] = pa_z[k] * gl_143[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pa_z, pb_z, gk_226, gk_227, gk_228, \
                         gk_229, gl_144, gl_145, gl_146, hk_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_13 * gk_226[k]
                   + pb_z[k] * hk_343[k];

        t_758[k] = f_14 * gk_227[k]
                   + pa_z[k] * gl_144[k];

        t_759[k] = f_15 * gk_228[k]
                   + pa_z[k] * gl_145[k];

        t_760[k] = f_16 * gk_229[k]
                   + pa_z[k] * gl_146[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pa_z, pb_y, gk_230, gk_231, gk_233, \
                         gk_258, gl_147, gl_148, gl_150, hk_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_0 * gk_230[k]
                   + pa_z[k] * gl_147[k];

        t_762[k] = f_17 * gk_231[k]
                   + pa_z[k] * gl_148[k];

        t_763[k] = f_16 * gk_258[k]
                   + pb_y[k] * hk_350[k];

        t_764[k] = f_18 * gk_233[k]
                   + pa_z[k] * gl_150[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, pb_x, pb_y, pb_z, gk_234, gk_259, \
                         hi0_108, hi0_109, hi1_108, hi1_109, hk_351, \
                         hk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_1 * hi0_108[k]
                   - f_2 * hi1_108[k]
                   + pb_x[k] * hk_351[k];

        t_766[k] = f_15 * gk_259[k]
                   + pb_y[k] * hk_351[k];

        t_767[k] = f_14 * gk_234[k]
                   + pb_z[k] * hk_351[k];

        t_768[k] = f_11 * hi0_109[k]
                   - f_12 * hi1_109[k]
                   + pb_x[k] * hk_353[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, pb_x, pb_y, gk_260, hi0_110, hi0_111, hi1_110, \
                         hi1_111, hk_352, hk_354, hk_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_15 * gk_260[k]
                   + pb_y[k] * hk_352[k];

        t_770[k] = f_11 * hi0_110[k]
                   - f_12 * hi1_110[k]
                   + pb_x[k] * hk_354[k];

        t_771[k] = f_9 * hi0_111[k]
                   - f_10 * hi1_111[k]
                   + pb_x[k] * hk_355[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pb_x, pb_y, pb_z, gk_236, gk_262, hi0_112, \
                         hi1_112, hk_353, hk_354, hk_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_14 * gk_236[k]
                   + pb_z[k] * hk_353[k];

        t_773[k] = f_15 * gk_262[k]
                   + pb_y[k] * hk_354[k];

        t_774[k] = f_9 * hi0_112[k]
                   - f_10 * hi1_112[k]
                   + pb_x[k] * hk_356[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, pb_x, pb_z, gk_238, hi0_113, hi0_114, hi1_113, \
                         hi1_114, hk_355, hk_357, hk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_7 * hi0_113[k]
                   - f_8 * hi1_113[k]
                   + pb_x[k] * hk_357[k];

        t_776[k] = f_14 * gk_238[k]
                   + pb_z[k] * hk_355[k];

        t_777[k] = f_7 * hi0_114[k]
                   - f_8 * hi1_114[k]
                   + pb_x[k] * hk_358[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, pb_x, pb_y, gk_264, hi0_115, hi0_116, hi1_115, \
                         hi1_116, hk_356, hk_359, hk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_15 * gk_264[k]
                   + pb_y[k] * hk_356[k];

        t_779[k] = f_7 * hi0_115[k]
                   - f_8 * hi1_115[k]
                   + pb_x[k] * hk_359[k];

        t_780[k] = f_5 * hi0_116[k]
                   - f_6 * hi1_116[k]
                   + pb_x[k] * hk_360[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, pb_x, pb_z, gk_240, hi0_117, hi0_118, hi1_117, \
                         hi1_118, hk_357, hk_361, hk_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_14 * gk_240[k]
                   + pb_z[k] * hk_357[k];

        t_782[k] = f_5 * hi0_117[k]
                   - f_6 * hi1_117[k]
                   + pb_x[k] * hk_361[k];

        t_783[k] = f_5 * hi0_118[k]
                   - f_6 * hi1_118[k]
                   + pb_x[k] * hk_362[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pb_x, pb_y, gk_267, hi0_119, hi0_120, hi1_119, \
                         hi1_120, hk_359, hk_363, hk_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_15 * gk_267[k]
                   + pb_y[k] * hk_359[k];

        t_785[k] = f_5 * hi0_119[k]
                   - f_6 * hi1_119[k]
                   + pb_x[k] * hk_363[k];

        t_786[k] = f_3 * hi0_120[k]
                   - f_4 * hi1_120[k]
                   + pb_x[k] * hk_364[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pb_x, pb_z, gk_243, hi0_121, hi0_122, hi1_121, \
                         hi1_122, hk_360, hk_365, hk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_14 * gk_243[k]
                   + pb_z[k] * hk_360[k];

        t_788[k] = f_3 * hi0_121[k]
                   - f_4 * hi1_121[k]
                   + pb_x[k] * hk_365[k];

        t_789[k] = f_3 * hi0_122[k]
                   - f_4 * hi1_122[k]
                   + pb_x[k] * hk_366[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pb_x, pb_y, gk_271, hi0_123, hi0_125, \
                         hi1_123, hi1_125, hk_363, hk_367, hk_368, \
                         hk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_3 * hi0_123[k]
                   - f_4 * hi1_123[k]
                   + pb_x[k] * hk_367[k];

        t_791[k] = f_15 * gk_271[k]
                   + pb_y[k] * hk_363[k];

        t_792[k] = f_3 * hi0_125[k]
                   - f_4 * hi1_125[k]
                   + pb_x[k] * hk_368[k];

        t_793[k] = pb_x[k] * hk_369[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, t_799, t_800, pb_x, hk_370, \
                         hk_371, hk_372, hk_373, hk_374, hk_375, \
                         hk_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = pb_x[k] * hk_370[k];

        t_795[k] = pb_x[k] * hk_371[k];

        t_796[k] = pb_x[k] * hk_372[k];

        t_797[k] = pb_x[k] * hk_373[k];

        t_798[k] = pb_x[k] * hk_374[k];

        t_799[k] = pb_x[k] * hk_375[k];

        t_800[k] = pb_x[k] * hk_376[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pa_z, pb_y, pb_z, fl0_5, fl1_5, gk_251, gk_279, \
                         gl_162, hi0_121, hi1_121, hk_369, hk_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_19 * fl0_5[k]
                   - f_20 * fl1_5[k]
                   + pa_z[k] * gl_162[k];

        t_802[k] = f_14 * gk_251[k]
                   + pb_z[k] * hk_369[k];

        t_803[k] = f_15 * gk_279[k]
                   + f_11 * hi0_121[k]
                   - f_12 * hi1_121[k]
                   + pb_y[k] * hk_371[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pb_y, gk_280, gk_281, gk_282, hi0_122, hi0_123, \
                         hi0_124, hi1_122, hi1_123, hi1_124, hk_372, hk_373, \
                         hk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_15 * gk_280[k]
                   + f_9 * hi0_122[k]
                   - f_10 * hi1_122[k]
                   + pb_y[k] * hk_372[k];

        t_805[k] = f_15 * gk_281[k]
                   + f_7 * hi0_123[k]
                   - f_8 * hi1_123[k]
                   + pb_y[k] * hk_373[k];

        t_806[k] = f_15 * gk_282[k]
                   + f_5 * hi0_124[k]
                   - f_6 * hi1_124[k]
                   + pb_y[k] * hk_374[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, pa_y, pb_y, fl0_7, fl1_7, gk_283, gk_284, \
                         gl_196, hi0_125, hi1_125, hk_375, hk_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_15 * gk_283[k]
                   + f_3 * hi0_125[k]
                   - f_4 * hi1_125[k]
                   + pb_y[k] * hk_375[k];

        t_808[k] = f_15 * gk_284[k]
                   + pb_y[k] * hk_376[k];

        t_809[k] = f_21 * fl0_7[k]
                   - f_22 * fl1_7[k]
                   + pa_y[k] * gl_196[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pb_x, pb_y, pb_z, gk_259, gk_285, \
                         hi0_126, hi0_127, hi1_126, hi1_127, hk_377, \
                         hk_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_1 * hi0_126[k]
                   - f_2 * hi1_126[k]
                   + pb_x[k] * hk_377[k];

        t_811[k] = f_14 * gk_285[k]
                   + pb_y[k] * hk_377[k];

        t_812[k] = f_15 * gk_259[k]
                   + pb_z[k] * hk_377[k];

        t_813[k] = f_11 * hi0_127[k]
                   - f_12 * hi1_127[k]
                   + pb_x[k] * hk_379[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pb_x, pb_y, gk_286, hi0_128, hi0_129, hi1_128, \
                         hi1_129, hk_378, hk_380, hk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_14 * gk_286[k]
                   + pb_y[k] * hk_378[k];

        t_815[k] = f_11 * hi0_128[k]
                   - f_12 * hi1_128[k]
                   + pb_x[k] * hk_380[k];

        t_816[k] = f_9 * hi0_129[k]
                   - f_10 * hi1_129[k]
                   + pb_x[k] * hk_381[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pb_x, pb_y, pb_z, gk_261, gk_288, hi0_130, \
                         hi1_130, hk_379, hk_380, hk_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_15 * gk_261[k]
                   + pb_z[k] * hk_379[k];

        t_818[k] = f_14 * gk_288[k]
                   + pb_y[k] * hk_380[k];

        t_819[k] = f_9 * hi0_130[k]
                   - f_10 * hi1_130[k]
                   + pb_x[k] * hk_382[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pb_x, pb_z, gk_263, hi0_131, hi0_132, hi1_131, \
                         hi1_132, hk_381, hk_383, hk_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_7 * hi0_131[k]
                   - f_8 * hi1_131[k]
                   + pb_x[k] * hk_383[k];

        t_821[k] = f_15 * gk_263[k]
                   + pb_z[k] * hk_381[k];

        t_822[k] = f_7 * hi0_132[k]
                   - f_8 * hi1_132[k]
                   + pb_x[k] * hk_384[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pb_x, pb_y, gk_290, hi0_133, hi0_134, hi1_133, \
                         hi1_134, hk_382, hk_385, hk_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_14 * gk_290[k]
                   + pb_y[k] * hk_382[k];

        t_824[k] = f_7 * hi0_133[k]
                   - f_8 * hi1_133[k]
                   + pb_x[k] * hk_385[k];

        t_825[k] = f_5 * hi0_134[k]
                   - f_6 * hi1_134[k]
                   + pb_x[k] * hk_386[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, pb_x, pb_z, gk_265, hi0_135, hi0_136, hi1_135, \
                         hi1_136, hk_383, hk_387, hk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_15 * gk_265[k]
                   + pb_z[k] * hk_383[k];

        t_827[k] = f_5 * hi0_135[k]
                   - f_6 * hi1_135[k]
                   + pb_x[k] * hk_387[k];

        t_828[k] = f_5 * hi0_136[k]
                   - f_6 * hi1_136[k]
                   + pb_x[k] * hk_388[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, pb_x, pb_y, gk_293, hi0_137, hi0_138, hi1_137, \
                         hi1_138, hk_385, hk_389, hk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_14 * gk_293[k]
                   + pb_y[k] * hk_385[k];

        t_830[k] = f_5 * hi0_137[k]
                   - f_6 * hi1_137[k]
                   + pb_x[k] * hk_389[k];

        t_831[k] = f_3 * hi0_138[k]
                   - f_4 * hi1_138[k]
                   + pb_x[k] * hk_390[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, pb_x, pb_z, gk_268, hi0_139, hi0_140, hi1_139, \
                         hi1_140, hk_386, hk_391, hk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_15 * gk_268[k]
                   + pb_z[k] * hk_386[k];

        t_833[k] = f_3 * hi0_139[k]
                   - f_4 * hi1_139[k]
                   + pb_x[k] * hk_391[k];

        t_834[k] = f_3 * hi0_140[k]
                   - f_4 * hi1_140[k]
                   + pb_x[k] * hk_392[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, pb_x, pb_y, gk_297, hi0_141, hi0_143, \
                         hi1_141, hi1_143, hk_389, hk_393, hk_394, \
                         hk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_3 * hi0_141[k]
                   - f_4 * hi1_141[k]
                   + pb_x[k] * hk_393[k];

        t_836[k] = f_14 * gk_297[k]
                   + pb_y[k] * hk_389[k];

        t_837[k] = f_3 * hi0_143[k]
                   - f_4 * hi1_143[k]
                   + pb_x[k] * hk_394[k];

        t_838[k] = pb_x[k] * hk_395[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, t_842, t_843, t_844, t_845, pb_x, hk_396, \
                         hk_397, hk_398, hk_399, hk_400, hk_401, \
                         hk_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = pb_x[k] * hk_396[k];

        t_840[k] = pb_x[k] * hk_397[k];

        t_841[k] = pb_x[k] * hk_398[k];

        t_842[k] = pb_x[k] * hk_399[k];

        t_843[k] = pb_x[k] * hk_400[k];

        t_844[k] = pb_x[k] * hk_401[k];

        t_845[k] = pb_x[k] * hk_402[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pa_z, pb_y, pb_z, fl0_6, fl1_6, gk_277, gk_304, \
                         gl_188, hi0_139, hi1_139, hk_395, hk_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_21 * fl0_6[k]
                   - f_22 * fl1_6[k]
                   + pa_z[k] * gl_188[k];

        t_847[k] = f_15 * gk_277[k]
                   + pb_z[k] * hk_395[k];

        t_848[k] = f_14 * gk_304[k]
                   + f_11 * hi0_139[k]
                   - f_12 * hi1_139[k]
                   + pb_y[k] * hk_397[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pb_y, gk_305, gk_306, gk_307, hi0_140, hi0_141, \
                         hi0_142, hi1_140, hi1_141, hi1_142, hk_398, hk_399, \
                         hk_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_14 * gk_305[k]
                   + f_9 * hi0_140[k]
                   - f_10 * hi1_140[k]
                   + pb_y[k] * hk_398[k];

        t_850[k] = f_14 * gk_306[k]
                   + f_7 * hi0_141[k]
                   - f_8 * hi1_141[k]
                   + pb_y[k] * hk_399[k];

        t_851[k] = f_14 * gk_307[k]
                   + f_5 * hi0_142[k]
                   - f_6 * hi1_142[k]
                   + pb_y[k] * hk_400[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pa_y, pb_y, fl0_8, fl1_8, gk_308, gk_309, \
                         gl_216, gl_217, hi0_143, hi1_143, hk_401, \
                         hk_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_14 * gk_308[k]
                   + f_3 * hi0_143[k]
                   - f_4 * hi1_143[k]
                   + pb_y[k] * hk_401[k];

        t_853[k] = f_14 * gk_309[k]
                   + pb_y[k] * hk_402[k];

        t_854[k] = f_19 * fl0_8[k]
                   - f_20 * fl1_8[k]
                   + pa_y[k] * gl_216[k];

        t_855[k] = pa_y[k] * gl_217[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, t_859, t_860, pa_y, pb_y, gk_310, gk_311, \
                         gk_312, gl_218, gl_219, gl_220, hk_403, \
                         hk_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_13 * gk_310[k]
                   + pb_y[k] * hk_403[k];

        t_857[k] = pa_y[k] * gl_218[k];

        t_858[k] = f_14 * gk_311[k]
                   + pa_y[k] * gl_219[k];

        t_859[k] = f_13 * gk_312[k]
                   + pb_y[k] * hk_404[k];

        t_860[k] = pa_y[k] * gl_220[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, t_864, pa_y, pb_y, pb_z, gk_287, gk_313, gk_314, \
                         gl_221, gl_222, hk_405, hk_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_15 * gk_313[k]
                   + pa_y[k] * gl_221[k];

        t_862[k] = f_16 * gk_287[k]
                   + pb_z[k] * hk_405[k];

        t_863[k] = f_13 * gk_314[k]
                   + pb_y[k] * hk_406[k];

        t_864[k] = pa_y[k] * gl_222[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pa_y, pb_y, pb_z, gk_289, gk_315, gk_316, \
                         gk_317, gl_223, gl_224, hk_407, hk_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_16 * gk_315[k]
                   + pa_y[k] * gl_223[k];

        t_866[k] = f_16 * gk_289[k]
                   + pb_z[k] * hk_407[k];

        t_867[k] = f_14 * gk_316[k]
                   + pa_y[k] * gl_224[k];

        t_868[k] = f_13 * gk_317[k]
                   + pb_y[k] * hk_408[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, t_873, pa_y, pb_z, gk_291, gk_318, \
                         gk_319, gk_320, gl_225, gl_226, gl_227, gl_228, \
                         hk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = pa_y[k] * gl_225[k];

        t_870[k] = f_0 * gk_318[k]
                   + pa_y[k] * gl_226[k];

        t_871[k] = f_16 * gk_291[k]
                   + pb_z[k] * hk_409[k];

        t_872[k] = f_15 * gk_319[k]
                   + pa_y[k] * gl_227[k];

        t_873[k] = f_14 * gk_320[k]
                   + pa_y[k] * gl_228[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, pa_y, pb_y, pb_z, gk_294, gk_321, gk_322, \
                         gl_229, gl_230, hk_410, hk_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_13 * gk_321[k]
                   + pb_y[k] * hk_410[k];

        t_875[k] = pa_y[k] * gl_229[k];

        t_876[k] = f_17 * gk_322[k]
                   + pa_y[k] * gl_230[k];

        t_877[k] = f_16 * gk_294[k]
                   + pb_z[k] * hk_411[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, t_882, pa_y, pb_y, gk_323, gk_324, \
                         gk_325, gk_326, gl_231, gl_232, gl_233, gl_234, \
                         hk_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_16 * gk_323[k]
                   + pa_y[k] * gl_231[k];

        t_879[k] = f_15 * gk_324[k]
                   + pa_y[k] * gl_232[k];

        t_880[k] = f_14 * gk_325[k]
                   + pa_y[k] * gl_233[k];

        t_881[k] = f_13 * gk_326[k]
                   + pb_y[k] * hk_412[k];

        t_882[k] = pa_y[k] * gl_234[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, t_887, t_888, t_889, pb_x, hk_413, \
                         hk_414, hk_415, hk_416, hk_417, hk_418, \
                         hk_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = pb_x[k] * hk_413[k];

        t_884[k] = pb_x[k] * hk_414[k];

        t_885[k] = pb_x[k] * hk_415[k];

        t_886[k] = pb_x[k] * hk_416[k];

        t_887[k] = pb_x[k] * hk_417[k];

        t_888[k] = pb_x[k] * hk_418[k];

        t_889[k] = pb_x[k] * hk_419[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pa_y, pb_x, pb_z, gk_302, gk_332, gk_334, \
                         gl_235, gl_237, hk_413, hk_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = pb_x[k] * hk_420[k];

        t_891[k] = f_18 * gk_332[k]
                   + pa_y[k] * gl_235[k];

        t_892[k] = f_16 * gk_302[k]
                   + pb_z[k] * hk_413[k];

        t_893[k] = f_17 * gk_334[k]
                   + pa_y[k] * gl_237[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pa_y, gk_335, gk_336, gk_337, gk_338, \
                         gl_238, gl_239, gl_240, gl_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_0 * gk_335[k]
                   + pa_y[k] * gl_238[k];

        t_895[k] = f_16 * gk_336[k]
                   + pa_y[k] * gl_239[k];

        t_896[k] = f_15 * gk_337[k]
                   + pa_y[k] * gl_240[k];

        t_897[k] = f_14 * gk_338[k]
                   + pa_y[k] * gl_241[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, t_902, pa_y, pb_x, pb_y, pb_z, gk_310, \
                         gk_339, gl_242, hi0_144, hi1_144, hk_420, \
                         hk_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_13 * gk_339[k]
                   + pb_y[k] * hk_420[k];

        t_899[k] = pa_y[k] * gl_242[k];

        t_900[k] = f_1 * hi0_144[k]
                   - f_2 * hi1_144[k]
                   + pb_x[k] * hk_421[k];

        t_901[k] = pb_y[k] * hk_421[k];

        t_902[k] = f_0 * gk_310[k]
                   + pb_z[k] * hk_421[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, pb_x, pb_y, hi0_145, hi0_146, hi0_147, \
                         hi1_145, hi1_146, hi1_147, hk_422, hk_423, hk_424, \
                         hk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_11 * hi0_145[k]
                   - f_12 * hi1_145[k]
                   + pb_x[k] * hk_423[k];

        t_904[k] = pb_y[k] * hk_422[k];

        t_905[k] = f_11 * hi0_146[k]
                   - f_12 * hi1_146[k]
                   + pb_x[k] * hk_424[k];

        t_906[k] = f_9 * hi0_147[k]
                   - f_10 * hi1_147[k]
                   + pb_x[k] * hk_425[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, pb_x, pb_y, pb_z, gk_313, hi0_148, \
                         hi0_149, hi1_148, hi1_149, hk_423, hk_424, hk_426, \
                         hk_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_0 * gk_313[k]
                   + pb_z[k] * hk_423[k];

        t_908[k] = pb_y[k] * hk_424[k];

        t_909[k] = f_9 * hi0_148[k]
                   - f_10 * hi1_148[k]
                   + pb_x[k] * hk_426[k];

        t_910[k] = f_7 * hi0_149[k]
                   - f_8 * hi1_149[k]
                   + pb_x[k] * hk_427[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, pb_x, pb_y, pb_z, gk_315, hi0_150, \
                         hi0_151, hi1_150, hi1_151, hk_425, hk_426, hk_428, \
                         hk_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_0 * gk_315[k]
                   + pb_z[k] * hk_425[k];

        t_912[k] = f_7 * hi0_150[k]
                   - f_8 * hi1_150[k]
                   + pb_x[k] * hk_428[k];

        t_913[k] = pb_y[k] * hk_426[k];

        t_914[k] = f_7 * hi0_151[k]
                   - f_8 * hi1_151[k]
                   + pb_x[k] * hk_429[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pb_x, pb_z, gk_318, hi0_152, hi0_153, hi1_152, \
                         hi1_153, hk_427, hk_430, hk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_5 * hi0_152[k]
                   - f_6 * hi1_152[k]
                   + pb_x[k] * hk_430[k];

        t_916[k] = f_0 * gk_318[k]
                   + pb_z[k] * hk_427[k];

        t_917[k] = f_5 * hi0_153[k]
                   - f_6 * hi1_153[k]
                   + pb_x[k] * hk_431[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, pb_x, pb_y, hi0_154, hi0_155, hi0_156, \
                         hi1_154, hi1_155, hi1_156, hk_429, hk_432, hk_433, \
                         hk_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_5 * hi0_154[k]
                   - f_6 * hi1_154[k]
                   + pb_x[k] * hk_432[k];

        t_919[k] = pb_y[k] * hk_429[k];

        t_920[k] = f_5 * hi0_155[k]
                   - f_6 * hi1_155[k]
                   + pb_x[k] * hk_433[k];

        t_921[k] = f_3 * hi0_156[k]
                   - f_4 * hi1_156[k]
                   + pb_x[k] * hk_434[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pb_x, pb_z, gk_322, hi0_157, hi0_158, hi1_157, \
                         hi1_158, hk_430, hk_435, hk_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * gk_322[k]
                   + pb_z[k] * hk_430[k];

        t_923[k] = f_3 * hi0_157[k]
                   - f_4 * hi1_157[k]
                   + pb_x[k] * hk_435[k];

        t_924[k] = f_3 * hi0_158[k]
                   - f_4 * hi1_158[k]
                   + pb_x[k] * hk_436[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, pb_x, pb_y, hi0_159, hi0_161, \
                         hi1_159, hi1_161, hk_433, hk_437, hk_438, hk_439, \
                         hk_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_3 * hi0_159[k]
                   - f_4 * hi1_159[k]
                   + pb_x[k] * hk_437[k];

        t_926[k] = pb_y[k] * hk_433[k];

        t_927[k] = f_3 * hi0_161[k]
                   - f_4 * hi1_161[k]
                   + pb_x[k] * hk_438[k];

        t_928[k] = pb_x[k] * hk_439[k];

        t_929[k] = pb_x[k] * hk_440[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, t_935, pb_x, hk_441, hk_442, \
                         hk_443, hk_444, hk_445, hk_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = pb_x[k] * hk_441[k];

        t_931[k] = pb_x[k] * hk_442[k];

        t_932[k] = pb_x[k] * hk_443[k];

        t_933[k] = pb_x[k] * hk_444[k];

        t_934[k] = pb_x[k] * hk_445[k];

        t_935[k] = pb_x[k] * hk_446[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pb_y, pb_z, gk_332, hi0_156, hi0_157, \
                         hi0_158, hi1_156, hi1_157, hi1_158, hk_439, hk_441, \
                         hk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * hi0_156[k]
                   - f_2 * hi1_156[k]
                   + pb_y[k] * hk_439[k];

        t_937[k] = f_0 * gk_332[k]
                   + pb_z[k] * hk_439[k];

        t_938[k] = f_11 * hi0_157[k]
                   - f_12 * hi1_157[k]
                   + pb_y[k] * hk_441[k];

        t_939[k] = f_9 * hi0_158[k]
                   - f_10 * hi1_158[k]
                   + pb_y[k] * hk_442[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pb_y, hi0_159, hi0_160, hi0_161, hi1_159, \
                         hi1_160, hi1_161, hk_443, hk_444, hk_445, \
                         hk_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_7 * hi0_159[k]
                   - f_8 * hi1_159[k]
                   + pb_y[k] * hk_443[k];

        t_941[k] = f_5 * hi0_160[k]
                   - f_6 * hi1_160[k]
                   + pb_y[k] * hk_444[k];

        t_942[k] = f_3 * hi0_161[k]
                   - f_4 * hi1_161[k]
                   + pb_y[k] * hk_445[k];

        t_943[k] = pb_y[k] * hk_446[k];
    }

#pragma omp simd aligned(t_944, pb_z, gk_339, hi0_161, hi1_161, \
                         hk_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_0 * gk_339[k]
                   + f_1 * hi0_161[k]
                   - f_2 * hi1_161[k]
                   + pb_z[k] * hk_446[k];
    }
}

auto
compute_prim_hl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
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
    const auto f_16 = 2.0 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_73 = buffer.data(gk + 73);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
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
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
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
    const auto *gk_178 = buffer.data(gk + 178);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_181 = buffer.data(gk + 181);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_1 = buffer.data(hi0 + 1);
    const auto *hi0_2 = buffer.data(hi0 + 2);
    const auto *hi0_3 = buffer.data(hi0 + 3);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_6 = buffer.data(hi0 + 6);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_9 = buffer.data(hi0 + 9);
    const auto *hi0_10 = buffer.data(hi0 + 10);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_15 = buffer.data(hi0 + 15);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_21 = buffer.data(hi0 + 21);
    const auto *hi0_22 = buffer.data(hi0 + 22);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_25 = buffer.data(hi0 + 25);
    const auto *hi0_26 = buffer.data(hi0 + 26);
    const auto *hi0_27 = buffer.data(hi0 + 27);
    const auto *hi0_28 = buffer.data(hi0 + 28);
    const auto *hi0_29 = buffer.data(hi0 + 29);
    const auto *hi0_30 = buffer.data(hi0 + 30);
    const auto *hi0_31 = buffer.data(hi0 + 31);
    const auto *hi0_32 = buffer.data(hi0 + 32);
    const auto *hi0_33 = buffer.data(hi0 + 33);
    const auto *hi0_34 = buffer.data(hi0 + 34);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_36 = buffer.data(hi0 + 36);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_38 = buffer.data(hi0 + 38);
    const auto *hi0_39 = buffer.data(hi0 + 39);
    const auto *hi0_40 = buffer.data(hi0 + 40);
    const auto *hi0_41 = buffer.data(hi0 + 41);
    const auto *hi0_42 = buffer.data(hi0 + 42);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_44 = buffer.data(hi0 + 44);
    const auto *hi0_45 = buffer.data(hi0 + 45);
    const auto *hi0_46 = buffer.data(hi0 + 46);
    const auto *hi0_47 = buffer.data(hi0 + 47);
    const auto *hi0_48 = buffer.data(hi0 + 48);
    const auto *hi0_49 = buffer.data(hi0 + 49);
    const auto *hi0_50 = buffer.data(hi0 + 50);
    const auto *hi0_51 = buffer.data(hi0 + 51);
    const auto *hi0_52 = buffer.data(hi0 + 52);
    const auto *hi0_53 = buffer.data(hi0 + 53);
    const auto *hi0_54 = buffer.data(hi0 + 54);
    const auto *hi0_55 = buffer.data(hi0 + 55);
    const auto *hi0_56 = buffer.data(hi0 + 56);
    const auto *hi0_57 = buffer.data(hi0 + 57);
    const auto *hi0_58 = buffer.data(hi0 + 58);
    const auto *hi0_59 = buffer.data(hi0 + 59);
    const auto *hi0_60 = buffer.data(hi0 + 60);
    const auto *hi0_61 = buffer.data(hi0 + 61);
    const auto *hi0_62 = buffer.data(hi0 + 62);
    const auto *hi0_63 = buffer.data(hi0 + 63);
    const auto *hi0_64 = buffer.data(hi0 + 64);
    const auto *hi0_65 = buffer.data(hi0 + 65);
    const auto *hi0_66 = buffer.data(hi0 + 66);
    const auto *hi0_67 = buffer.data(hi0 + 67);
    const auto *hi0_68 = buffer.data(hi0 + 68);
    const auto *hi0_69 = buffer.data(hi0 + 69);
    const auto *hi0_70 = buffer.data(hi0 + 70);
    const auto *hi0_71 = buffer.data(hi0 + 71);
    const auto *hi0_72 = buffer.data(hi0 + 72);
    const auto *hi0_73 = buffer.data(hi0 + 73);
    const auto *hi0_74 = buffer.data(hi0 + 74);
    const auto *hi0_75 = buffer.data(hi0 + 75);
    const auto *hi0_76 = buffer.data(hi0 + 76);
    const auto *hi0_77 = buffer.data(hi0 + 77);
    const auto *hi0_78 = buffer.data(hi0 + 78);
    const auto *hi0_79 = buffer.data(hi0 + 79);
    const auto *hi0_80 = buffer.data(hi0 + 80);
    const auto *hi0_81 = buffer.data(hi0 + 81);
    const auto *hi0_82 = buffer.data(hi0 + 82);
    const auto *hi0_83 = buffer.data(hi0 + 83);
    const auto *hi0_84 = buffer.data(hi0 + 84);
    const auto *hi0_85 = buffer.data(hi0 + 85);
    const auto *hi0_86 = buffer.data(hi0 + 86);
    const auto *hi0_87 = buffer.data(hi0 + 87);
    const auto *hi0_88 = buffer.data(hi0 + 88);
    const auto *hi0_89 = buffer.data(hi0 + 89);
    const auto *hi0_90 = buffer.data(hi0 + 90);
    const auto *hi0_91 = buffer.data(hi0 + 91);
    const auto *hi0_97 = buffer.data(hi0 + 97);
    const auto *hi0_98 = buffer.data(hi0 + 98);
    const auto *hi0_99 = buffer.data(hi0 + 99);
    const auto *hi0_100 = buffer.data(hi0 + 100);
    const auto *hi0_101 = buffer.data(hi0 + 101);
    const auto *hi0_102 = buffer.data(hi0 + 102);
    const auto *hi0_103 = buffer.data(hi0 + 103);
    const auto *hi0_104 = buffer.data(hi0 + 104);
    const auto *hi0_105 = buffer.data(hi0 + 105);
    const auto *hi0_106 = buffer.data(hi0 + 106);
    const auto *hi0_107 = buffer.data(hi0 + 107);
    const auto *hi0_108 = buffer.data(hi0 + 108);
    const auto *hi0_109 = buffer.data(hi0 + 109);
    const auto *hi0_110 = buffer.data(hi0 + 110);
    const auto *hi0_111 = buffer.data(hi0 + 111);
    const auto *hi0_112 = buffer.data(hi0 + 112);
    const auto *hi0_113 = buffer.data(hi0 + 113);
    const auto *hi0_114 = buffer.data(hi0 + 114);
    const auto *hi0_116 = buffer.data(hi0 + 116);
    const auto *hi0_117 = buffer.data(hi0 + 117);
    const auto *hi0_118 = buffer.data(hi0 + 118);
    const auto *hi0_119 = buffer.data(hi0 + 119);
    const auto *hi0_120 = buffer.data(hi0 + 120);
    const auto *hi0_121 = buffer.data(hi0 + 121);
    const auto *hi0_122 = buffer.data(hi0 + 122);
    const auto *hi0_123 = buffer.data(hi0 + 123);
    const auto *hi0_124 = buffer.data(hi0 + 124);
    const auto *hi0_125 = buffer.data(hi0 + 125);
    const auto *hi0_126 = buffer.data(hi0 + 126);
    const auto *hi0_127 = buffer.data(hi0 + 127);
    const auto *hi0_128 = buffer.data(hi0 + 128);
    const auto *hi0_129 = buffer.data(hi0 + 129);
    const auto *hi0_130 = buffer.data(hi0 + 130);
    const auto *hi0_131 = buffer.data(hi0 + 131);
    const auto *hi0_132 = buffer.data(hi0 + 132);
    const auto *hi0_133 = buffer.data(hi0 + 133);
    const auto *hi0_134 = buffer.data(hi0 + 134);
    const auto *hi0_135 = buffer.data(hi0 + 135);
    const auto *hi0_136 = buffer.data(hi0 + 136);
    const auto *hi0_137 = buffer.data(hi0 + 137);
    const auto *hi0_138 = buffer.data(hi0 + 138);
    const auto *hi0_139 = buffer.data(hi0 + 139);
    const auto *hi0_140 = buffer.data(hi0 + 140);
    const auto *hi0_141 = buffer.data(hi0 + 141);
    const auto *hi0_142 = buffer.data(hi0 + 142);
    const auto *hi0_143 = buffer.data(hi0 + 143);
    const auto *hi0_144 = buffer.data(hi0 + 144);
    const auto *hi0_145 = buffer.data(hi0 + 145);
    const auto *hi0_146 = buffer.data(hi0 + 146);
    const auto *hi0_147 = buffer.data(hi0 + 147);
    const auto *hi0_148 = buffer.data(hi0 + 148);
    const auto *hi0_149 = buffer.data(hi0 + 149);
    const auto *hi0_150 = buffer.data(hi0 + 150);
    const auto *hi0_151 = buffer.data(hi0 + 151);
    const auto *hi0_153 = buffer.data(hi0 + 153);
    const auto *hi0_154 = buffer.data(hi0 + 154);
    const auto *hi0_155 = buffer.data(hi0 + 155);
    const auto *hi0_156 = buffer.data(hi0 + 156);
    const auto *hi0_157 = buffer.data(hi0 + 157);
    const auto *hi0_158 = buffer.data(hi0 + 158);
    const auto *hi0_159 = buffer.data(hi0 + 159);
    const auto *hi0_160 = buffer.data(hi0 + 160);
    const auto *hi0_161 = buffer.data(hi0 + 161);
    const auto *hi0_162 = buffer.data(hi0 + 162);
    const auto *hi0_163 = buffer.data(hi0 + 163);
    const auto *hi0_164 = buffer.data(hi0 + 164);
    const auto *hi0_165 = buffer.data(hi0 + 165);
    const auto *hi0_166 = buffer.data(hi0 + 166);
    const auto *hi0_167 = buffer.data(hi0 + 167);
    const auto *hi0_168 = buffer.data(hi0 + 168);
    const auto *hi0_169 = buffer.data(hi0 + 169);
    const auto *hi0_170 = buffer.data(hi0 + 170);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_1 = buffer.data(hi1 + 1);
    const auto *hi1_2 = buffer.data(hi1 + 2);
    const auto *hi1_3 = buffer.data(hi1 + 3);
    const auto *hi1_4 = buffer.data(hi1 + 4);
    const auto *hi1_5 = buffer.data(hi1 + 5);
    const auto *hi1_6 = buffer.data(hi1 + 6);
    const auto *hi1_7 = buffer.data(hi1 + 7);
    const auto *hi1_8 = buffer.data(hi1 + 8);
    const auto *hi1_9 = buffer.data(hi1 + 9);
    const auto *hi1_10 = buffer.data(hi1 + 10);
    const auto *hi1_11 = buffer.data(hi1 + 11);
    const auto *hi1_12 = buffer.data(hi1 + 12);
    const auto *hi1_14 = buffer.data(hi1 + 14);
    const auto *hi1_15 = buffer.data(hi1 + 15);
    const auto *hi1_16 = buffer.data(hi1 + 16);
    const auto *hi1_17 = buffer.data(hi1 + 17);
    const auto *hi1_18 = buffer.data(hi1 + 18);
    const auto *hi1_32 = buffer.data(hi1 + 32);
    const auto *hi1_33 = buffer.data(hi1 + 33);
    const auto *hi1_34 = buffer.data(hi1 + 34);
    const auto *hi1_35 = buffer.data(hi1 + 35);
    const auto *hi1_36 = buffer.data(hi1 + 36);
    const auto *hi1_37 = buffer.data(hi1 + 37);
    const auto *hi1_38 = buffer.data(hi1 + 38);
    const auto *hi1_39 = buffer.data(hi1 + 39);
    const auto *hi1_40 = buffer.data(hi1 + 40);
    const auto *hi1_41 = buffer.data(hi1 + 41);
    const auto *hi1_42 = buffer.data(hi1 + 42);
    const auto *hi1_43 = buffer.data(hi1 + 43);
    const auto *hi1_44 = buffer.data(hi1 + 44);
    const auto *hi1_45 = buffer.data(hi1 + 45);
    const auto *hi1_46 = buffer.data(hi1 + 46);
    const auto *hi1_47 = buffer.data(hi1 + 47);
    const auto *hi1_48 = buffer.data(hi1 + 48);
    const auto *hi1_49 = buffer.data(hi1 + 49);
    const auto *hi1_52 = buffer.data(hi1 + 52);
    const auto *hi1_53 = buffer.data(hi1 + 53);
    const auto *hi1_54 = buffer.data(hi1 + 54);
    const auto *hi1_55 = buffer.data(hi1 + 55);
    const auto *hi1_56 = buffer.data(hi1 + 56);
    const auto *hi1_57 = buffer.data(hi1 + 57);
    const auto *hi1_58 = buffer.data(hi1 + 58);
    const auto *hi1_59 = buffer.data(hi1 + 59);
    const auto *hi1_60 = buffer.data(hi1 + 60);
    const auto *hi1_61 = buffer.data(hi1 + 61);
    const auto *hi1_62 = buffer.data(hi1 + 62);
    const auto *hi1_63 = buffer.data(hi1 + 63);
    const auto *hi1_64 = buffer.data(hi1 + 64);
    const auto *hi1_65 = buffer.data(hi1 + 65);
    const auto *hi1_66 = buffer.data(hi1 + 66);
    const auto *hi1_67 = buffer.data(hi1 + 67);
    const auto *hi1_68 = buffer.data(hi1 + 68);
    const auto *hi1_69 = buffer.data(hi1 + 69);
    const auto *hi1_70 = buffer.data(hi1 + 70);
    const auto *hi1_71 = buffer.data(hi1 + 71);
    const auto *hi1_72 = buffer.data(hi1 + 72);
    const auto *hi1_73 = buffer.data(hi1 + 73);
    const auto *hi1_74 = buffer.data(hi1 + 74);
    const auto *hi1_75 = buffer.data(hi1 + 75);
    const auto *hi1_76 = buffer.data(hi1 + 76);
    const auto *hi1_77 = buffer.data(hi1 + 77);
    const auto *hi1_78 = buffer.data(hi1 + 78);
    const auto *hi1_79 = buffer.data(hi1 + 79);
    const auto *hi1_80 = buffer.data(hi1 + 80);
    const auto *hi1_81 = buffer.data(hi1 + 81);
    const auto *hi1_82 = buffer.data(hi1 + 82);
    const auto *hi1_83 = buffer.data(hi1 + 83);
    const auto *hi1_84 = buffer.data(hi1 + 84);
    const auto *hi1_85 = buffer.data(hi1 + 85);
    const auto *hi1_86 = buffer.data(hi1 + 86);
    const auto *hi1_87 = buffer.data(hi1 + 87);
    const auto *hi1_93 = buffer.data(hi1 + 93);
    const auto *hi1_94 = buffer.data(hi1 + 94);
    const auto *hi1_95 = buffer.data(hi1 + 95);
    const auto *hi1_96 = buffer.data(hi1 + 96);
    const auto *hi1_97 = buffer.data(hi1 + 97);
    const auto *hi1_98 = buffer.data(hi1 + 98);
    const auto *hi1_99 = buffer.data(hi1 + 99);
    const auto *hi1_100 = buffer.data(hi1 + 100);
    const auto *hi1_101 = buffer.data(hi1 + 101);
    const auto *hi1_102 = buffer.data(hi1 + 102);
    const auto *hi1_103 = buffer.data(hi1 + 103);
    const auto *hi1_104 = buffer.data(hi1 + 104);
    const auto *hi1_105 = buffer.data(hi1 + 105);
    const auto *hi1_106 = buffer.data(hi1 + 106);
    const auto *hi1_107 = buffer.data(hi1 + 107);
    const auto *hi1_108 = buffer.data(hi1 + 108);
    const auto *hi1_109 = buffer.data(hi1 + 109);
    const auto *hi1_110 = buffer.data(hi1 + 110);
    const auto *hi1_131 = buffer.data(hi1 + 131);
    const auto *hi1_133 = buffer.data(hi1 + 133);
    const auto *hi1_134 = buffer.data(hi1 + 134);
    const auto *hi1_135 = buffer.data(hi1 + 135);
    const auto *hi1_136 = buffer.data(hi1 + 136);
    const auto *hi1_137 = buffer.data(hi1 + 137);
    const auto *hi1_138 = buffer.data(hi1 + 138);
    const auto *hi1_139 = buffer.data(hi1 + 139);
    const auto *hi1_140 = buffer.data(hi1 + 140);
    const auto *hi1_141 = buffer.data(hi1 + 141);
    const auto *hi1_142 = buffer.data(hi1 + 142);
    const auto *hi1_143 = buffer.data(hi1 + 143);
    const auto *hi1_144 = buffer.data(hi1 + 144);
    const auto *hi1_145 = buffer.data(hi1 + 145);
    const auto *hi1_146 = buffer.data(hi1 + 146);
    const auto *hi1_147 = buffer.data(hi1 + 147);
    const auto *hi1_148 = buffer.data(hi1 + 148);
    const auto *hi1_149 = buffer.data(hi1 + 149);
    const auto *hi1_160 = buffer.data(hi1 + 160);
    const auto *hi1_161 = buffer.data(hi1 + 161);
    const auto *hi1_162 = buffer.data(hi1 + 162);
    const auto *hi1_163 = buffer.data(hi1 + 163);
    const auto *hi1_164 = buffer.data(hi1 + 164);
    const auto *hi1_165 = buffer.data(hi1 + 165);
    const auto *hi1_166 = buffer.data(hi1 + 166);
    const auto *hi1_167 = buffer.data(hi1 + 167);
    const auto *hi1_168 = buffer.data(hi1 + 168);
    const auto *hi1_169 = buffer.data(hi1 + 169);
    const auto *hi1_170 = buffer.data(hi1 + 170);
    const auto *hi1_171 = buffer.data(hi1 + 171);
    const auto *hi1_172 = buffer.data(hi1 + 172);
    const auto *hi1_173 = buffer.data(hi1 + 173);
    const auto *hi1_174 = buffer.data(hi1 + 174);
    const auto *hi1_175 = buffer.data(hi1 + 175);
    const auto *hi1_176 = buffer.data(hi1 + 176);
    const auto *hi1_177 = buffer.data(hi1 + 177);
    const auto *hi1_178 = buffer.data(hi1 + 178);
    const auto *hi1_179 = buffer.data(hi1 + 179);
    const auto *hi1_180 = buffer.data(hi1 + 180);
    const auto *hi1_181 = buffer.data(hi1 + 181);
    const auto *hi1_182 = buffer.data(hi1 + 182);
    const auto *hi1_183 = buffer.data(hi1 + 183);
    const auto *hi1_184 = buffer.data(hi1 + 184);
    const auto *hi1_185 = buffer.data(hi1 + 185);
    const auto *hi1_186 = buffer.data(hi1 + 186);
    const auto *hi1_187 = buffer.data(hi1 + 187);
    const auto *hi1_188 = buffer.data(hi1 + 188);
    const auto *hi1_189 = buffer.data(hi1 + 189);
    const auto *hi1_190 = buffer.data(hi1 + 190);
    const auto *hi1_191 = buffer.data(hi1 + 191);
    const auto *hi1_192 = buffer.data(hi1 + 192);
    const auto *hi1_193 = buffer.data(hi1 + 193);
    const auto *hi1_194 = buffer.data(hi1 + 194);
    const auto *hi1_195 = buffer.data(hi1 + 195);
    const auto *hi1_206 = buffer.data(hi1 + 206);
    const auto *hi1_208 = buffer.data(hi1 + 208);
    const auto *hi1_209 = buffer.data(hi1 + 209);
    const auto *hi1_210 = buffer.data(hi1 + 210);
    const auto *hi1_211 = buffer.data(hi1 + 211);
    const auto *hi1_212 = buffer.data(hi1 + 212);
    const auto *hi1_213 = buffer.data(hi1 + 213);
    const auto *hi1_214 = buffer.data(hi1 + 214);
    const auto *hi1_215 = buffer.data(hi1 + 215);
    const auto *hi1_216 = buffer.data(hi1 + 216);
    const auto *hi1_217 = buffer.data(hi1 + 217);
    const auto *hi1_218 = buffer.data(hi1 + 218);
    const auto *hi1_219 = buffer.data(hi1 + 219);
    const auto *hi1_220 = buffer.data(hi1 + 220);
    const auto *hi1_221 = buffer.data(hi1 + 221);
    const auto *hi1_222 = buffer.data(hi1 + 222);
    const auto *hi1_223 = buffer.data(hi1 + 223);
    const auto *hi1_224 = buffer.data(hi1 + 224);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_24 = buffer.data(hk + 24);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_40 = buffer.data(hk + 40);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_47 = buffer.data(hk + 47);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_95 = buffer.data(hk + 95);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);
    const auto *hk_135 = buffer.data(hk + 135);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_144 = buffer.data(hk + 144);
    const auto *hk_145 = buffer.data(hk + 145);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_148 = buffer.data(hk + 148);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_160 = buffer.data(hk + 160);
    const auto *hk_166 = buffer.data(hk + 166);
    const auto *hk_167 = buffer.data(hk + 167);
    const auto *hk_168 = buffer.data(hk + 168);
    const auto *hk_169 = buffer.data(hk + 169);
    const auto *hk_170 = buffer.data(hk + 170);
    const auto *hk_171 = buffer.data(hk + 171);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_184 = buffer.data(hk + 184);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_187 = buffer.data(hk + 187);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_196 = buffer.data(hk + 196);
    const auto *hk_203 = buffer.data(hk + 203);
    const auto *hk_204 = buffer.data(hk + 204);
    const auto *hk_205 = buffer.data(hk + 205);
    const auto *hk_206 = buffer.data(hk + 206);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_217 = buffer.data(hk + 217);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_220 = buffer.data(hk + 220);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_224 = buffer.data(hk + 224);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_229 = buffer.data(hk + 229);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_235 = buffer.data(hk + 235);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_238 = buffer.data(hk + 238);
    const auto *hk_239 = buffer.data(hk + 239);
    const auto *hk_240 = buffer.data(hk + 240);
    const auto *hk_241 = buffer.data(hk + 241);
    const auto *hk_242 = buffer.data(hk + 242);
    const auto *hk_243 = buffer.data(hk + 243);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_253 = buffer.data(hk + 253);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_265 = buffer.data(hk + 265);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_268 = buffer.data(hk + 268);
    const auto *hk_269 = buffer.data(hk + 269);
    const auto *hk_270 = buffer.data(hk + 270);
    const auto *hk_271 = buffer.data(hk + 271);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_273 = buffer.data(hk + 273);
    const auto *hk_274 = buffer.data(hk + 274);
    const auto *hk_275 = buffer.data(hk + 275);
    const auto *hk_276 = buffer.data(hk + 276);
    const auto *hk_277 = buffer.data(hk + 277);
    const auto *hk_278 = buffer.data(hk + 278);
    const auto *hk_279 = buffer.data(hk + 279);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_289 = buffer.data(hk + 289);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gk_0, hi0_0, hi0_1, hi1_0, \
                         hi1_1, hk_0, hk_1, hk_2, hk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = f_3 * hi0_0[k]
                 - f_4 * hi1_0[k]
                 + pb_y[k] * hk_1[k];

        t_2[k] = f_3 * hi0_0[k]
                 - f_4 * hi1_0[k]
                 + pb_z[k] * hk_2[k];

        t_3[k] = f_5 * hi0_1[k]
                 - f_6 * hi1_1[k]
                 + pb_y[k] * hk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, hi0_2, hi0_3, hi0_4, hi1_2, hi1_3, \
                         hi1_4, hk_4, hk_5, hk_6, hk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hi0_2[k]
                 - f_6 * hi1_2[k]
                 + pb_z[k] * hk_4[k];

        t_5[k] = f_7 * hi0_3[k]
                 - f_8 * hi1_3[k]
                 + pb_y[k] * hk_5[k];

        t_6[k] = f_3 * hi0_4[k]
                 - f_4 * hi1_4[k]
                 + pb_y[k] * hk_6[k];

        t_7[k] = f_7 * hi0_4[k]
                 - f_8 * hi1_4[k]
                 + pb_z[k] * hk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, hi0_5, hi0_6, hi0_7, hi1_5, hi1_6, \
                         hi1_7, hk_8, hk_9, hk_10, hk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * hi0_5[k]
                 - f_10 * hi1_5[k]
                 + pb_y[k] * hk_8[k];

        t_9[k] = f_5 * hi0_6[k]
                 - f_6 * hi1_6[k]
                 + pb_y[k] * hk_9[k];

        t_10[k] = f_3 * hi0_7[k]
                  - f_4 * hi1_7[k]
                  + pb_y[k] * hk_10[k];

        t_11[k] = f_9 * hi0_7[k]
                  - f_10 * hi1_7[k]
                  + pb_z[k] * hk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, hi0_8, hi0_9, hi0_10, hi1_8, hi1_9, hi1_10, \
                         hk_12, hk_13, hk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * hi0_8[k]
                  - f_12 * hi1_8[k]
                  + pb_y[k] * hk_12[k];

        t_13[k] = f_7 * hi0_9[k]
                  - f_8 * hi1_9[k]
                  + pb_y[k] * hk_13[k];

        t_14[k] = f_5 * hi0_10[k]
                  - f_6 * hi1_10[k]
                  + pb_y[k] * hk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gk_17, gk_23, hi0_11, \
                         hi1_11, hk_15, hk_16, hk_17, hk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * hi0_11[k]
                  - f_4 * hi1_11[k]
                  + pb_y[k] * hk_15[k];

        t_16[k] = f_11 * hi0_11[k]
                  - f_12 * hi1_11[k]
                  + pb_z[k] * hk_16[k];

        t_17[k] = f_0 * gk_17[k]
                  + pb_x[k] * hk_17[k];

        t_18[k] = f_0 * gk_23[k]
                  + pb_x[k] * hk_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, hi0_12, hi0_13, hi0_14, hi1_12, hi1_14, \
                         hi1_15, hk_17, hk_18, hk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * hi0_12[k]
                  - f_2 * hi1_12[k]
                  + pb_y[k] * hk_17[k];

        t_20[k] = f_11 * hi0_13[k]
                  - f_12 * hi1_14[k]
                  + pb_y[k] * hk_18[k];

        t_21[k] = f_9 * hi0_14[k]
                  - f_10 * hi1_15[k]
                  + pb_y[k] * hk_19[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_y, pb_z, hi0_15, hi0_16, hi0_17, hi1_16, \
                         hi1_17, hi1_18, hk_20, hk_21, hk_22, hk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hi0_15[k]
                  - f_8 * hi1_16[k]
                  + pb_y[k] * hk_20[k];

        t_23[k] = f_5 * hi0_16[k]
                  - f_6 * hi1_17[k]
                  + pb_y[k] * hk_21[k];

        t_24[k] = f_3 * hi0_17[k]
                  - f_4 * hi1_18[k]
                  + pb_y[k] * hk_22[k];

        t_25[k] = f_1 * hi0_17[k]
                  - f_2 * hi1_18[k]
                  + pb_z[k] * hk_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pb_y, gk_0, gk_1, gk_3, gk_5, \
                         gl_0, gl_1, gl_3, gl_5, hk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * gl_0[k];

        t_27[k] = f_13 * gk_0[k]
                  + pb_y[k] * hk_24[k];

        t_28[k] = f_14 * gk_1[k]
                  + pa_y[k] * gl_1[k];

        t_29[k] = f_15 * gk_3[k]
                  + pa_y[k] * gl_3[k];

        t_30[k] = f_16 * gk_5[k]
                  + pa_y[k] * gl_5[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pb_x, gk_8, gk_12, gk_17, gk_29, gl_8, \
                         gl_12, gl_17, hk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * gk_8[k]
                  + pa_y[k] * gl_8[k];

        t_32[k] = f_17 * gk_12[k]
                  + pa_y[k] * gl_12[k];

        t_33[k] = f_16 * gk_29[k]
                  + pb_x[k] * hk_29[k];

        t_34[k] = f_18 * gk_17[k]
                  + pa_y[k] * gl_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, pb_z, gk_0, gk_2, gk_4, gk_6, \
                         gl_0, gl_2, gl_4, gl_6, hk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * gl_0[k];

        t_36[k] = f_13 * gk_0[k]
                  + pb_z[k] * hk_30[k];

        t_37[k] = f_14 * gk_2[k]
                  + pa_z[k] * gl_2[k];

        t_38[k] = f_15 * gk_4[k]
                  + pa_z[k] * gl_4[k];

        t_39[k] = f_14 * gk_6[k]
                  + pa_z[k] * gl_6[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_z, gk_7, gk_9, gk_10, gk_11, gk_13, \
                         gl_7, gl_9, gl_10, gl_11, gl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_16 * gk_7[k]
                  + pa_z[k] * gl_7[k];

        t_41[k] = f_14 * gk_9[k]
                  + pa_z[k] * gl_9[k];

        t_42[k] = f_15 * gk_10[k]
                  + pa_z[k] * gl_10[k];

        t_43[k] = f_0 * gk_11[k]
                  + pa_z[k] * gl_11[k];

        t_44[k] = f_14 * gk_13[k]
                  + pa_z[k] * gl_13[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_z, pb_x, gk_14, gk_15, gk_16, gk_40, \
                         gl_14, gl_15, gl_16, hk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * gk_14[k]
                  + pa_z[k] * gl_14[k];

        t_46[k] = f_16 * gk_15[k]
                  + pa_z[k] * gl_15[k];

        t_47[k] = f_17 * gk_16[k]
                  + pa_z[k] * gl_16[k];

        t_48[k] = f_16 * gk_40[k]
                  + pb_x[k] * hk_40[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_z, gk_18, gk_19, gk_20, gk_21, \
                         gk_22, gl_18, gl_19, gl_20, gl_21, gl_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_14 * gk_18[k]
                  + pa_z[k] * gl_18[k];

        t_50[k] = f_15 * gk_19[k]
                  + pa_z[k] * gl_19[k];

        t_51[k] = f_16 * gk_20[k]
                  + pa_z[k] * gl_20[k];

        t_52[k] = f_0 * gk_21[k]
                  + pa_z[k] * gl_21[k];

        t_53[k] = f_17 * gk_22[k]
                  + pa_z[k] * gl_22[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pa_z, pb_y, fl0_0, fl1_0, gk_23, gk_24, \
                         gl_23, gl_24, hk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_18 * gk_23[k]
                  + pa_z[k] * gl_23[k];

        t_55[k] = f_19 * fl0_0[k]
                  - f_20 * fl1_0[k]
                  + pa_y[k] * gl_24[k];

        t_56[k] = f_14 * gk_24[k]
                  + pb_y[k] * hk_41[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, pb_z, gk_42, gk_44, hi0_20, hi0_22, hi0_24, \
                         hi1_32, hi1_34, hi1_36, hk_42, hk_43, hk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_15 * gk_42[k]
                  + f_11 * hi0_22[k]
                  - f_12 * hi1_34[k]
                  + pb_x[k] * hk_43[k];

        t_58[k] = f_3 * hi0_20[k]
                  - f_4 * hi1_32[k]
                  + pb_z[k] * hk_42[k];

        t_59[k] = f_15 * gk_44[k]
                  + f_9 * hi0_24[k]
                  - f_10 * hi1_36[k]
                  + pb_x[k] * hk_45[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_x, pb_z, gk_46, hi0_21, hi0_22, hi0_27, hi1_33, \
                         hi1_34, hi1_39, hk_44, hk_46, hk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * hi0_21[k]
                  - f_6 * hi1_33[k]
                  + pb_z[k] * hk_44[k];

        t_61[k] = f_15 * gk_46[k]
                  + f_7 * hi0_27[k]
                  - f_8 * hi1_39[k]
                  + pb_x[k] * hk_48[k];

        t_62[k] = f_3 * hi0_22[k]
                  - f_4 * hi1_34[k]
                  + pb_z[k] * hk_46[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, pb_z, gk_48, hi0_23, hi0_24, hi0_31, hi1_35, \
                         hi1_36, hi1_43, hk_47, hk_49, hk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * hi0_23[k]
                  - f_8 * hi1_35[k]
                  + pb_z[k] * hk_47[k];

        t_64[k] = f_15 * gk_48[k]
                  + f_5 * hi0_31[k]
                  - f_6 * hi1_43[k]
                  + pb_x[k] * hk_52[k];

        t_65[k] = f_3 * hi0_24[k]
                  - f_4 * hi1_36[k]
                  + pb_z[k] * hk_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, gk_50, hi0_25, hi0_26, hi0_32, hi1_37, \
                         hi1_38, hi1_44, hk_50, hk_51, hk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * hi0_25[k]
                  - f_6 * hi1_37[k]
                  + pb_z[k] * hk_50[k];

        t_67[k] = f_9 * hi0_26[k]
                  - f_10 * hi1_38[k]
                  + pb_z[k] * hk_51[k];

        t_68[k] = f_15 * gk_50[k]
                  + f_3 * hi0_32[k]
                  - f_4 * hi1_44[k]
                  + pb_x[k] * hk_57[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, hi0_27, hi0_28, hi0_29, hi1_39, hi1_40, \
                         hi1_41, hk_53, hk_54, hk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * hi0_27[k]
                  - f_4 * hi1_39[k]
                  + pb_z[k] * hk_53[k];

        t_70[k] = f_5 * hi0_28[k]
                  - f_6 * hi1_40[k]
                  + pb_z[k] * hk_54[k];

        t_71[k] = f_7 * hi0_29[k]
                  - f_8 * hi1_41[k]
                  + pb_z[k] * hk_55[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_x, pb_x, pb_z, fl0_3, fl1_3, gk_51, gl_32, \
                         hi0_30, hi1_42, hk_56, hk_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * hi0_30[k]
                  - f_12 * hi1_42[k]
                  + pb_z[k] * hk_56[k];

        t_73[k] = f_15 * gk_51[k]
                  + pb_x[k] * hk_58[k];

        t_74[k] = f_21 * fl0_3[k]
                  - f_22 * fl1_3[k]
                  + pa_x[k] * gl_32[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, hi0_32, hi0_33, hi0_34, hi1_44, hi1_45, \
                         hi1_46, hk_59, hk_60, hk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * hi0_32[k]
                  - f_4 * hi1_44[k]
                  + pb_z[k] * hk_59[k];

        t_76[k] = f_5 * hi0_33[k]
                  - f_6 * hi1_45[k]
                  + pb_z[k] * hk_60[k];

        t_77[k] = f_7 * hi0_34[k]
                  - f_8 * hi1_46[k]
                  + pb_z[k] * hk_61[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_z, hi0_35, hi0_36, hi0_37, hi1_47, hi1_48, \
                         hi1_49, hk_62, hk_63, hk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * hi0_35[k]
                  - f_10 * hi1_47[k]
                  + pb_z[k] * hk_62[k];

        t_79[k] = f_11 * hi0_36[k]
                  - f_12 * hi1_48[k]
                  + pb_z[k] * hk_63[k];

        t_80[k] = f_1 * hi0_37[k]
                  - f_2 * hi1_49[k]
                  + pb_z[k] * hk_64[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, pb_z, fl0_0, fl1_0, gk_30, gl_25, \
                         hi0_38, hi1_52, hk_65, hk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_19 * fl0_0[k]
                  - f_20 * fl1_0[k]
                  + pa_z[k] * gl_25[k];

        t_82[k] = f_14 * gk_30[k]
                  + pb_z[k] * hk_65[k];

        t_83[k] = f_3 * hi0_38[k]
                  - f_4 * hi1_52[k]
                  + pb_y[k] * hk_66[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, gk_60, gk_62, hi0_39, hi0_41, hi0_44, \
                         hi1_53, hi1_55, hi1_58, hk_68, hk_69, hk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_15 * gk_60[k]
                  + f_11 * hi0_41[k]
                  - f_12 * hi1_55[k]
                  + pb_x[k] * hk_69[k];

        t_85[k] = f_5 * hi0_39[k]
                  - f_6 * hi1_53[k]
                  + pb_y[k] * hk_68[k];

        t_86[k] = f_15 * gk_62[k]
                  + f_9 * hi0_44[k]
                  - f_10 * hi1_58[k]
                  + pb_x[k] * hk_72[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, gk_64, hi0_40, hi0_41, hi0_48, hi1_54, \
                         hi1_55, hi1_62, hk_70, hk_71, hk_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_7 * hi0_40[k]
                  - f_8 * hi1_54[k]
                  + pb_y[k] * hk_70[k];

        t_88[k] = f_3 * hi0_41[k]
                  - f_4 * hi1_55[k]
                  + pb_y[k] * hk_71[k];

        t_89[k] = f_15 * gk_64[k]
                  + f_7 * hi0_48[k]
                  - f_8 * hi1_62[k]
                  + pb_x[k] * hk_76[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_y, hi0_42, hi0_43, hi0_44, hi1_56, hi1_57, \
                         hi1_58, hk_73, hk_74, hk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * hi0_42[k]
                  - f_10 * hi1_56[k]
                  + pb_y[k] * hk_73[k];

        t_91[k] = f_5 * hi0_43[k]
                  - f_6 * hi1_57[k]
                  + pb_y[k] * hk_74[k];

        t_92[k] = f_3 * hi0_44[k]
                  - f_4 * hi1_58[k]
                  + pb_y[k] * hk_75[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, pb_y, gk_66, hi0_45, hi0_46, hi0_49, hi1_59, \
                         hi1_60, hi1_63, hk_77, hk_78, hk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_15 * gk_66[k]
                  + f_5 * hi0_49[k]
                  - f_6 * hi1_63[k]
                  + pb_x[k] * hk_81[k];

        t_94[k] = f_11 * hi0_45[k]
                  - f_12 * hi1_59[k]
                  + pb_y[k] * hk_77[k];

        t_95[k] = f_7 * hi0_46[k]
                  - f_8 * hi1_60[k]
                  + pb_y[k] * hk_78[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_x, pb_y, gk_67, hi0_47, hi0_48, hi0_55, hi1_61, \
                         hi1_62, hi1_69, hk_79, hk_80, hk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * hi0_47[k]
                  - f_6 * hi1_61[k]
                  + pb_y[k] * hk_79[k];

        t_97[k] = f_3 * hi0_48[k]
                  - f_4 * hi1_62[k]
                  + pb_y[k] * hk_80[k];

        t_98[k] = f_15 * gk_67[k]
                  + f_3 * hi0_55[k]
                  - f_4 * hi1_69[k]
                  + pb_x[k] * hk_82[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_y, gk_73, hi0_50, hi0_51, hi1_64, \
                         hi1_65, hk_83, hk_84, hk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_15 * gk_73[k]
                  + pb_x[k] * hk_89[k];

        t_100[k] = f_1 * hi0_50[k]
                   - f_2 * hi1_64[k]
                   + pb_y[k] * hk_83[k];

        t_101[k] = f_11 * hi0_51[k]
                   - f_12 * hi1_65[k]
                   + pb_y[k] * hk_84[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_y, hi0_52, hi0_53, hi0_54, hi1_66, hi1_67, \
                         hi1_68, hk_85, hk_86, hk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * hi0_52[k]
                   - f_10 * hi1_66[k]
                   + pb_y[k] * hk_85[k];

        t_103[k] = f_7 * hi0_53[k]
                   - f_8 * hi1_67[k]
                   + pb_y[k] * hk_86[k];

        t_104[k] = f_5 * hi0_54[k]
                   - f_6 * hi1_68[k]
                   + pb_y[k] * hk_87[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pa_y, pb_y, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gl_26, gl_39, hi0_55, hi1_69, hk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * hi0_55[k]
                   - f_4 * hi1_69[k]
                   + pb_y[k] * hk_88[k];

        t_106[k] = f_21 * fl0_4[k]
                   - f_22 * fl1_4[k]
                   + pa_x[k] * gl_39[k];

        t_107[k] = f_21 * fl0_1[k]
                   - f_22 * fl1_1[k]
                   + pa_y[k] * gl_26[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_x, pb_y, pb_z, gk_41, gk_75, hi0_56, hi0_58, \
                         hi1_70, hi1_72, hk_90, hk_91, hk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * gk_41[k]
                   + pb_y[k] * hk_90[k];

        t_109[k] = f_14 * gk_75[k]
                   + f_11 * hi0_58[k]
                   - f_12 * hi1_72[k]
                   + pb_x[k] * hk_92[k];

        t_110[k] = f_3 * hi0_56[k]
                   - f_4 * hi1_70[k]
                   + pb_z[k] * hk_91[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pb_z, gk_76, gk_77, hi0_57, hi0_60, \
                         hi0_63, hi1_71, hi1_74, hi1_77, hk_93, hk_94, \
                         hk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_14 * gk_76[k]
                   + f_9 * hi0_60[k]
                   - f_10 * hi1_74[k]
                   + pb_x[k] * hk_94[k];

        t_112[k] = f_5 * hi0_57[k]
                   - f_6 * hi1_71[k]
                   + pb_z[k] * hk_93[k];

        t_113[k] = f_14 * gk_77[k]
                   + f_7 * hi0_63[k]
                   - f_8 * hi1_77[k]
                   + pb_x[k] * hk_97[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, pb_z, gk_78, hi0_58, hi0_59, hi0_67, \
                         hi1_72, hi1_73, hi1_81, hk_95, hk_96, hk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * hi0_58[k]
                   - f_4 * hi1_72[k]
                   + pb_z[k] * hk_95[k];

        t_115[k] = f_7 * hi0_59[k]
                   - f_8 * hi1_73[k]
                   + pb_z[k] * hk_96[k];

        t_116[k] = f_14 * gk_78[k]
                   + f_5 * hi0_67[k]
                   - f_6 * hi1_81[k]
                   + pb_x[k] * hk_101[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_z, hi0_60, hi0_61, hi0_62, hi1_74, hi1_75, \
                         hi1_76, hk_98, hk_99, hk_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_3 * hi0_60[k]
                   - f_4 * hi1_74[k]
                   + pb_z[k] * hk_98[k];

        t_118[k] = f_5 * hi0_61[k]
                   - f_6 * hi1_75[k]
                   + pb_z[k] * hk_99[k];

        t_119[k] = f_9 * hi0_62[k]
                   - f_10 * hi1_76[k]
                   + pb_z[k] * hk_100[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pb_x, pb_z, gk_79, hi0_63, hi0_64, hi0_68, \
                         hi1_77, hi1_78, hi1_82, hk_102, hk_103, \
                         hk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_14 * gk_79[k]
                   + f_3 * hi0_68[k]
                   - f_4 * hi1_82[k]
                   + pb_x[k] * hk_106[k];

        t_121[k] = f_3 * hi0_63[k]
                   - f_4 * hi1_77[k]
                   + pb_z[k] * hk_102[k];

        t_122[k] = f_5 * hi0_64[k]
                   - f_6 * hi1_78[k]
                   + pb_z[k] * hk_103[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_x, pb_z, gk_80, hi0_65, hi0_66, hi1_79, \
                         hi1_80, hk_104, hk_105, hk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_7 * hi0_65[k]
                   - f_8 * hi1_79[k]
                   + pb_z[k] * hk_104[k];

        t_124[k] = f_11 * hi0_66[k]
                   - f_12 * hi1_80[k]
                   + pb_z[k] * hk_105[k];

        t_125[k] = f_14 * gk_80[k]
                   + pb_x[k] * hk_107[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_z, fl0_5, fl1_5, gl_40, hi0_68, hi0_69, \
                         hi1_82, hi1_83, hk_108, hk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_19 * fl0_5[k]
                   - f_20 * fl1_5[k]
                   + pa_x[k] * gl_40[k];

        t_127[k] = f_3 * hi0_68[k]
                   - f_4 * hi1_82[k]
                   + pb_z[k] * hk_108[k];

        t_128[k] = f_5 * hi0_69[k]
                   - f_6 * hi1_83[k]
                   + pb_z[k] * hk_109[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pb_z, hi0_70, hi0_71, hi0_72, hi1_84, hi1_85, \
                         hi1_86, hk_110, hk_111, hk_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_7 * hi0_70[k]
                   - f_8 * hi1_84[k]
                   + pb_z[k] * hk_110[k];

        t_130[k] = f_9 * hi0_71[k]
                   - f_10 * hi1_85[k]
                   + pb_z[k] * hk_111[k];

        t_131[k] = f_11 * hi0_72[k]
                   - f_12 * hi1_86[k]
                   + pb_z[k] * hk_112[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, pa_z, pb_z, gl_27, gl_28, \
                         gl_29, gl_30, gl_31, hi0_73, hi1_87, hk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * hi0_73[k]
                   - f_2 * hi1_87[k]
                   + pb_z[k] * hk_113[k];

        t_133[k] = pa_z[k] * gl_27[k];

        t_134[k] = pa_z[k] * gl_28[k];

        t_135[k] = pa_z[k] * gl_29[k];

        t_136[k] = pa_z[k] * gl_30[k];

        t_137[k] = pa_z[k] * gl_31[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, pa_y, gl_33, gl_34, gl_35, \
                         gl_36, gl_37, gl_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pa_y[k] * gl_33[k];

        t_139[k] = pa_y[k] * gl_34[k];

        t_140[k] = pa_y[k] * gl_35[k];

        t_141[k] = pa_y[k] * gl_36[k];

        t_142[k] = pa_y[k] * gl_37[k];

        t_143[k] = pa_y[k] * gl_38[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pa_z, pb_y, pb_z, fl0_2, fl1_2, gk_57, gl_33, \
                         hi0_74, hi1_93, hk_123, hk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_21 * fl0_2[k]
                   - f_22 * fl1_2[k]
                   + pa_z[k] * gl_33[k];

        t_145[k] = f_15 * gk_57[k]
                   + pb_z[k] * hk_123[k];

        t_146[k] = f_3 * hi0_74[k]
                   - f_4 * hi1_93[k]
                   + pb_y[k] * hk_124[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_x, pb_y, gk_83, gk_84, hi0_75, hi0_77, \
                         hi0_80, hi1_94, hi1_96, hi1_99, hk_126, hk_127, \
                         hk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_14 * gk_83[k]
                   + f_11 * hi0_77[k]
                   - f_12 * hi1_96[k]
                   + pb_x[k] * hk_127[k];

        t_148[k] = f_5 * hi0_75[k]
                   - f_6 * hi1_94[k]
                   + pb_y[k] * hk_126[k];

        t_149[k] = f_14 * gk_84[k]
                   + f_9 * hi0_80[k]
                   - f_10 * hi1_99[k]
                   + pb_x[k] * hk_130[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_x, pb_y, gk_85, hi0_76, hi0_77, hi0_84, \
                         hi1_95, hi1_96, hi1_103, hk_128, hk_129, \
                         hk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_7 * hi0_76[k]
                   - f_8 * hi1_95[k]
                   + pb_y[k] * hk_128[k];

        t_151[k] = f_3 * hi0_77[k]
                   - f_4 * hi1_96[k]
                   + pb_y[k] * hk_129[k];

        t_152[k] = f_14 * gk_85[k]
                   + f_7 * hi0_84[k]
                   - f_8 * hi1_103[k]
                   + pb_x[k] * hk_134[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_y, hi0_78, hi0_79, hi0_80, hi1_97, hi1_98, \
                         hi1_99, hk_131, hk_132, hk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_9 * hi0_78[k]
                   - f_10 * hi1_97[k]
                   + pb_y[k] * hk_131[k];

        t_154[k] = f_5 * hi0_79[k]
                   - f_6 * hi1_98[k]
                   + pb_y[k] * hk_132[k];

        t_155[k] = f_3 * hi0_80[k]
                   - f_4 * hi1_99[k]
                   + pb_y[k] * hk_133[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_y, gk_86, hi0_81, hi0_82, hi0_85, \
                         hi1_100, hi1_101, hi1_104, hk_135, hk_136, \
                         hk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_14 * gk_86[k]
                   + f_5 * hi0_85[k]
                   - f_6 * hi1_104[k]
                   + pb_x[k] * hk_139[k];

        t_157[k] = f_11 * hi0_81[k]
                   - f_12 * hi1_100[k]
                   + pb_y[k] * hk_135[k];

        t_158[k] = f_7 * hi0_82[k]
                   - f_8 * hi1_101[k]
                   + pb_y[k] * hk_136[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_x, pb_y, gk_87, hi0_83, hi0_84, hi0_91, \
                         hi1_102, hi1_103, hi1_110, hk_137, hk_138, \
                         hk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * hi0_83[k]
                   - f_6 * hi1_102[k]
                   + pb_y[k] * hk_137[k];

        t_160[k] = f_3 * hi0_84[k]
                   - f_4 * hi1_103[k]
                   + pb_y[k] * hk_138[k];

        t_161[k] = f_14 * gk_87[k]
                   + f_3 * hi0_91[k]
                   - f_4 * hi1_110[k]
                   + pb_x[k] * hk_140[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_x, pb_y, gk_88, hi0_86, hi0_87, hi1_105, \
                         hi1_106, hk_141, hk_142, hk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_14 * gk_88[k]
                   + pb_x[k] * hk_147[k];

        t_163[k] = f_1 * hi0_86[k]
                   - f_2 * hi1_105[k]
                   + pb_y[k] * hk_141[k];

        t_164[k] = f_11 * hi0_87[k]
                   - f_12 * hi1_106[k]
                   + pb_y[k] * hk_142[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pb_y, hi0_88, hi0_89, hi0_90, hi1_107, hi1_108, \
                         hi1_109, hk_143, hk_144, hk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_9 * hi0_88[k]
                   - f_10 * hi1_107[k]
                   + pb_y[k] * hk_143[k];

        t_166[k] = f_7 * hi0_89[k]
                   - f_8 * hi1_108[k]
                   + pb_y[k] * hk_144[k];

        t_167[k] = f_5 * hi0_90[k]
                   - f_6 * hi1_109[k]
                   + pb_y[k] * hk_145[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pb_y, fl0_8, fl1_8, gk_74, gk_89, \
                         gl_41, gl_42, hi0_91, hi1_110, hk_146, \
                         hk_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_3 * hi0_91[k]
                   - f_4 * hi1_110[k]
                   + pb_y[k] * hk_146[k];

        t_169[k] = f_19 * fl0_8[k]
                   - f_20 * fl1_8[k]
                   + pa_x[k] * gl_41[k];

        t_170[k] = f_18 * gk_89[k]
                   + pa_x[k] * gl_42[k];

        t_171[k] = f_16 * gk_74[k]
                   + pb_y[k] * hk_148[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, pa_x, gk_91, gk_93, gk_96, gk_100, \
                         gk_105, gl_43, gl_45, gl_47, gl_50, gl_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_17 * gk_91[k]
                   + pa_x[k] * gl_43[k];

        t_173[k] = f_0 * gk_93[k]
                   + pa_x[k] * gl_45[k];

        t_174[k] = f_16 * gk_96[k]
                   + pa_x[k] * gl_47[k];

        t_175[k] = f_15 * gk_100[k]
                   + pa_x[k] * gl_50[k];

        t_176[k] = f_14 * gk_105[k]
                   + pa_x[k] * gl_54[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, pa_x, pb_x, gk_106, gl_59, \
                         gl_67, gl_68, gl_69, gl_70, hk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * gk_106[k]
                   + pb_x[k] * hk_153[k];

        t_178[k] = pa_x[k] * gl_59[k];

        t_179[k] = pa_x[k] * gl_67[k];

        t_180[k] = pa_x[k] * gl_68[k];

        t_181[k] = pa_x[k] * gl_69[k];

        t_182[k] = pa_x[k] * gl_70[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pa_x, pb_z, gk_81, gk_158, gl_71, \
                         gl_72, gl_73, gl_75, hk_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_x[k] * gl_71[k];

        t_184[k] = pa_x[k] * gl_72[k];

        t_185[k] = pa_x[k] * gl_73[k];

        t_186[k] = f_18 * gk_158[k]
                   + pa_x[k] * gl_75[k];

        t_187[k] = f_16 * gk_81[k]
                   + pb_z[k] * hk_160[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pa_x, gk_162, gk_165, gk_169, \
                         gk_174, gk_175, gl_77, gl_79, gl_82, gl_86, \
                         gl_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_17 * gk_162[k]
                   + pa_x[k] * gl_77[k];

        t_189[k] = f_0 * gk_165[k]
                   + pa_x[k] * gl_79[k];

        t_190[k] = f_16 * gk_169[k]
                   + pa_x[k] * gl_82[k];

        t_191[k] = f_15 * gk_174[k]
                   + pa_x[k] * gl_86[k];

        t_192[k] = f_14 * gk_175[k]
                   + pa_x[k] * gl_91[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_x, pb_x, pb_y, gk_89, gk_183, gl_98, \
                         hi0_97, hi1_131, hk_166, hk_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_13 * gk_183[k]
                   + pb_x[k] * hk_166[k];

        t_194[k] = pa_x[k] * gl_98[k];

        t_195[k] = f_1 * hi0_97[k]
                   - f_2 * hi1_131[k]
                   + pb_x[k] * hk_167[k];

        t_196[k] = f_0 * gk_89[k]
                   + pb_y[k] * hk_167[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_x, hi0_98, hi0_99, hi0_100, hi1_133, hi1_134, \
                         hi1_135, hk_168, hk_169, hk_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_11 * hi0_98[k]
                   - f_12 * hi1_133[k]
                   + pb_x[k] * hk_168[k];

        t_198[k] = f_11 * hi0_99[k]
                   - f_12 * hi1_134[k]
                   + pb_x[k] * hk_169[k];

        t_199[k] = f_9 * hi0_100[k]
                   - f_10 * hi1_135[k]
                   + pb_x[k] * hk_170[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, hi0_101, hi0_102, hi0_103, hi1_136, \
                         hi1_137, hi1_138, hk_171, hk_172, hk_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_9 * hi0_101[k]
                   - f_10 * hi1_136[k]
                   + pb_x[k] * hk_171[k];

        t_201[k] = f_7 * hi0_102[k]
                   - f_8 * hi1_137[k]
                   + pb_x[k] * hk_172[k];

        t_202[k] = f_7 * hi0_103[k]
                   - f_8 * hi1_138[k]
                   + pb_x[k] * hk_173[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pb_x, hi0_104, hi0_105, hi0_106, hi1_139, \
                         hi1_140, hi1_141, hk_174, hk_175, hk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_7 * hi0_104[k]
                   - f_8 * hi1_139[k]
                   + pb_x[k] * hk_174[k];

        t_204[k] = f_5 * hi0_105[k]
                   - f_6 * hi1_140[k]
                   + pb_x[k] * hk_175[k];

        t_205[k] = f_5 * hi0_106[k]
                   - f_6 * hi1_141[k]
                   + pb_x[k] * hk_176[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pb_x, hi0_107, hi0_108, hi0_109, hi1_142, \
                         hi1_143, hi1_144, hk_177, hk_178, hk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_5 * hi0_107[k]
                   - f_6 * hi1_142[k]
                   + pb_x[k] * hk_177[k];

        t_207[k] = f_5 * hi0_108[k]
                   - f_6 * hi1_143[k]
                   + pb_x[k] * hk_178[k];

        t_208[k] = f_3 * hi0_109[k]
                   - f_4 * hi1_144[k]
                   + pb_x[k] * hk_179[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_x, hi0_111, hi0_112, hi0_113, hi1_146, \
                         hi1_147, hi1_148, hk_180, hk_181, hk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_3 * hi0_111[k]
                   - f_4 * hi1_146[k]
                   + pb_x[k] * hk_180[k];

        t_210[k] = f_3 * hi0_112[k]
                   - f_4 * hi1_147[k]
                   + pb_x[k] * hk_181[k];

        t_211[k] = f_3 * hi0_113[k]
                   - f_4 * hi1_148[k]
                   + pb_x[k] * hk_182[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pb_x, pb_y, pb_z, gk_106, hi0_109, hi0_114, \
                         hi1_144, hi1_149, hk_183, hk_184, hk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * hi0_114[k]
                   - f_4 * hi1_149[k]
                   + pb_x[k] * hk_183[k];

        t_213[k] = f_0 * gk_106[k]
                   + f_1 * hi0_109[k]
                   - f_2 * hi1_144[k]
                   + pb_y[k] * hk_184[k];

        t_214[k] = f_3 * hi0_109[k]
                   - f_4 * hi1_144[k]
                   + pb_z[k] * hk_185[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pb_z, hi0_110, hi0_111, hi0_112, hi1_145, \
                         hi1_146, hi1_147, hk_186, hk_187, hk_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_5 * hi0_110[k]
                   - f_6 * hi1_145[k]
                   + pb_z[k] * hk_186[k];

        t_216[k] = f_7 * hi0_111[k]
                   - f_8 * hi1_146[k]
                   + pb_z[k] * hk_187[k];

        t_217[k] = f_9 * hi0_112[k]
                   - f_10 * hi1_147[k]
                   + pb_z[k] * hk_188[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pb_y, pb_z, gk_90, gk_113, gl_44, \
                         hi0_113, hi0_114, hi1_148, hi1_149, hk_189, \
                         hk_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * hi0_113[k]
                   - f_12 * hi1_148[k]
                   + pb_z[k] * hk_189[k];

        t_219[k] = f_0 * gk_113[k]
                   + pb_y[k] * hk_191[k];

        t_220[k] = f_1 * hi0_114[k]
                   - f_2 * hi1_149[k]
                   + pb_z[k] * hk_191[k];

        t_221[k] = f_14 * gk_90[k]
                   + pa_z[k] * gl_44[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pa_z, gk_92, gk_94, gk_95, gk_97, \
                         gk_98, gl_46, gl_48, gl_49, gl_51, gl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_15 * gk_92[k]
                   + pa_z[k] * gl_46[k];

        t_223[k] = f_14 * gk_94[k]
                   + pa_z[k] * gl_48[k];

        t_224[k] = f_16 * gk_95[k]
                   + pa_z[k] * gl_49[k];

        t_225[k] = f_14 * gk_97[k]
                   + pa_z[k] * gl_51[k];

        t_226[k] = f_15 * gk_98[k]
                   + pa_z[k] * gl_52[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, pa_z, gk_99, gk_101, gk_102, \
                         gk_103, gk_104, gl_53, gl_55, gl_56, gl_57, \
                         gl_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_0 * gk_99[k]
                   + pa_z[k] * gl_53[k];

        t_228[k] = f_14 * gk_101[k]
                   + pa_z[k] * gl_55[k];

        t_229[k] = f_15 * gk_102[k]
                   + pa_z[k] * gl_56[k];

        t_230[k] = f_16 * gk_103[k]
                   + pa_z[k] * gl_57[k];

        t_231[k] = f_17 * gk_104[k]
                   + pa_z[k] * gl_58[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pa_z, pb_z, gk_106, gk_107, \
                         gk_108, gk_109, gl_59, gl_60, gl_61, gl_62, \
                         hk_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pa_z[k] * gl_59[k];

        t_233[k] = f_13 * gk_106[k]
                   + pb_z[k] * hk_196[k];

        t_234[k] = f_14 * gk_107[k]
                   + pa_z[k] * gl_60[k];

        t_235[k] = f_15 * gk_108[k]
                   + pa_z[k] * gl_61[k];

        t_236[k] = f_16 * gk_109[k]
                   + pa_z[k] * gl_62[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_z, pb_y, gk_110, gk_111, gk_113, \
                         gk_125, gl_63, gl_64, gl_65, hk_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_0 * gk_110[k]
                   + pa_z[k] * gl_63[k];

        t_238[k] = f_17 * gk_111[k]
                   + pa_z[k] * gl_64[k];

        t_239[k] = f_16 * gk_125[k]
                   + pb_y[k] * hk_203[k];

        t_240[k] = f_18 * gk_113[k]
                   + pa_z[k] * gl_65[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, pb_x, hi0_116, hi0_117, hi0_118, hi1_160, \
                         hi1_161, hi1_162, hk_204, hk_205, hk_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_1 * hi0_116[k]
                   - f_2 * hi1_160[k]
                   + pb_x[k] * hk_204[k];

        t_242[k] = f_11 * hi0_117[k]
                   - f_12 * hi1_161[k]
                   + pb_x[k] * hk_205[k];

        t_243[k] = f_11 * hi0_118[k]
                   - f_12 * hi1_162[k]
                   + pb_x[k] * hk_206[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pb_x, hi0_119, hi0_120, hi0_121, hi1_163, \
                         hi1_164, hi1_165, hk_207, hk_208, hk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_9 * hi0_119[k]
                   - f_10 * hi1_163[k]
                   + pb_x[k] * hk_207[k];

        t_245[k] = f_9 * hi0_120[k]
                   - f_10 * hi1_164[k]
                   + pb_x[k] * hk_208[k];

        t_246[k] = f_7 * hi0_121[k]
                   - f_8 * hi1_165[k]
                   + pb_x[k] * hk_209[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pb_x, hi0_122, hi0_123, hi0_124, hi1_166, \
                         hi1_167, hi1_168, hk_210, hk_211, hk_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_7 * hi0_122[k]
                   - f_8 * hi1_166[k]
                   + pb_x[k] * hk_210[k];

        t_248[k] = f_7 * hi0_123[k]
                   - f_8 * hi1_167[k]
                   + pb_x[k] * hk_211[k];

        t_249[k] = f_5 * hi0_124[k]
                   - f_6 * hi1_168[k]
                   + pb_x[k] * hk_212[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pb_x, hi0_125, hi0_126, hi0_127, hi1_169, \
                         hi1_170, hi1_171, hk_213, hk_214, hk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_5 * hi0_125[k]
                   - f_6 * hi1_169[k]
                   + pb_x[k] * hk_213[k];

        t_251[k] = f_5 * hi0_126[k]
                   - f_6 * hi1_170[k]
                   + pb_x[k] * hk_214[k];

        t_252[k] = f_5 * hi0_127[k]
                   - f_6 * hi1_171[k]
                   + pb_x[k] * hk_215[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, pb_x, hi0_128, hi0_129, hi0_130, hi1_172, \
                         hi1_173, hi1_174, hk_216, hk_217, hk_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_3 * hi0_128[k]
                   - f_4 * hi1_172[k]
                   + pb_x[k] * hk_216[k];

        t_254[k] = f_3 * hi0_129[k]
                   - f_4 * hi1_173[k]
                   + pb_x[k] * hk_217[k];

        t_255[k] = f_3 * hi0_130[k]
                   - f_4 * hi1_174[k]
                   + pb_x[k] * hk_218[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pa_z, pb_x, fl0_5, fl1_5, gl_66, hi0_131, \
                         hi0_133, hi1_175, hi1_177, hk_219, hk_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_3 * hi0_131[k]
                   - f_4 * hi1_175[k]
                   + pb_x[k] * hk_219[k];

        t_257[k] = f_3 * hi0_133[k]
                   - f_4 * hi1_177[k]
                   + pb_x[k] * hk_220[k];

        t_258[k] = f_19 * fl0_5[k]
                   - f_20 * fl1_5[k]
                   + pa_z[k] * gl_66[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pb_y, pb_z, gk_118, gk_140, gk_141, hi0_129, \
                         hi0_130, hi1_173, hi1_174, hk_221, hk_223, \
                         hk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_14 * gk_118[k]
                   + pb_z[k] * hk_221[k];

        t_260[k] = f_15 * gk_140[k]
                   + f_11 * hi0_129[k]
                   - f_12 * hi1_173[k]
                   + pb_y[k] * hk_223[k];

        t_261[k] = f_15 * gk_141[k]
                   + f_9 * hi0_130[k]
                   - f_10 * hi1_174[k]
                   + pb_y[k] * hk_224[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, pb_y, gk_142, gk_143, gk_144, hi0_131, hi0_132, \
                         hi0_133, hi1_175, hi1_176, hi1_177, hk_225, hk_226, \
                         hk_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_15 * gk_142[k]
                   + f_7 * hi0_131[k]
                   - f_8 * hi1_175[k]
                   + pb_y[k] * hk_225[k];

        t_263[k] = f_15 * gk_143[k]
                   + f_5 * hi0_132[k]
                   - f_6 * hi1_176[k]
                   + pb_y[k] * hk_226[k];

        t_264[k] = f_15 * gk_144[k]
                   + f_3 * hi0_133[k]
                   - f_4 * hi1_177[k]
                   + pb_y[k] * hk_227[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pa_y, pb_x, pb_y, fl0_7, fl1_7, gk_145, gl_73, \
                         hi0_134, hi1_178, hk_228, hk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_15 * gk_145[k]
                   + pb_y[k] * hk_228[k];

        t_266[k] = f_21 * fl0_7[k]
                   - f_22 * fl1_7[k]
                   + pa_y[k] * gl_73[k];

        t_267[k] = f_1 * hi0_134[k]
                   - f_2 * hi1_178[k]
                   + pb_x[k] * hk_229[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pb_x, hi0_135, hi0_136, hi0_137, hi1_179, \
                         hi1_180, hi1_181, hk_230, hk_231, hk_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_11 * hi0_135[k]
                   - f_12 * hi1_179[k]
                   + pb_x[k] * hk_230[k];

        t_269[k] = f_11 * hi0_136[k]
                   - f_12 * hi1_180[k]
                   + pb_x[k] * hk_231[k];

        t_270[k] = f_9 * hi0_137[k]
                   - f_10 * hi1_181[k]
                   + pb_x[k] * hk_232[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pb_x, hi0_138, hi0_139, hi0_140, hi1_182, \
                         hi1_183, hi1_184, hk_233, hk_234, hk_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_9 * hi0_138[k]
                   - f_10 * hi1_182[k]
                   + pb_x[k] * hk_233[k];

        t_272[k] = f_7 * hi0_139[k]
                   - f_8 * hi1_183[k]
                   + pb_x[k] * hk_234[k];

        t_273[k] = f_7 * hi0_140[k]
                   - f_8 * hi1_184[k]
                   + pb_x[k] * hk_235[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pb_x, hi0_141, hi0_142, hi0_143, hi1_185, \
                         hi1_186, hi1_187, hk_236, hk_237, hk_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_7 * hi0_141[k]
                   - f_8 * hi1_185[k]
                   + pb_x[k] * hk_236[k];

        t_275[k] = f_5 * hi0_142[k]
                   - f_6 * hi1_186[k]
                   + pb_x[k] * hk_237[k];

        t_276[k] = f_5 * hi0_143[k]
                   - f_6 * hi1_187[k]
                   + pb_x[k] * hk_238[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pb_x, hi0_144, hi0_145, hi0_146, hi1_188, \
                         hi1_189, hi1_190, hk_239, hk_240, hk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_5 * hi0_144[k]
                   - f_6 * hi1_188[k]
                   + pb_x[k] * hk_239[k];

        t_278[k] = f_5 * hi0_145[k]
                   - f_6 * hi1_189[k]
                   + pb_x[k] * hk_240[k];

        t_279[k] = f_3 * hi0_146[k]
                   - f_4 * hi1_190[k]
                   + pb_x[k] * hk_241[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_x, hi0_147, hi0_148, hi0_149, hi1_191, \
                         hi1_192, hi1_193, hk_242, hk_243, hk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_3 * hi0_147[k]
                   - f_4 * hi1_191[k]
                   + pb_x[k] * hk_242[k];

        t_281[k] = f_3 * hi0_148[k]
                   - f_4 * hi1_192[k]
                   + pb_x[k] * hk_243[k];

        t_282[k] = f_3 * hi0_149[k]
                   - f_4 * hi1_193[k]
                   + pb_x[k] * hk_244[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pa_z, pb_x, pb_z, fl0_6, fl1_6, gk_138, gl_67, \
                         hi0_151, hi1_195, hk_245, hk_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_3 * hi0_151[k]
                   - f_4 * hi1_195[k]
                   + pb_x[k] * hk_245[k];

        t_284[k] = f_21 * fl0_6[k]
                   - f_22 * fl1_6[k]
                   + pa_z[k] * gl_67[k];

        t_285[k] = f_15 * gk_138[k]
                   + pb_z[k] * hk_246[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, pb_y, gk_152, gk_153, gk_154, hi0_147, hi0_148, \
                         hi0_149, hi1_191, hi1_192, hi1_193, hk_248, hk_249, \
                         hk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * gk_152[k]
                   + f_11 * hi0_147[k]
                   - f_12 * hi1_191[k]
                   + pb_y[k] * hk_248[k];

        t_287[k] = f_14 * gk_153[k]
                   + f_9 * hi0_148[k]
                   - f_10 * hi1_192[k]
                   + pb_y[k] * hk_249[k];

        t_288[k] = f_14 * gk_154[k]
                   + f_7 * hi0_149[k]
                   - f_8 * hi1_193[k]
                   + pb_y[k] * hk_250[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, pb_y, gk_155, gk_156, gk_157, hi0_150, hi0_151, \
                         hi1_194, hi1_195, hk_251, hk_252, hk_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_14 * gk_155[k]
                   + f_5 * hi0_150[k]
                   - f_6 * hi1_194[k]
                   + pb_y[k] * hk_251[k];

        t_290[k] = f_14 * gk_156[k]
                   + f_3 * hi0_151[k]
                   - f_4 * hi1_195[k]
                   + pb_y[k] * hk_252[k];

        t_291[k] = f_14 * gk_157[k]
                   + pb_y[k] * hk_253[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, fl0_8, fl1_8, gk_159, gk_161, \
                         gk_163, gl_74, gl_76, gl_78, gl_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_19 * fl0_8[k]
                   - f_20 * fl1_8[k]
                   + pa_y[k] * gl_74[k];

        t_293[k] = f_14 * gk_159[k]
                   + pa_y[k] * gl_76[k];

        t_294[k] = f_15 * gk_161[k]
                   + pa_y[k] * gl_78[k];

        t_295[k] = f_16 * gk_163[k]
                   + pa_y[k] * gl_80[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pa_y, gk_164, gk_166, gk_167, \
                         gk_168, gk_170, gl_81, gl_83, gl_84, gl_85, \
                         gl_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_14 * gk_164[k]
                   + pa_y[k] * gl_81[k];

        t_297[k] = f_0 * gk_166[k]
                   + pa_y[k] * gl_83[k];

        t_298[k] = f_15 * gk_167[k]
                   + pa_y[k] * gl_84[k];

        t_299[k] = f_14 * gk_168[k]
                   + pa_y[k] * gl_85[k];

        t_300[k] = f_17 * gk_170[k]
                   + pa_y[k] * gl_87[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_y, gk_171, gk_172, gk_173, gk_176, \
                         gl_88, gl_89, gl_90, gl_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_16 * gk_171[k]
                   + pa_y[k] * gl_88[k];

        t_302[k] = f_15 * gk_172[k]
                   + pa_y[k] * gl_89[k];

        t_303[k] = f_14 * gk_173[k]
                   + pa_y[k] * gl_90[k];

        t_304[k] = f_18 * gk_176[k]
                   + pa_y[k] * gl_92[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pb_z, gk_150, gk_178, gk_179, \
                         gk_180, gl_93, gl_94, gl_95, hk_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_16 * gk_150[k]
                   + pb_z[k] * hk_258[k];

        t_306[k] = f_17 * gk_178[k]
                   + pa_y[k] * gl_93[k];

        t_307[k] = f_0 * gk_179[k]
                   + pa_y[k] * gl_94[k];

        t_308[k] = f_16 * gk_180[k]
                   + pa_y[k] * gl_95[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_y, pb_y, gk_181, gk_182, gk_183, \
                         gl_96, gl_97, gl_98, hk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_15 * gk_181[k]
                   + pa_y[k] * gl_96[k];

        t_310[k] = f_14 * gk_182[k]
                   + pa_y[k] * gl_97[k];

        t_311[k] = f_13 * gk_183[k]
                   + pb_y[k] * hk_265[k];

        t_312[k] = pa_y[k] * gl_98[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pb_x, pb_z, gk_158, hi0_153, hi0_154, \
                         hi0_155, hi1_206, hi1_208, hi1_209, hk_266, hk_268, \
                         hk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_1 * hi0_153[k]
                   - f_2 * hi1_206[k]
                   + pb_x[k] * hk_266[k];

        t_314[k] = f_0 * gk_158[k]
                   + pb_z[k] * hk_266[k];

        t_315[k] = f_11 * hi0_154[k]
                   - f_12 * hi1_208[k]
                   + pb_x[k] * hk_268[k];

        t_316[k] = f_11 * hi0_155[k]
                   - f_12 * hi1_209[k]
                   + pb_x[k] * hk_269[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pb_x, hi0_156, hi0_157, hi0_158, hi1_210, \
                         hi1_211, hi1_212, hk_270, hk_271, hk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_9 * hi0_156[k]
                   - f_10 * hi1_210[k]
                   + pb_x[k] * hk_270[k];

        t_318[k] = f_9 * hi0_157[k]
                   - f_10 * hi1_211[k]
                   + pb_x[k] * hk_271[k];

        t_319[k] = f_7 * hi0_158[k]
                   - f_8 * hi1_212[k]
                   + pb_x[k] * hk_272[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pb_x, hi0_159, hi0_160, hi0_161, hi1_213, \
                         hi1_214, hi1_215, hk_273, hk_274, hk_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_7 * hi0_159[k]
                   - f_8 * hi1_213[k]
                   + pb_x[k] * hk_273[k];

        t_321[k] = f_7 * hi0_160[k]
                   - f_8 * hi1_214[k]
                   + pb_x[k] * hk_274[k];

        t_322[k] = f_5 * hi0_161[k]
                   - f_6 * hi1_215[k]
                   + pb_x[k] * hk_275[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, pb_x, hi0_162, hi0_163, hi0_164, hi1_216, \
                         hi1_217, hi1_218, hk_276, hk_277, hk_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_5 * hi0_162[k]
                   - f_6 * hi1_216[k]
                   + pb_x[k] * hk_276[k];

        t_324[k] = f_5 * hi0_163[k]
                   - f_6 * hi1_217[k]
                   + pb_x[k] * hk_277[k];

        t_325[k] = f_5 * hi0_164[k]
                   - f_6 * hi1_218[k]
                   + pb_x[k] * hk_278[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, pb_x, hi0_165, hi0_166, hi0_167, hi1_219, \
                         hi1_220, hi1_221, hk_279, hk_280, hk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_3 * hi0_165[k]
                   - f_4 * hi1_219[k]
                   + pb_x[k] * hk_279[k];

        t_327[k] = f_3 * hi0_166[k]
                   - f_4 * hi1_220[k]
                   + pb_x[k] * hk_280[k];

        t_328[k] = f_3 * hi0_167[k]
                   - f_4 * hi1_221[k]
                   + pb_x[k] * hk_281[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_x, pb_y, hi0_165, hi0_168, hi0_170, hi1_219, \
                         hi1_222, hi1_224, hk_282, hk_283, hk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_3 * hi0_168[k]
                   - f_4 * hi1_222[k]
                   + pb_x[k] * hk_282[k];

        t_330[k] = f_3 * hi0_170[k]
                   - f_4 * hi1_224[k]
                   + pb_x[k] * hk_283[k];

        t_331[k] = f_1 * hi0_165[k]
                   - f_2 * hi1_219[k]
                   + pb_y[k] * hk_284[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_y, pb_z, gk_176, hi0_166, hi0_167, hi1_220, \
                         hi1_221, hk_284, hk_286, hk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_0 * gk_176[k]
                   + pb_z[k] * hk_284[k];

        t_333[k] = f_11 * hi0_166[k]
                   - f_12 * hi1_220[k]
                   + pb_y[k] * hk_286[k];

        t_334[k] = f_9 * hi0_167[k]
                   - f_10 * hi1_221[k]
                   + pb_y[k] * hk_287[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_y, hi0_168, hi0_169, hi0_170, hi1_222, \
                         hi1_223, hi1_224, hk_288, hk_289, hk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_7 * hi0_168[k]
                   - f_8 * hi1_222[k]
                   + pb_y[k] * hk_288[k];

        t_336[k] = f_5 * hi0_169[k]
                   - f_6 * hi1_223[k]
                   + pb_y[k] * hk_289[k];

        t_337[k] = f_3 * hi0_170[k]
                   - f_4 * hi1_224[k]
                   + pb_y[k] * hk_290[k];
    }

#pragma omp simd aligned(t_338, pb_z, gk_183, hi0_170, hi1_224, \
                         hk_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_0 * gk_183[k]
                   + f_1 * hi0_170[k]
                   - f_2 * hi1_224[k]
                   + pb_z[k] * hk_291[k];
    }
}

auto
compute_prim_hl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
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
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);
    const auto f_18 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_44 = buffer.data(gk + 44);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_6 = buffer.data(hi0 + 6);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_10 = buffer.data(hi0 + 10);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_18 = buffer.data(hi0 + 18);
    const auto *hi0_19 = buffer.data(hi0 + 19);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_22 = buffer.data(hi0 + 22);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_25 = buffer.data(hi0 + 25);
    const auto *hi0_26 = buffer.data(hi0 + 26);
    const auto *hi0_32 = buffer.data(hi0 + 32);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_36 = buffer.data(hi0 + 36);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_38 = buffer.data(hi0 + 38);
    const auto *hi0_39 = buffer.data(hi0 + 39);
    const auto *hi0_41 = buffer.data(hi0 + 41);
    const auto *hi0_42 = buffer.data(hi0 + 42);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_44 = buffer.data(hi0 + 44);
    const auto *hi0_45 = buffer.data(hi0 + 45);
    const auto *hi0_47 = buffer.data(hi0 + 47);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_4 = buffer.data(hi1 + 4);
    const auto *hi1_5 = buffer.data(hi1 + 5);
    const auto *hi1_6 = buffer.data(hi1 + 6);
    const auto *hi1_7 = buffer.data(hi1 + 7);
    const auto *hi1_8 = buffer.data(hi1 + 8);
    const auto *hi1_10 = buffer.data(hi1 + 10);
    const auto *hi1_11 = buffer.data(hi1 + 11);
    const auto *hi1_12 = buffer.data(hi1 + 12);
    const auto *hi1_13 = buffer.data(hi1 + 13);
    const auto *hi1_14 = buffer.data(hi1 + 14);
    const auto *hi1_16 = buffer.data(hi1 + 16);
    const auto *hi1_17 = buffer.data(hi1 + 17);
    const auto *hi1_18 = buffer.data(hi1 + 18);
    const auto *hi1_19 = buffer.data(hi1 + 19);
    const auto *hi1_20 = buffer.data(hi1 + 20);
    const auto *hi1_22 = buffer.data(hi1 + 22);
    const auto *hi1_23 = buffer.data(hi1 + 23);
    const auto *hi1_24 = buffer.data(hi1 + 24);
    const auto *hi1_25 = buffer.data(hi1 + 25);
    const auto *hi1_26 = buffer.data(hi1 + 26);
    const auto *hi1_32 = buffer.data(hi1 + 32);
    const auto *hi1_35 = buffer.data(hi1 + 35);
    const auto *hi1_36 = buffer.data(hi1 + 36);
    const auto *hi1_37 = buffer.data(hi1 + 37);
    const auto *hi1_38 = buffer.data(hi1 + 38);
    const auto *hi1_39 = buffer.data(hi1 + 39);
    const auto *hi1_41 = buffer.data(hi1 + 41);
    const auto *hi1_42 = buffer.data(hi1 + 42);
    const auto *hi1_43 = buffer.data(hi1 + 43);
    const auto *hi1_44 = buffer.data(hi1 + 44);
    const auto *hi1_45 = buffer.data(hi1 + 45);
    const auto *hi1_47 = buffer.data(hi1 + 47);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_65 = buffer.data(hk + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, fl0_0, fl1_0, gk_0, gl_0, gl_1, \
                         hi0_0, hi1_0, hk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = pa_y[k] * gl_0[k];

        t_2[k] = pa_z[k] * gl_0[k];

        t_3[k] = f_3 * fl0_0[k]
                 - f_4 * fl1_0[k]
                 + pa_y[k] * gl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, gk_4, gk_5, gk_6, hi0_4, hi0_5, hi0_6, hi1_4, \
                         hi1_5, hi1_6, hk_4, hk_5, hk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gk_4[k]
                 + f_6 * hi0_4[k]
                 - f_7 * hi1_4[k]
                 + pb_x[k] * hk_4[k];

        t_5[k] = f_5 * gk_5[k]
                 + f_8 * hi0_5[k]
                 - f_9 * hi1_5[k]
                 + pb_x[k] * hk_5[k];

        t_6[k] = f_5 * gk_6[k]
                 + f_10 * hi0_6[k]
                 - f_11 * hi1_6[k]
                 + pb_x[k] * hk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, fl0_3, fl1_3, gk_7, gk_8, gl_9, hi0_7, \
                         hi0_8, hi1_7, hi1_8, hk_7, hk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * gk_7[k]
                 + f_12 * hi0_7[k]
                 - f_13 * hi1_7[k]
                 + pb_x[k] * hk_7[k];

        t_8[k] = f_5 * gk_8[k]
                 + f_14 * hi0_8[k]
                 - f_15 * hi1_8[k]
                 + pb_x[k] * hk_8[k];

        t_9[k] = f_16 * fl0_3[k]
                 - f_17 * fl1_3[k]
                 + pa_x[k] * gl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, fl0_0, fl1_0, gk_11, gk_12, gl_2, \
                         hi0_10, hi0_11, hi1_10, hi1_11, hk_11, hk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fl0_0[k]
                  - f_4 * fl1_0[k]
                  + pa_z[k] * gl_2[k];

        t_11[k] = f_5 * gk_11[k]
                  + f_6 * hi0_10[k]
                  - f_7 * hi1_10[k]
                  + pb_x[k] * hk_11[k];

        t_12[k] = f_5 * gk_12[k]
                  + f_8 * hi0_11[k]
                  - f_9 * hi1_11[k]
                  + pb_x[k] * hk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, gk_13, gk_14, gk_15, hi0_12, hi0_13, hi0_14, \
                         hi1_12, hi1_13, hi1_14, hk_13, hk_14, hk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gk_13[k]
                  + f_10 * hi0_12[k]
                  - f_11 * hi1_12[k]
                  + pb_x[k] * hk_13[k];

        t_14[k] = f_5 * gk_14[k]
                  + f_12 * hi0_13[k]
                  - f_13 * hi1_13[k]
                  + pb_x[k] * hk_14[k];

        t_15[k] = f_5 * gk_15[k]
                  + f_14 * hi0_14[k]
                  - f_15 * hi1_14[k]
                  + pb_x[k] * hk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_x, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gk_17, gl_3, gl_16, hi0_16, hi1_16, hk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_16 * fl0_4[k]
                  - f_17 * fl1_4[k]
                  + pa_x[k] * gl_16[k];

        t_17[k] = f_16 * fl0_1[k]
                  - f_17 * fl1_1[k]
                  + pa_y[k] * gl_3[k];

        t_18[k] = f_18 * gk_17[k]
                  + f_6 * hi0_16[k]
                  - f_7 * hi1_16[k]
                  + pb_x[k] * hk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gk_18, gk_19, gk_20, hi0_17, hi0_18, hi0_19, \
                         hi1_17, hi1_18, hi1_19, hk_19, hk_20, hk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_18 * gk_18[k]
                  + f_8 * hi0_17[k]
                  - f_9 * hi1_17[k]
                  + pb_x[k] * hk_19[k];

        t_20[k] = f_18 * gk_19[k]
                  + f_10 * hi0_18[k]
                  - f_11 * hi1_18[k]
                  + pb_x[k] * hk_20[k];

        t_21[k] = f_18 * gk_20[k]
                  + f_12 * hi0_19[k]
                  - f_13 * hi1_19[k]
                  + pb_x[k] * hk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, pb_x, fl0_5, fl1_5, gk_21, gl_4, \
                         gl_5, gl_17, hi0_20, hi1_20, hk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_18 * gk_21[k]
                  + f_14 * hi0_20[k]
                  - f_15 * hi1_20[k]
                  + pb_x[k] * hk_22[k];

        t_23[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_x[k] * gl_17[k];

        t_24[k] = pa_z[k] * gl_4[k];

        t_25[k] = pa_z[k] * gl_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, gl_6, gl_7, \
                         gl_8, gl_10, gl_11, gl_12, gl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * gl_6[k];

        t_27[k] = pa_z[k] * gl_7[k];

        t_28[k] = pa_z[k] * gl_8[k];

        t_29[k] = pa_y[k] * gl_10[k];

        t_30[k] = pa_y[k] * gl_11[k];

        t_31[k] = pa_y[k] * gl_12[k];

        t_32[k] = pa_y[k] * gl_13[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pa_z, pb_x, fl0_2, fl1_2, gk_23, gl_10, \
                         gl_14, gl_15, hi0_22, hi1_22, hk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * gl_14[k];

        t_34[k] = pa_y[k] * gl_15[k];

        t_35[k] = f_16 * fl0_2[k]
                  - f_17 * fl1_2[k]
                  + pa_z[k] * gl_10[k];

        t_36[k] = f_18 * gk_23[k]
                  + f_6 * hi0_22[k]
                  - f_7 * hi1_22[k]
                  + pb_x[k] * hk_34[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, gk_24, gk_25, gk_26, hi0_23, hi0_24, hi0_25, \
                         hi1_23, hi1_24, hi1_25, hk_35, hk_36, hk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_18 * gk_24[k]
                  + f_8 * hi0_23[k]
                  - f_9 * hi1_23[k]
                  + pb_x[k] * hk_35[k];

        t_38[k] = f_18 * gk_25[k]
                  + f_10 * hi0_24[k]
                  - f_11 * hi1_24[k]
                  + pb_x[k] * hk_36[k];

        t_39[k] = f_18 * gk_26[k]
                  + f_12 * hi0_25[k]
                  - f_13 * hi1_25[k]
                  + pb_x[k] * hk_37[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, fl0_8, fl1_8, gk_27, gl_18, \
                         gl_19, gl_21, hi0_26, hi1_26, hk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_18 * gk_27[k]
                  + f_14 * hi0_26[k]
                  - f_15 * hi1_26[k]
                  + pb_x[k] * hk_38[k];

        t_41[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_x[k] * gl_18[k];

        t_42[k] = pa_x[k] * gl_19[k];

        t_43[k] = pa_x[k] * gl_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, pa_x, gl_22, gl_23, gl_24, \
                         gl_25, gl_26, gl_27, gl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gl_22[k];

        t_45[k] = pa_x[k] * gl_23[k];

        t_46[k] = pa_x[k] * gl_24[k];

        t_47[k] = pa_x[k] * gl_25[k];

        t_48[k] = pa_x[k] * gl_26[k];

        t_49[k] = pa_x[k] * gl_27[k];

        t_50[k] = pa_x[k] * gl_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, fl0_5, fl1_5, gk_29, gl_19, gl_20, \
                         hi0_32, hi1_32, hk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * gk_29[k]
                  + f_1 * hi0_32[k]
                  - f_2 * hi1_32[k]
                  + pb_y[k] * hk_48[k];

        t_52[k] = pa_z[k] * gl_19[k];

        t_53[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_z[k] * gl_20[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, gk_32, gk_33, gk_34, hi0_35, hi0_36, hi0_37, \
                         hi1_35, hi1_36, hi1_37, hk_51, hk_52, hk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * gk_32[k]
                  + f_6 * hi0_35[k]
                  - f_7 * hi1_35[k]
                  + pb_y[k] * hk_51[k];

        t_55[k] = f_5 * gk_33[k]
                  + f_8 * hi0_36[k]
                  - f_9 * hi1_36[k]
                  + pb_y[k] * hk_52[k];

        t_56[k] = f_5 * gk_34[k]
                  + f_10 * hi0_37[k]
                  - f_11 * hi1_37[k]
                  + pb_y[k] * hk_53[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pb_y, fl0_7, fl1_7, gk_35, gk_36, gl_27, \
                         hi0_38, hi0_39, hi1_38, hi1_39, hk_54, hk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * gk_35[k]
                  + f_12 * hi0_38[k]
                  - f_13 * hi1_38[k]
                  + pb_y[k] * hk_54[k];

        t_58[k] = f_5 * gk_36[k]
                  + f_14 * hi0_39[k]
                  - f_15 * hi1_39[k]
                  + pb_y[k] * hk_55[k];

        t_59[k] = f_16 * fl0_7[k]
                  - f_17 * fl1_7[k]
                  + pa_y[k] * gl_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_z, pb_y, fl0_6, fl1_6, gk_38, gk_39, gl_21, \
                         hi0_41, hi0_42, hi1_41, hi1_42, hk_58, hk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_16 * fl0_6[k]
                  - f_17 * fl1_6[k]
                  + pa_z[k] * gl_21[k];

        t_61[k] = f_18 * gk_38[k]
                  + f_6 * hi0_41[k]
                  - f_7 * hi1_41[k]
                  + pb_y[k] * hk_58[k];

        t_62[k] = f_18 * gk_39[k]
                  + f_8 * hi0_42[k]
                  - f_9 * hi1_42[k]
                  + pb_y[k] * hk_59[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, gk_40, gk_41, gk_42, hi0_43, hi0_44, hi0_45, \
                         hi1_43, hi1_44, hi1_45, hk_60, hk_61, hk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_18 * gk_40[k]
                  + f_10 * hi0_43[k]
                  - f_11 * hi1_43[k]
                  + pb_y[k] * hk_60[k];

        t_64[k] = f_18 * gk_41[k]
                  + f_12 * hi0_44[k]
                  - f_13 * hi1_44[k]
                  + pb_y[k] * hk_61[k];

        t_65[k] = f_18 * gk_42[k]
                  + f_14 * hi0_45[k]
                  - f_15 * hi1_45[k]
                  + pb_y[k] * hk_62[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pb_z, fl0_8, fl1_8, gk_44, gl_28, gl_29, \
                         hi0_47, hi1_47, hk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_y[k] * gl_28[k];

        t_67[k] = pa_y[k] * gl_29[k];

        t_68[k] = f_0 * gk_44[k]
                  + f_1 * hi0_47[k]
                  - f_2 * hi1_47[k]
                  + pb_z[k] * hk_65[k];
    }
}

auto
compute_prim_hl_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
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
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);
    const auto f_18 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_44 = buffer.data(gk + 44);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_6 = buffer.data(hi0 + 6);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_10 = buffer.data(hi0 + 10);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_18 = buffer.data(hi0 + 18);
    const auto *hi0_19 = buffer.data(hi0 + 19);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_22 = buffer.data(hi0 + 22);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_25 = buffer.data(hi0 + 25);
    const auto *hi0_26 = buffer.data(hi0 + 26);
    const auto *hi0_32 = buffer.data(hi0 + 32);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_36 = buffer.data(hi0 + 36);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_38 = buffer.data(hi0 + 38);
    const auto *hi0_39 = buffer.data(hi0 + 39);
    const auto *hi0_41 = buffer.data(hi0 + 41);
    const auto *hi0_42 = buffer.data(hi0 + 42);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_44 = buffer.data(hi0 + 44);
    const auto *hi0_45 = buffer.data(hi0 + 45);
    const auto *hi0_47 = buffer.data(hi0 + 47);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_33 = buffer.data(hi1 + 33);
    const auto *hi1_35 = buffer.data(hi1 + 35);
    const auto *hi1_37 = buffer.data(hi1 + 37);
    const auto *hi1_39 = buffer.data(hi1 + 39);
    const auto *hi1_40 = buffer.data(hi1 + 40);
    const auto *hi1_51 = buffer.data(hi1 + 51);
    const auto *hi1_53 = buffer.data(hi1 + 53);
    const auto *hi1_55 = buffer.data(hi1 + 55);
    const auto *hi1_56 = buffer.data(hi1 + 56);
    const auto *hi1_62 = buffer.data(hi1 + 62);
    const auto *hi1_65 = buffer.data(hi1 + 65);
    const auto *hi1_67 = buffer.data(hi1 + 67);
    const auto *hi1_69 = buffer.data(hi1 + 69);
    const auto *hi1_71 = buffer.data(hi1 + 71);
    const auto *hi1_72 = buffer.data(hi1 + 72);
    const auto *hi1_87 = buffer.data(hi1 + 87);
    const auto *hi1_89 = buffer.data(hi1 + 89);
    const auto *hi1_91 = buffer.data(hi1 + 91);
    const auto *hi1_92 = buffer.data(hi1 + 92);
    const auto *hi1_98 = buffer.data(hi1 + 98);
    const auto *hi1_136 = buffer.data(hi1 + 136);
    const auto *hi1_162 = buffer.data(hi1 + 162);
    const auto *hi1_163 = buffer.data(hi1 + 163);
    const auto *hi1_164 = buffer.data(hi1 + 164);
    const auto *hi1_165 = buffer.data(hi1 + 165);
    const auto *hi1_166 = buffer.data(hi1 + 166);
    const auto *hi1_180 = buffer.data(hi1 + 180);
    const auto *hi1_181 = buffer.data(hi1 + 181);
    const auto *hi1_182 = buffer.data(hi1 + 182);
    const auto *hi1_183 = buffer.data(hi1 + 183);
    const auto *hi1_184 = buffer.data(hi1 + 184);
    const auto *hi1_215 = buffer.data(hi1 + 215);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_251 = buffer.data(hk + 251);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, fl0_0, fl1_0, gk_0, gl_0, gl_1, \
                         hi0_0, hi1_0, hk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = pa_y[k] * gl_0[k];

        t_2[k] = pa_z[k] * gl_0[k];

        t_3[k] = f_3 * fl0_0[k]
                 - f_4 * fl1_0[k]
                 + pa_y[k] * gl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, gk_4, gk_5, gk_6, hi0_4, hi0_5, hi0_6, hi1_33, \
                         hi1_35, hi1_37, hk_42, hk_44, hk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gk_4[k]
                 + f_6 * hi0_4[k]
                 - f_7 * hi1_33[k]
                 + pb_x[k] * hk_42[k];

        t_5[k] = f_5 * gk_5[k]
                 + f_8 * hi0_5[k]
                 - f_9 * hi1_35[k]
                 + pb_x[k] * hk_44[k];

        t_6[k] = f_5 * gk_6[k]
                 + f_10 * hi0_6[k]
                 - f_11 * hi1_37[k]
                 + pb_x[k] * hk_46[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, fl0_3, fl1_3, gk_7, gk_8, gl_9, hi0_7, \
                         hi0_8, hi1_39, hi1_40, hk_48, hk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * gk_7[k]
                 + f_12 * hi0_7[k]
                 - f_13 * hi1_39[k]
                 + pb_x[k] * hk_48[k];

        t_8[k] = f_5 * gk_8[k]
                 + f_14 * hi0_8[k]
                 - f_15 * hi1_40[k]
                 + pb_x[k] * hk_50[k];

        t_9[k] = f_16 * fl0_3[k]
                 - f_17 * fl1_3[k]
                 + pa_x[k] * gl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, fl0_0, fl1_0, gk_11, gk_12, gl_2, \
                         hi0_10, hi0_11, hi1_51, hi1_53, hk_60, hk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fl0_0[k]
                  - f_4 * fl1_0[k]
                  + pa_z[k] * gl_2[k];

        t_11[k] = f_5 * gk_11[k]
                  + f_6 * hi0_10[k]
                  - f_7 * hi1_51[k]
                  + pb_x[k] * hk_60[k];

        t_12[k] = f_5 * gk_12[k]
                  + f_8 * hi0_11[k]
                  - f_9 * hi1_53[k]
                  + pb_x[k] * hk_62[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, gk_13, gk_14, gk_15, hi0_12, hi0_13, hi0_14, \
                         hi1_55, hi1_56, hi1_62, hk_64, hk_66, hk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gk_13[k]
                  + f_10 * hi0_12[k]
                  - f_11 * hi1_55[k]
                  + pb_x[k] * hk_64[k];

        t_14[k] = f_5 * gk_14[k]
                  + f_12 * hi0_13[k]
                  - f_13 * hi1_56[k]
                  + pb_x[k] * hk_66[k];

        t_15[k] = f_5 * gk_15[k]
                  + f_14 * hi0_14[k]
                  - f_15 * hi1_62[k]
                  + pb_x[k] * hk_67[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_x, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gk_17, gl_3, gl_16, hi0_16, hi1_65, hk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_16 * fl0_4[k]
                  - f_17 * fl1_4[k]
                  + pa_x[k] * gl_16[k];

        t_17[k] = f_16 * fl0_1[k]
                  - f_17 * fl1_1[k]
                  + pa_y[k] * gl_3[k];

        t_18[k] = f_18 * gk_17[k]
                  + f_6 * hi0_16[k]
                  - f_7 * hi1_65[k]
                  + pb_x[k] * hk_75[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gk_18, gk_19, gk_20, hi0_17, hi0_18, hi0_19, \
                         hi1_67, hi1_69, hi1_71, hk_77, hk_79, hk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_18 * gk_18[k]
                  + f_8 * hi0_17[k]
                  - f_9 * hi1_67[k]
                  + pb_x[k] * hk_77[k];

        t_20[k] = f_18 * gk_19[k]
                  + f_10 * hi0_18[k]
                  - f_11 * hi1_69[k]
                  + pb_x[k] * hk_79[k];

        t_21[k] = f_18 * gk_20[k]
                  + f_12 * hi0_19[k]
                  - f_13 * hi1_71[k]
                  + pb_x[k] * hk_81[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, pb_x, fl0_5, fl1_5, gk_21, gl_4, \
                         gl_5, gl_17, hi0_20, hi1_72, hk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_18 * gk_21[k]
                  + f_14 * hi0_20[k]
                  - f_15 * hi1_72[k]
                  + pb_x[k] * hk_83[k];

        t_23[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_x[k] * gl_17[k];

        t_24[k] = pa_z[k] * gl_4[k];

        t_25[k] = pa_z[k] * gl_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, gl_6, gl_7, \
                         gl_8, gl_10, gl_11, gl_12, gl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * gl_6[k];

        t_27[k] = pa_z[k] * gl_7[k];

        t_28[k] = pa_z[k] * gl_8[k];

        t_29[k] = pa_y[k] * gl_10[k];

        t_30[k] = pa_y[k] * gl_11[k];

        t_31[k] = pa_y[k] * gl_12[k];

        t_32[k] = pa_y[k] * gl_13[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pa_z, pb_x, fl0_2, fl1_2, gk_23, gl_10, \
                         gl_14, gl_15, hi0_22, hi1_87, hk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * gl_14[k];

        t_34[k] = pa_y[k] * gl_15[k];

        t_35[k] = f_16 * fl0_2[k]
                  - f_17 * fl1_2[k]
                  + pa_z[k] * gl_10[k];

        t_36[k] = f_18 * gk_23[k]
                  + f_6 * hi0_22[k]
                  - f_7 * hi1_87[k]
                  + pb_x[k] * hk_102[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, gk_24, gk_25, gk_26, hi0_23, hi0_24, hi0_25, \
                         hi1_89, hi1_91, hi1_92, hk_104, hk_106, \
                         hk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_18 * gk_24[k]
                  + f_8 * hi0_23[k]
                  - f_9 * hi1_89[k]
                  + pb_x[k] * hk_104[k];

        t_38[k] = f_18 * gk_25[k]
                  + f_10 * hi0_24[k]
                  - f_11 * hi1_91[k]
                  + pb_x[k] * hk_106[k];

        t_39[k] = f_18 * gk_26[k]
                  + f_12 * hi0_25[k]
                  - f_13 * hi1_92[k]
                  + pb_x[k] * hk_108[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, fl0_8, fl1_8, gk_27, gl_18, \
                         gl_19, gl_21, hi0_26, hi1_98, hk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_18 * gk_27[k]
                  + f_14 * hi0_26[k]
                  - f_15 * hi1_98[k]
                  + pb_x[k] * hk_109[k];

        t_41[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_x[k] * gl_18[k];

        t_42[k] = pa_x[k] * gl_19[k];

        t_43[k] = pa_x[k] * gl_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, pa_x, gl_22, gl_23, gl_24, \
                         gl_25, gl_26, gl_27, gl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gl_22[k];

        t_45[k] = pa_x[k] * gl_23[k];

        t_46[k] = pa_x[k] * gl_24[k];

        t_47[k] = pa_x[k] * gl_25[k];

        t_48[k] = pa_x[k] * gl_26[k];

        t_49[k] = pa_x[k] * gl_27[k];

        t_50[k] = pa_x[k] * gl_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, fl0_5, fl1_5, gk_29, gl_19, gl_20, \
                         hi0_32, hi1_136, hk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * gk_29[k]
                  + f_1 * hi0_32[k]
                  - f_2 * hi1_136[k]
                  + pb_y[k] * hk_154[k];

        t_52[k] = pa_z[k] * gl_19[k];

        t_53[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_z[k] * gl_20[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, gk_32, gk_33, gk_34, hi0_35, hi0_36, hi0_37, \
                         hi1_162, hi1_163, hi1_164, hk_188, hk_189, \
                         hk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * gk_32[k]
                  + f_6 * hi0_35[k]
                  - f_7 * hi1_162[k]
                  + pb_y[k] * hk_188[k];

        t_55[k] = f_5 * gk_33[k]
                  + f_8 * hi0_36[k]
                  - f_9 * hi1_163[k]
                  + pb_y[k] * hk_189[k];

        t_56[k] = f_5 * gk_34[k]
                  + f_10 * hi0_37[k]
                  - f_11 * hi1_164[k]
                  + pb_y[k] * hk_190[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pb_y, fl0_7, fl1_7, gk_35, gk_36, gl_27, \
                         hi0_38, hi0_39, hi1_165, hi1_166, hk_191, \
                         hk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * gk_35[k]
                  + f_12 * hi0_38[k]
                  - f_13 * hi1_165[k]
                  + pb_y[k] * hk_191[k];

        t_58[k] = f_5 * gk_36[k]
                  + f_14 * hi0_39[k]
                  - f_15 * hi1_166[k]
                  + pb_y[k] * hk_192[k];

        t_59[k] = f_16 * fl0_7[k]
                  - f_17 * fl1_7[k]
                  + pa_y[k] * gl_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_z, pb_y, fl0_6, fl1_6, gk_38, gk_39, gl_21, \
                         hi0_41, hi0_42, hi1_180, hi1_181, hk_208, \
                         hk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_16 * fl0_6[k]
                  - f_17 * fl1_6[k]
                  + pa_z[k] * gl_21[k];

        t_61[k] = f_18 * gk_38[k]
                  + f_6 * hi0_41[k]
                  - f_7 * hi1_180[k]
                  + pb_y[k] * hk_208[k];

        t_62[k] = f_18 * gk_39[k]
                  + f_8 * hi0_42[k]
                  - f_9 * hi1_181[k]
                  + pb_y[k] * hk_209[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, gk_40, gk_41, gk_42, hi0_43, hi0_44, hi0_45, \
                         hi1_182, hi1_183, hi1_184, hk_210, hk_211, \
                         hk_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_18 * gk_40[k]
                  + f_10 * hi0_43[k]
                  - f_11 * hi1_182[k]
                  + pb_y[k] * hk_210[k];

        t_64[k] = f_18 * gk_41[k]
                  + f_12 * hi0_44[k]
                  - f_13 * hi1_183[k]
                  + pb_y[k] * hk_211[k];

        t_65[k] = f_18 * gk_42[k]
                  + f_14 * hi0_45[k]
                  - f_15 * hi1_184[k]
                  + pb_y[k] * hk_212[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pb_z, fl0_8, fl1_8, gk_44, gl_28, gl_29, \
                         hi0_47, hi1_215, hk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_y[k] * gl_28[k];

        t_67[k] = pa_y[k] * gl_29[k];

        t_68[k] = f_0 * gk_44[k]
                  + f_1 * hi0_47[k]
                  - f_2 * hi1_215[k]
                  + pb_z[k] * hk_251[k];
    }
}

auto
compute_prim_hl_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
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
    const auto f_15 = 1.5 / p;
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);
    const auto f_18 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_86 = buffer.data(gk + 86);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_1 = buffer.data(hi0 + 1);
    const auto *hi0_2 = buffer.data(hi0 + 2);
    const auto *hi0_3 = buffer.data(hi0 + 3);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_9 = buffer.data(hi0 + 9);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_18 = buffer.data(hi0 + 18);
    const auto *hi0_19 = buffer.data(hi0 + 19);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_33 = buffer.data(hi0 + 33);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_39 = buffer.data(hi0 + 39);
    const auto *hi0_40 = buffer.data(hi0 + 40);
    const auto *hi0_51 = buffer.data(hi0 + 51);
    const auto *hi0_53 = buffer.data(hi0 + 53);
    const auto *hi0_55 = buffer.data(hi0 + 55);
    const auto *hi0_56 = buffer.data(hi0 + 56);
    const auto *hi0_62 = buffer.data(hi0 + 62);
    const auto *hi0_65 = buffer.data(hi0 + 65);
    const auto *hi0_67 = buffer.data(hi0 + 67);
    const auto *hi0_69 = buffer.data(hi0 + 69);
    const auto *hi0_71 = buffer.data(hi0 + 71);
    const auto *hi0_72 = buffer.data(hi0 + 72);
    const auto *hi0_87 = buffer.data(hi0 + 87);
    const auto *hi0_89 = buffer.data(hi0 + 89);
    const auto *hi0_91 = buffer.data(hi0 + 91);
    const auto *hi0_92 = buffer.data(hi0 + 92);
    const auto *hi0_98 = buffer.data(hi0 + 98);
    const auto *hi0_121 = buffer.data(hi0 + 121);
    const auto *hi0_123 = buffer.data(hi0 + 123);
    const auto *hi0_124 = buffer.data(hi0 + 124);
    const auto *hi0_125 = buffer.data(hi0 + 125);
    const auto *hi0_127 = buffer.data(hi0 + 127);
    const auto *hi0_128 = buffer.data(hi0 + 128);
    const auto *hi0_130 = buffer.data(hi0 + 130);
    const auto *hi0_131 = buffer.data(hi0 + 131);
    const auto *hi0_132 = buffer.data(hi0 + 132);
    const auto *hi0_133 = buffer.data(hi0 + 133);
    const auto *hi0_134 = buffer.data(hi0 + 134);
    const auto *hi0_135 = buffer.data(hi0 + 135);
    const auto *hi0_136 = buffer.data(hi0 + 136);
    const auto *hi0_137 = buffer.data(hi0 + 137);
    const auto *hi0_138 = buffer.data(hi0 + 138);
    const auto *hi0_139 = buffer.data(hi0 + 139);
    const auto *hi0_140 = buffer.data(hi0 + 140);
    const auto *hi0_141 = buffer.data(hi0 + 141);
    const auto *hi0_162 = buffer.data(hi0 + 162);
    const auto *hi0_163 = buffer.data(hi0 + 163);
    const auto *hi0_164 = buffer.data(hi0 + 164);
    const auto *hi0_165 = buffer.data(hi0 + 165);
    const auto *hi0_166 = buffer.data(hi0 + 166);
    const auto *hi0_180 = buffer.data(hi0 + 180);
    const auto *hi0_181 = buffer.data(hi0 + 181);
    const auto *hi0_182 = buffer.data(hi0 + 182);
    const auto *hi0_183 = buffer.data(hi0 + 183);
    const auto *hi0_184 = buffer.data(hi0 + 184);
    const auto *hi0_195 = buffer.data(hi0 + 195);
    const auto *hi0_197 = buffer.data(hi0 + 197);
    const auto *hi0_198 = buffer.data(hi0 + 198);
    const auto *hi0_199 = buffer.data(hi0 + 199);
    const auto *hi0_201 = buffer.data(hi0 + 201);
    const auto *hi0_202 = buffer.data(hi0 + 202);
    const auto *hi0_203 = buffer.data(hi0 + 203);
    const auto *hi0_205 = buffer.data(hi0 + 205);
    const auto *hi0_206 = buffer.data(hi0 + 206);
    const auto *hi0_207 = buffer.data(hi0 + 207);
    const auto *hi0_208 = buffer.data(hi0 + 208);
    const auto *hi0_209 = buffer.data(hi0 + 209);
    const auto *hi0_210 = buffer.data(hi0 + 210);
    const auto *hi0_211 = buffer.data(hi0 + 211);
    const auto *hi0_212 = buffer.data(hi0 + 212);
    const auto *hi0_213 = buffer.data(hi0 + 213);
    const auto *hi0_214 = buffer.data(hi0 + 214);
    const auto *hi0_215 = buffer.data(hi0 + 215);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_1 = buffer.data(hi1 + 1);
    const auto *hi1_2 = buffer.data(hi1 + 2);
    const auto *hi1_3 = buffer.data(hi1 + 3);
    const auto *hi1_4 = buffer.data(hi1 + 4);
    const auto *hi1_5 = buffer.data(hi1 + 5);
    const auto *hi1_6 = buffer.data(hi1 + 6);
    const auto *hi1_7 = buffer.data(hi1 + 7);
    const auto *hi1_8 = buffer.data(hi1 + 8);
    const auto *hi1_9 = buffer.data(hi1 + 9);
    const auto *hi1_10 = buffer.data(hi1 + 10);
    const auto *hi1_11 = buffer.data(hi1 + 11);
    const auto *hi1_12 = buffer.data(hi1 + 12);
    const auto *hi1_14 = buffer.data(hi1 + 14);
    const auto *hi1_15 = buffer.data(hi1 + 15);
    const auto *hi1_16 = buffer.data(hi1 + 16);
    const auto *hi1_17 = buffer.data(hi1 + 17);
    const auto *hi1_18 = buffer.data(hi1 + 18);
    const auto *hi1_25 = buffer.data(hi1 + 25);
    const auto *hi1_27 = buffer.data(hi1 + 27);
    const auto *hi1_29 = buffer.data(hi1 + 29);
    const auto *hi1_31 = buffer.data(hi1 + 31);
    const auto *hi1_32 = buffer.data(hi1 + 32);
    const auto *hi1_41 = buffer.data(hi1 + 41);
    const auto *hi1_43 = buffer.data(hi1 + 43);
    const auto *hi1_45 = buffer.data(hi1 + 45);
    const auto *hi1_46 = buffer.data(hi1 + 46);
    const auto *hi1_52 = buffer.data(hi1 + 52);
    const auto *hi1_55 = buffer.data(hi1 + 55);
    const auto *hi1_57 = buffer.data(hi1 + 57);
    const auto *hi1_59 = buffer.data(hi1 + 59);
    const auto *hi1_61 = buffer.data(hi1 + 61);
    const auto *hi1_62 = buffer.data(hi1 + 62);
    const auto *hi1_72 = buffer.data(hi1 + 72);
    const auto *hi1_74 = buffer.data(hi1 + 74);
    const auto *hi1_76 = buffer.data(hi1 + 76);
    const auto *hi1_77 = buffer.data(hi1 + 77);
    const auto *hi1_83 = buffer.data(hi1 + 83);
    const auto *hi1_101 = buffer.data(hi1 + 101);
    const auto *hi1_103 = buffer.data(hi1 + 103);
    const auto *hi1_104 = buffer.data(hi1 + 104);
    const auto *hi1_105 = buffer.data(hi1 + 105);
    const auto *hi1_106 = buffer.data(hi1 + 106);
    const auto *hi1_107 = buffer.data(hi1 + 107);
    const auto *hi1_108 = buffer.data(hi1 + 108);
    const auto *hi1_109 = buffer.data(hi1 + 109);
    const auto *hi1_110 = buffer.data(hi1 + 110);
    const auto *hi1_111 = buffer.data(hi1 + 111);
    const auto *hi1_112 = buffer.data(hi1 + 112);
    const auto *hi1_113 = buffer.data(hi1 + 113);
    const auto *hi1_114 = buffer.data(hi1 + 114);
    const auto *hi1_115 = buffer.data(hi1 + 115);
    const auto *hi1_116 = buffer.data(hi1 + 116);
    const auto *hi1_117 = buffer.data(hi1 + 117);
    const auto *hi1_118 = buffer.data(hi1 + 118);
    const auto *hi1_119 = buffer.data(hi1 + 119);
    const auto *hi1_135 = buffer.data(hi1 + 135);
    const auto *hi1_136 = buffer.data(hi1 + 136);
    const auto *hi1_137 = buffer.data(hi1 + 137);
    const auto *hi1_138 = buffer.data(hi1 + 138);
    const auto *hi1_139 = buffer.data(hi1 + 139);
    const auto *hi1_153 = buffer.data(hi1 + 153);
    const auto *hi1_154 = buffer.data(hi1 + 154);
    const auto *hi1_155 = buffer.data(hi1 + 155);
    const auto *hi1_156 = buffer.data(hi1 + 156);
    const auto *hi1_157 = buffer.data(hi1 + 157);
    const auto *hi1_164 = buffer.data(hi1 + 164);
    const auto *hi1_166 = buffer.data(hi1 + 166);
    const auto *hi1_167 = buffer.data(hi1 + 167);
    const auto *hi1_168 = buffer.data(hi1 + 168);
    const auto *hi1_169 = buffer.data(hi1 + 169);
    const auto *hi1_170 = buffer.data(hi1 + 170);
    const auto *hi1_171 = buffer.data(hi1 + 171);
    const auto *hi1_172 = buffer.data(hi1 + 172);
    const auto *hi1_173 = buffer.data(hi1 + 173);
    const auto *hi1_174 = buffer.data(hi1 + 174);
    const auto *hi1_175 = buffer.data(hi1 + 175);
    const auto *hi1_176 = buffer.data(hi1 + 176);
    const auto *hi1_177 = buffer.data(hi1 + 177);
    const auto *hi1_178 = buffer.data(hi1 + 178);
    const auto *hi1_179 = buffer.data(hi1 + 179);
    const auto *hi1_180 = buffer.data(hi1 + 180);
    const auto *hi1_181 = buffer.data(hi1 + 181);
    const auto *hi1_182 = buffer.data(hi1 + 182);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_116 = buffer.data(hk + 116);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_121 = buffer.data(hk + 121);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gk_0, hi0_0, hi0_1, hi1_0, \
                         hi1_1, hk_0, hk_1, hk_2, hk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = f_3 * hi0_0[k]
                 - f_4 * hi1_0[k]
                 + pb_y[k] * hk_1[k];

        t_2[k] = f_3 * hi0_0[k]
                 - f_4 * hi1_0[k]
                 + pb_z[k] * hk_2[k];

        t_3[k] = f_5 * hi0_1[k]
                 - f_6 * hi1_1[k]
                 + pb_y[k] * hk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, hi0_2, hi0_3, hi0_4, hi1_2, hi1_3, \
                         hi1_4, hk_4, hk_5, hk_6, hk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hi0_2[k]
                 - f_6 * hi1_2[k]
                 + pb_z[k] * hk_4[k];

        t_5[k] = f_7 * hi0_3[k]
                 - f_8 * hi1_3[k]
                 + pb_y[k] * hk_5[k];

        t_6[k] = f_3 * hi0_4[k]
                 - f_4 * hi1_4[k]
                 + pb_y[k] * hk_6[k];

        t_7[k] = f_7 * hi0_4[k]
                 - f_8 * hi1_4[k]
                 + pb_z[k] * hk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, hi0_5, hi0_7, hi0_8, hi1_5, hi1_6, \
                         hi1_7, hk_8, hk_9, hk_10, hk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * hi0_5[k]
                 - f_10 * hi1_5[k]
                 + pb_y[k] * hk_8[k];

        t_9[k] = f_5 * hi0_7[k]
                 - f_6 * hi1_6[k]
                 + pb_y[k] * hk_9[k];

        t_10[k] = f_3 * hi0_8[k]
                  - f_4 * hi1_7[k]
                  + pb_y[k] * hk_10[k];

        t_11[k] = f_9 * hi0_8[k]
                  - f_10 * hi1_7[k]
                  + pb_z[k] * hk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, hi0_9, hi0_11, hi0_12, hi1_8, hi1_9, hi1_10, \
                         hk_12, hk_13, hk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * hi0_9[k]
                  - f_12 * hi1_8[k]
                  + pb_y[k] * hk_12[k];

        t_13[k] = f_7 * hi0_11[k]
                  - f_8 * hi1_9[k]
                  + pb_y[k] * hk_13[k];

        t_14[k] = f_5 * hi0_12[k]
                  - f_6 * hi1_10[k]
                  + pb_y[k] * hk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, hi0_13, hi0_14, hi0_16, hi1_11, \
                         hi1_12, hi1_14, hk_15, hk_16, hk_17, hk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * hi0_13[k]
                  - f_4 * hi1_11[k]
                  + pb_y[k] * hk_15[k];

        t_16[k] = f_11 * hi0_13[k]
                  - f_12 * hi1_11[k]
                  + pb_z[k] * hk_16[k];

        t_17[k] = f_1 * hi0_14[k]
                  - f_2 * hi1_12[k]
                  + pb_y[k] * hk_17[k];

        t_18[k] = f_11 * hi0_16[k]
                  - f_12 * hi1_14[k]
                  + pb_y[k] * hk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, hi0_17, hi0_18, hi0_19, hi1_15, hi1_16, \
                         hi1_17, hk_19, hk_20, hk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_9 * hi0_17[k]
                  - f_10 * hi1_15[k]
                  + pb_y[k] * hk_19[k];

        t_20[k] = f_7 * hi0_18[k]
                  - f_8 * hi1_16[k]
                  + pb_y[k] * hk_20[k];

        t_21[k] = f_5 * hi0_19[k]
                  - f_6 * hi1_17[k]
                  + pb_y[k] * hk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pa_z, pb_y, pb_z, gl_0, hi0_20, hi1_18, \
                         hk_22, hk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * hi0_20[k]
                  - f_4 * hi1_18[k]
                  + pb_y[k] * hk_22[k];

        t_23[k] = f_1 * hi0_20[k]
                  - f_2 * hi1_18[k]
                  + pb_z[k] * hk_23[k];

        t_24[k] = pa_y[k] * gl_0[k];

        t_25[k] = pa_z[k] * gl_0[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, fl0_0, fl1_0, gk_18, gk_19, gl_1, \
                         hi0_33, hi0_35, hi1_25, hi1_27, hk_27, hk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_13 * fl0_0[k]
                  - f_14 * fl1_0[k]
                  + pa_y[k] * gl_1[k];

        t_27[k] = f_15 * gk_18[k]
                  + f_11 * hi0_33[k]
                  - f_12 * hi1_25[k]
                  + pb_x[k] * hk_27[k];

        t_28[k] = f_15 * gk_19[k]
                  + f_9 * hi0_35[k]
                  - f_10 * hi1_27[k]
                  + pb_x[k] * hk_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, gk_20, gk_21, gk_22, hi0_37, hi0_39, hi0_40, \
                         hi1_29, hi1_31, hi1_32, hk_29, hk_30, hk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_15 * gk_20[k]
                  + f_7 * hi0_37[k]
                  - f_8 * hi1_29[k]
                  + pb_x[k] * hk_29[k];

        t_30[k] = f_15 * gk_21[k]
                  + f_5 * hi0_39[k]
                  - f_6 * hi1_31[k]
                  + pb_x[k] * hk_30[k];

        t_31[k] = f_15 * gk_22[k]
                  + f_3 * hi0_40[k]
                  - f_4 * hi1_32[k]
                  + pb_x[k] * hk_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pa_z, pb_x, fl0_0, fl0_3, fl1_0, fl1_3, \
                         gk_25, gl_2, gl_9, hi0_51, hi1_41, hk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_16 * fl0_3[k]
                  - f_17 * fl1_3[k]
                  + pa_x[k] * gl_9[k];

        t_33[k] = f_13 * fl0_0[k]
                  - f_14 * fl1_0[k]
                  + pa_z[k] * gl_2[k];

        t_34[k] = f_15 * gk_25[k]
                  + f_11 * hi0_51[k]
                  - f_12 * hi1_41[k]
                  + pb_x[k] * hk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, gk_26, gk_27, gk_28, hi0_53, hi0_55, hi0_56, \
                         hi1_43, hi1_45, hi1_46, hk_35, hk_36, hk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_15 * gk_26[k]
                  + f_9 * hi0_53[k]
                  - f_10 * hi1_43[k]
                  + pb_x[k] * hk_35[k];

        t_36[k] = f_15 * gk_27[k]
                  + f_7 * hi0_55[k]
                  - f_8 * hi1_45[k]
                  + pb_x[k] * hk_36[k];

        t_37[k] = f_15 * gk_28[k]
                  + f_5 * hi0_56[k]
                  - f_6 * hi1_46[k]
                  + pb_x[k] * hk_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pa_y, pb_x, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gk_29, gl_3, gl_16, hi0_62, hi1_52, hk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_15 * gk_29[k]
                  + f_3 * hi0_62[k]
                  - f_4 * hi1_52[k]
                  + pb_x[k] * hk_38[k];

        t_39[k] = f_16 * fl0_4[k]
                  - f_17 * fl1_4[k]
                  + pa_x[k] * gl_16[k];

        t_40[k] = f_16 * fl0_1[k]
                  - f_17 * fl1_1[k]
                  + pa_y[k] * gl_3[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_x, gk_31, gk_32, gk_33, hi0_65, hi0_67, hi0_69, \
                         hi1_55, hi1_57, hi1_59, hk_41, hk_42, hk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_18 * gk_31[k]
                  + f_11 * hi0_65[k]
                  - f_12 * hi1_55[k]
                  + pb_x[k] * hk_41[k];

        t_42[k] = f_18 * gk_32[k]
                  + f_9 * hi0_67[k]
                  - f_10 * hi1_57[k]
                  + pb_x[k] * hk_42[k];

        t_43[k] = f_18 * gk_33[k]
                  + f_7 * hi0_69[k]
                  - f_8 * hi1_59[k]
                  + pb_x[k] * hk_43[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_x, pb_x, fl0_5, fl1_5, gk_34, gk_35, gl_17, \
                         hi0_71, hi0_72, hi1_61, hi1_62, hk_44, hk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_18 * gk_34[k]
                  + f_5 * hi0_71[k]
                  - f_6 * hi1_61[k]
                  + pb_x[k] * hk_44[k];

        t_45[k] = f_18 * gk_35[k]
                  + f_3 * hi0_72[k]
                  - f_4 * hi1_62[k]
                  + pb_x[k] * hk_45[k];

        t_46[k] = f_13 * fl0_5[k]
                  - f_14 * fl1_5[k]
                  + pa_x[k] * gl_17[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, pa_y, pa_z, gl_4, gl_5, \
                         gl_6, gl_7, gl_8, gl_10, gl_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_z[k] * gl_4[k];

        t_48[k] = pa_z[k] * gl_5[k];

        t_49[k] = pa_z[k] * gl_6[k];

        t_50[k] = pa_z[k] * gl_7[k];

        t_51[k] = pa_z[k] * gl_8[k];

        t_52[k] = pa_y[k] * gl_10[k];

        t_53[k] = pa_y[k] * gl_11[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_y, pa_z, fl0_2, fl1_2, gl_10, gl_12, \
                         gl_13, gl_14, gl_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_y[k] * gl_12[k];

        t_55[k] = pa_y[k] * gl_13[k];

        t_56[k] = pa_y[k] * gl_14[k];

        t_57[k] = pa_y[k] * gl_15[k];

        t_58[k] = f_16 * fl0_2[k]
                  - f_17 * fl1_2[k]
                  + pa_z[k] * gl_10[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_x, gk_37, gk_38, gk_39, hi0_87, hi0_89, hi0_91, \
                         hi1_72, hi1_74, hi1_76, hk_57, hk_58, hk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_18 * gk_37[k]
                  + f_11 * hi0_87[k]
                  - f_12 * hi1_72[k]
                  + pb_x[k] * hk_57[k];

        t_60[k] = f_18 * gk_38[k]
                  + f_9 * hi0_89[k]
                  - f_10 * hi1_74[k]
                  + pb_x[k] * hk_58[k];

        t_61[k] = f_18 * gk_39[k]
                  + f_7 * hi0_91[k]
                  - f_8 * hi1_76[k]
                  + pb_x[k] * hk_59[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_x, pb_x, fl0_8, fl1_8, gk_40, gk_41, gl_18, \
                         hi0_92, hi0_98, hi1_77, hi1_83, hk_60, hk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_18 * gk_40[k]
                  + f_5 * hi0_92[k]
                  - f_6 * hi1_77[k]
                  + pb_x[k] * hk_60[k];

        t_63[k] = f_18 * gk_41[k]
                  + f_3 * hi0_98[k]
                  - f_4 * hi1_83[k]
                  + pb_x[k] * hk_61[k];

        t_64[k] = f_13 * fl0_8[k]
                  - f_14 * fl1_8[k]
                  + pa_x[k] * gl_18[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, t_71, pa_x, gl_19, gl_21, gl_22, \
                         gl_23, gl_24, gl_25, gl_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_x[k] * gl_19[k];

        t_66[k] = pa_x[k] * gl_21[k];

        t_67[k] = pa_x[k] * gl_22[k];

        t_68[k] = pa_x[k] * gl_23[k];

        t_69[k] = pa_x[k] * gl_24[k];

        t_70[k] = pa_x[k] * gl_25[k];

        t_71[k] = pa_x[k] * gl_26[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_x, pb_x, gl_27, gl_29, hi0_121, hi0_123, \
                         hi1_101, hi1_103, hk_71, hk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_x[k] * gl_27[k];

        t_73[k] = pa_x[k] * gl_29[k];

        t_74[k] = f_1 * hi0_121[k]
                  - f_2 * hi1_101[k]
                  + pb_x[k] * hk_71[k];

        t_75[k] = f_11 * hi0_123[k]
                  - f_12 * hi1_103[k]
                  + pb_x[k] * hk_72[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_x, hi0_124, hi0_125, hi0_127, hi1_104, hi1_105, \
                         hi1_106, hk_73, hk_74, hk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_11 * hi0_124[k]
                  - f_12 * hi1_104[k]
                  + pb_x[k] * hk_73[k];

        t_77[k] = f_9 * hi0_125[k]
                  - f_10 * hi1_105[k]
                  + pb_x[k] * hk_74[k];

        t_78[k] = f_9 * hi0_127[k]
                  - f_10 * hi1_106[k]
                  + pb_x[k] * hk_75[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_x, hi0_128, hi0_130, hi0_131, hi1_107, hi1_108, \
                         hi1_109, hk_76, hk_77, hk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_7 * hi0_128[k]
                  - f_8 * hi1_107[k]
                  + pb_x[k] * hk_76[k];

        t_80[k] = f_7 * hi0_130[k]
                  - f_8 * hi1_108[k]
                  + pb_x[k] * hk_77[k];

        t_81[k] = f_7 * hi0_131[k]
                  - f_8 * hi1_109[k]
                  + pb_x[k] * hk_78[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pb_x, hi0_132, hi0_133, hi0_134, hi1_110, hi1_111, \
                         hi1_112, hk_79, hk_80, hk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_5 * hi0_132[k]
                  - f_6 * hi1_110[k]
                  + pb_x[k] * hk_79[k];

        t_83[k] = f_5 * hi0_133[k]
                  - f_6 * hi1_111[k]
                  + pb_x[k] * hk_80[k];

        t_84[k] = f_5 * hi0_134[k]
                  - f_6 * hi1_112[k]
                  + pb_x[k] * hk_81[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pb_x, hi0_135, hi0_136, hi0_138, hi1_113, hi1_114, \
                         hi1_116, hk_82, hk_83, hk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_5 * hi0_135[k]
                  - f_6 * hi1_113[k]
                  + pb_x[k] * hk_82[k];

        t_86[k] = f_3 * hi0_136[k]
                  - f_4 * hi1_114[k]
                  + pb_x[k] * hk_83[k];

        t_87[k] = f_3 * hi0_138[k]
                  - f_4 * hi1_116[k]
                  + pb_x[k] * hk_84[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pb_x, hi0_139, hi0_140, hi0_141, hi1_117, hi1_118, \
                         hi1_119, hk_85, hk_86, hk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * hi0_139[k]
                  - f_4 * hi1_117[k]
                  + pb_x[k] * hk_85[k];

        t_89[k] = f_3 * hi0_140[k]
                  - f_4 * hi1_118[k]
                  + pb_x[k] * hk_86[k];

        t_90[k] = f_3 * hi0_141[k]
                  - f_4 * hi1_119[k]
                  + pb_x[k] * hk_87[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_y, pb_z, gk_52, hi0_136, hi0_137, hi1_114, \
                         hi1_115, hk_88, hk_89, hk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_0 * gk_52[k]
                  + f_1 * hi0_136[k]
                  - f_2 * hi1_114[k]
                  + pb_y[k] * hk_88[k];

        t_92[k] = f_3 * hi0_136[k]
                  - f_4 * hi1_114[k]
                  + pb_z[k] * hk_89[k];

        t_93[k] = f_5 * hi0_137[k]
                  - f_6 * hi1_115[k]
                  + pb_z[k] * hk_90[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_z, hi0_138, hi0_139, hi0_140, hi1_116, hi1_117, \
                         hi1_118, hk_91, hk_92, hk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_7 * hi0_138[k]
                  - f_8 * hi1_116[k]
                  + pb_z[k] * hk_91[k];

        t_95[k] = f_9 * hi0_139[k]
                  - f_10 * hi1_117[k]
                  + pb_z[k] * hk_92[k];

        t_96[k] = f_11 * hi0_140[k]
                  - f_12 * hi1_118[k]
                  + pb_z[k] * hk_93[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pa_z, pb_z, fl0_5, fl1_5, gl_19, gl_20, hi0_141, \
                         hi1_119, hk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * hi0_141[k]
                  - f_2 * hi1_119[k]
                  + pb_z[k] * hk_94[k];

        t_98[k] = pa_z[k] * gl_19[k];

        t_99[k] = f_13 * fl0_5[k]
                  - f_14 * fl1_5[k]
                  + pa_z[k] * gl_20[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pb_y, gk_60, gk_61, gk_62, hi0_162, hi0_163, \
                         hi0_164, hi1_135, hi1_136, hi1_137, hk_97, hk_98, \
                         hk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_15 * gk_60[k]
                   + f_11 * hi0_162[k]
                   - f_12 * hi1_135[k]
                   + pb_y[k] * hk_97[k];

        t_101[k] = f_15 * gk_61[k]
                   + f_9 * hi0_163[k]
                   - f_10 * hi1_136[k]
                   + pb_y[k] * hk_98[k];

        t_102[k] = f_15 * gk_62[k]
                   + f_7 * hi0_164[k]
                   - f_8 * hi1_137[k]
                   + pb_y[k] * hk_99[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_y, pb_y, fl0_7, fl1_7, gk_63, gk_64, gl_27, \
                         hi0_165, hi0_166, hi1_138, hi1_139, hk_100, \
                         hk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_15 * gk_63[k]
                   + f_5 * hi0_165[k]
                   - f_6 * hi1_138[k]
                   + pb_y[k] * hk_100[k];

        t_104[k] = f_15 * gk_64[k]
                   + f_3 * hi0_166[k]
                   - f_4 * hi1_139[k]
                   + pb_y[k] * hk_101[k];

        t_105[k] = f_16 * fl0_7[k]
                   - f_17 * fl1_7[k]
                   + pa_y[k] * gl_27[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pa_z, pb_y, fl0_6, fl1_6, gk_66, gk_67, gl_21, \
                         hi0_180, hi0_181, hi1_153, hi1_154, hk_104, \
                         hk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_16 * fl0_6[k]
                   - f_17 * fl1_6[k]
                   + pa_z[k] * gl_21[k];

        t_107[k] = f_18 * gk_66[k]
                   + f_11 * hi0_180[k]
                   - f_12 * hi1_153[k]
                   + pb_y[k] * hk_104[k];

        t_108[k] = f_18 * gk_67[k]
                   + f_9 * hi0_181[k]
                   - f_10 * hi1_154[k]
                   + pb_y[k] * hk_105[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, pb_y, gk_68, gk_69, gk_70, hi0_182, hi0_183, \
                         hi0_184, hi1_155, hi1_156, hi1_157, hk_106, hk_107, \
                         hk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_18 * gk_68[k]
                   + f_7 * hi0_182[k]
                   - f_8 * hi1_155[k]
                   + pb_y[k] * hk_106[k];

        t_110[k] = f_18 * gk_69[k]
                   + f_5 * hi0_183[k]
                   - f_6 * hi1_156[k]
                   + pb_y[k] * hk_107[k];

        t_111[k] = f_18 * gk_70[k]
                   + f_3 * hi0_184[k]
                   - f_4 * hi1_157[k]
                   + pb_y[k] * hk_108[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_y, pb_x, fl0_8, fl1_8, gl_28, gl_29, \
                         hi0_195, hi0_197, hi1_164, hi1_166, hk_111, \
                         hk_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_13 * fl0_8[k]
                   - f_14 * fl1_8[k]
                   + pa_y[k] * gl_28[k];

        t_113[k] = pa_y[k] * gl_29[k];

        t_114[k] = f_1 * hi0_195[k]
                   - f_2 * hi1_164[k]
                   + pb_x[k] * hk_111[k];

        t_115[k] = f_11 * hi0_197[k]
                   - f_12 * hi1_166[k]
                   + pb_x[k] * hk_112[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_x, hi0_198, hi0_199, hi0_201, hi1_167, \
                         hi1_168, hi1_169, hk_113, hk_114, hk_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * hi0_198[k]
                   - f_12 * hi1_167[k]
                   + pb_x[k] * hk_113[k];

        t_117[k] = f_9 * hi0_199[k]
                   - f_10 * hi1_168[k]
                   + pb_x[k] * hk_114[k];

        t_118[k] = f_9 * hi0_201[k]
                   - f_10 * hi1_169[k]
                   + pb_x[k] * hk_115[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pb_x, hi0_202, hi0_203, hi0_205, hi1_170, \
                         hi1_171, hi1_172, hk_116, hk_117, hk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_7 * hi0_202[k]
                   - f_8 * hi1_170[k]
                   + pb_x[k] * hk_116[k];

        t_120[k] = f_7 * hi0_203[k]
                   - f_8 * hi1_171[k]
                   + pb_x[k] * hk_117[k];

        t_121[k] = f_7 * hi0_205[k]
                   - f_8 * hi1_172[k]
                   + pb_x[k] * hk_118[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pb_x, hi0_206, hi0_207, hi0_208, hi1_173, \
                         hi1_174, hi1_175, hk_119, hk_120, hk_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_5 * hi0_206[k]
                   - f_6 * hi1_173[k]
                   + pb_x[k] * hk_119[k];

        t_123[k] = f_5 * hi0_207[k]
                   - f_6 * hi1_174[k]
                   + pb_x[k] * hk_120[k];

        t_124[k] = f_5 * hi0_208[k]
                   - f_6 * hi1_175[k]
                   + pb_x[k] * hk_121[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pb_x, hi0_209, hi0_210, hi0_211, hi1_176, \
                         hi1_177, hi1_178, hk_122, hk_123, hk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_5 * hi0_209[k]
                   - f_6 * hi1_176[k]
                   + pb_x[k] * hk_122[k];

        t_126[k] = f_3 * hi0_210[k]
                   - f_4 * hi1_177[k]
                   + pb_x[k] * hk_123[k];

        t_127[k] = f_3 * hi0_211[k]
                   - f_4 * hi1_178[k]
                   + pb_x[k] * hk_124[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, hi0_212, hi0_213, hi0_215, hi1_179, \
                         hi1_180, hi1_182, hk_125, hk_126, hk_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * hi0_212[k]
                   - f_4 * hi1_179[k]
                   + pb_x[k] * hk_125[k];

        t_129[k] = f_3 * hi0_213[k]
                   - f_4 * hi1_180[k]
                   + pb_x[k] * hk_126[k];

        t_130[k] = f_3 * hi0_215[k]
                   - f_4 * hi1_182[k]
                   + pb_x[k] * hk_127[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pb_y, hi0_210, hi0_211, hi0_212, hi1_177, \
                         hi1_178, hi1_179, hk_128, hk_129, hk_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_1 * hi0_210[k]
                   - f_2 * hi1_177[k]
                   + pb_y[k] * hk_128[k];

        t_132[k] = f_11 * hi0_211[k]
                   - f_12 * hi1_178[k]
                   + pb_y[k] * hk_129[k];

        t_133[k] = f_9 * hi0_212[k]
                   - f_10 * hi1_179[k]
                   + pb_y[k] * hk_130[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pb_y, hi0_213, hi0_214, hi0_215, hi1_180, \
                         hi1_181, hi1_182, hk_131, hk_132, hk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_7 * hi0_213[k]
                   - f_8 * hi1_180[k]
                   + pb_y[k] * hk_131[k];

        t_135[k] = f_5 * hi0_214[k]
                   - f_6 * hi1_181[k]
                   + pb_y[k] * hk_132[k];

        t_136[k] = f_3 * hi0_215[k]
                   - f_4 * hi1_182[k]
                   + pb_y[k] * hk_133[k];
    }

#pragma omp simd aligned(t_137, pb_z, gk_86, hi0_215, hi1_182, hk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_0 * gk_86[k]
                   + f_1 * hi0_215[k]
                   - f_2 * hi1_182[k]
                   + pb_z[k] * hk_134[k];
    }
}

auto
compute_prim_hl_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
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
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);
    const auto f_18 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_56 = buffer.data(gk + 56);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_6 = buffer.data(hi0 + 6);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_10 = buffer.data(hi0 + 10);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_18 = buffer.data(hi0 + 18);
    const auto *hi0_19 = buffer.data(hi0 + 19);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_22 = buffer.data(hi0 + 22);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_25 = buffer.data(hi0 + 25);
    const auto *hi0_26 = buffer.data(hi0 + 26);
    const auto *hi0_32 = buffer.data(hi0 + 32);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_36 = buffer.data(hi0 + 36);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_38 = buffer.data(hi0 + 38);
    const auto *hi0_39 = buffer.data(hi0 + 39);
    const auto *hi0_41 = buffer.data(hi0 + 41);
    const auto *hi0_42 = buffer.data(hi0 + 42);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_44 = buffer.data(hi0 + 44);
    const auto *hi0_45 = buffer.data(hi0 + 45);
    const auto *hi0_47 = buffer.data(hi0 + 47);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_4 = buffer.data(hi1 + 4);
    const auto *hi1_5 = buffer.data(hi1 + 5);
    const auto *hi1_6 = buffer.data(hi1 + 6);
    const auto *hi1_7 = buffer.data(hi1 + 7);
    const auto *hi1_8 = buffer.data(hi1 + 8);
    const auto *hi1_10 = buffer.data(hi1 + 10);
    const auto *hi1_11 = buffer.data(hi1 + 11);
    const auto *hi1_12 = buffer.data(hi1 + 12);
    const auto *hi1_13 = buffer.data(hi1 + 13);
    const auto *hi1_14 = buffer.data(hi1 + 14);
    const auto *hi1_16 = buffer.data(hi1 + 16);
    const auto *hi1_17 = buffer.data(hi1 + 17);
    const auto *hi1_18 = buffer.data(hi1 + 18);
    const auto *hi1_19 = buffer.data(hi1 + 19);
    const auto *hi1_20 = buffer.data(hi1 + 20);
    const auto *hi1_23 = buffer.data(hi1 + 23);
    const auto *hi1_24 = buffer.data(hi1 + 24);
    const auto *hi1_25 = buffer.data(hi1 + 25);
    const auto *hi1_26 = buffer.data(hi1 + 26);
    const auto *hi1_27 = buffer.data(hi1 + 27);
    const auto *hi1_43 = buffer.data(hi1 + 43);
    const auto *hi1_46 = buffer.data(hi1 + 46);
    const auto *hi1_47 = buffer.data(hi1 + 47);
    const auto *hi1_48 = buffer.data(hi1 + 48);
    const auto *hi1_49 = buffer.data(hi1 + 49);
    const auto *hi1_50 = buffer.data(hi1 + 50);
    const auto *hi1_52 = buffer.data(hi1 + 52);
    const auto *hi1_53 = buffer.data(hi1 + 53);
    const auto *hi1_54 = buffer.data(hi1 + 54);
    const auto *hi1_55 = buffer.data(hi1 + 55);
    const auto *hi1_56 = buffer.data(hi1 + 56);
    const auto *hi1_62 = buffer.data(hi1 + 62);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_80 = buffer.data(hk + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, fl0_0, fl1_0, gk_0, gl_0, gl_1, \
                         hi0_0, hi1_0, hk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = pa_y[k] * gl_0[k];

        t_2[k] = pa_z[k] * gl_0[k];

        t_3[k] = f_3 * fl0_0[k]
                 - f_4 * fl1_0[k]
                 + pa_y[k] * gl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, gk_4, gk_5, gk_6, hi0_4, hi0_5, hi0_6, hi1_4, \
                         hi1_5, hi1_6, hk_4, hk_5, hk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gk_4[k]
                 + f_6 * hi0_4[k]
                 - f_7 * hi1_4[k]
                 + pb_x[k] * hk_4[k];

        t_5[k] = f_5 * gk_5[k]
                 + f_8 * hi0_5[k]
                 - f_9 * hi1_5[k]
                 + pb_x[k] * hk_5[k];

        t_6[k] = f_5 * gk_6[k]
                 + f_10 * hi0_6[k]
                 - f_11 * hi1_6[k]
                 + pb_x[k] * hk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, fl0_3, fl1_3, gk_7, gk_8, gl_9, hi0_7, \
                         hi0_8, hi1_7, hi1_8, hk_7, hk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * gk_7[k]
                 + f_12 * hi0_7[k]
                 - f_13 * hi1_7[k]
                 + pb_x[k] * hk_7[k];

        t_8[k] = f_5 * gk_8[k]
                 + f_14 * hi0_8[k]
                 - f_15 * hi1_8[k]
                 + pb_x[k] * hk_8[k];

        t_9[k] = f_16 * fl0_3[k]
                 - f_17 * fl1_3[k]
                 + pa_x[k] * gl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, fl0_0, fl1_0, gk_11, gk_12, gl_2, \
                         hi0_10, hi0_11, hi1_10, hi1_11, hk_11, hk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fl0_0[k]
                  - f_4 * fl1_0[k]
                  + pa_z[k] * gl_2[k];

        t_11[k] = f_5 * gk_11[k]
                  + f_6 * hi0_10[k]
                  - f_7 * hi1_10[k]
                  + pb_x[k] * hk_11[k];

        t_12[k] = f_5 * gk_12[k]
                  + f_8 * hi0_11[k]
                  - f_9 * hi1_11[k]
                  + pb_x[k] * hk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, gk_13, gk_14, gk_15, hi0_12, hi0_13, hi0_14, \
                         hi1_12, hi1_13, hi1_14, hk_13, hk_14, hk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gk_13[k]
                  + f_10 * hi0_12[k]
                  - f_11 * hi1_12[k]
                  + pb_x[k] * hk_13[k];

        t_14[k] = f_5 * gk_14[k]
                  + f_12 * hi0_13[k]
                  - f_13 * hi1_13[k]
                  + pb_x[k] * hk_14[k];

        t_15[k] = f_5 * gk_15[k]
                  + f_14 * hi0_14[k]
                  - f_15 * hi1_14[k]
                  + pb_x[k] * hk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_x, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gk_17, gl_3, gl_16, hi0_16, hi1_16, hk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_16 * fl0_4[k]
                  - f_17 * fl1_4[k]
                  + pa_x[k] * gl_16[k];

        t_17[k] = f_16 * fl0_1[k]
                  - f_17 * fl1_1[k]
                  + pa_y[k] * gl_3[k];

        t_18[k] = f_18 * gk_17[k]
                  + f_6 * hi0_16[k]
                  - f_7 * hi1_16[k]
                  + pb_x[k] * hk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gk_18, gk_19, gk_20, hi0_17, hi0_18, hi0_19, \
                         hi1_17, hi1_18, hi1_19, hk_19, hk_20, hk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_18 * gk_18[k]
                  + f_8 * hi0_17[k]
                  - f_9 * hi1_17[k]
                  + pb_x[k] * hk_19[k];

        t_20[k] = f_18 * gk_19[k]
                  + f_10 * hi0_18[k]
                  - f_11 * hi1_18[k]
                  + pb_x[k] * hk_20[k];

        t_21[k] = f_18 * gk_20[k]
                  + f_12 * hi0_19[k]
                  - f_13 * hi1_19[k]
                  + pb_x[k] * hk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, pb_x, fl0_5, fl1_5, gk_21, gl_4, \
                         gl_5, gl_17, hi0_20, hi1_20, hk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_18 * gk_21[k]
                  + f_14 * hi0_20[k]
                  - f_15 * hi1_20[k]
                  + pb_x[k] * hk_22[k];

        t_23[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_x[k] * gl_17[k];

        t_24[k] = pa_z[k] * gl_4[k];

        t_25[k] = pa_z[k] * gl_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, gl_6, gl_7, \
                         gl_8, gl_10, gl_11, gl_12, gl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * gl_6[k];

        t_27[k] = pa_z[k] * gl_7[k];

        t_28[k] = pa_z[k] * gl_8[k];

        t_29[k] = pa_y[k] * gl_10[k];

        t_30[k] = pa_y[k] * gl_11[k];

        t_31[k] = pa_y[k] * gl_12[k];

        t_32[k] = pa_y[k] * gl_13[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pa_z, pb_x, fl0_2, fl1_2, gk_23, gl_10, \
                         gl_14, gl_15, hi0_22, hi1_23, hk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * gl_14[k];

        t_34[k] = pa_y[k] * gl_15[k];

        t_35[k] = f_16 * fl0_2[k]
                  - f_17 * fl1_2[k]
                  + pa_z[k] * gl_10[k];

        t_36[k] = f_18 * gk_23[k]
                  + f_6 * hi0_22[k]
                  - f_7 * hi1_23[k]
                  + pb_x[k] * hk_34[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, gk_24, gk_25, gk_26, hi0_23, hi0_24, hi0_25, \
                         hi1_24, hi1_25, hi1_26, hk_35, hk_36, hk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_18 * gk_24[k]
                  + f_8 * hi0_23[k]
                  - f_9 * hi1_24[k]
                  + pb_x[k] * hk_35[k];

        t_38[k] = f_18 * gk_25[k]
                  + f_10 * hi0_24[k]
                  - f_11 * hi1_25[k]
                  + pb_x[k] * hk_36[k];

        t_39[k] = f_18 * gk_26[k]
                  + f_12 * hi0_25[k]
                  - f_13 * hi1_26[k]
                  + pb_x[k] * hk_37[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, fl0_8, fl1_8, gk_27, gl_18, \
                         gl_19, gl_21, hi0_26, hi1_27, hk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_18 * gk_27[k]
                  + f_14 * hi0_26[k]
                  - f_15 * hi1_27[k]
                  + pb_x[k] * hk_38[k];

        t_41[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_x[k] * gl_18[k];

        t_42[k] = pa_x[k] * gl_19[k];

        t_43[k] = pa_x[k] * gl_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, pa_x, gl_22, gl_23, gl_24, \
                         gl_25, gl_26, gl_27, gl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gl_22[k];

        t_45[k] = pa_x[k] * gl_23[k];

        t_46[k] = pa_x[k] * gl_24[k];

        t_47[k] = pa_x[k] * gl_25[k];

        t_48[k] = pa_x[k] * gl_26[k];

        t_49[k] = pa_x[k] * gl_27[k];

        t_50[k] = pa_x[k] * gl_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, fl0_5, fl1_5, gk_33, gl_19, gl_20, \
                         hi0_32, hi1_43, hk_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * gk_33[k]
                  + f_1 * hi0_32[k]
                  - f_2 * hi1_43[k]
                  + pb_y[k] * hk_58[k];

        t_52[k] = pa_z[k] * gl_19[k];

        t_53[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_z[k] * gl_20[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, gk_36, gk_37, gk_38, hi0_35, hi0_36, hi0_37, \
                         hi1_46, hi1_47, hi1_48, hk_61, hk_62, hk_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * gk_36[k]
                  + f_6 * hi0_35[k]
                  - f_7 * hi1_46[k]
                  + pb_y[k] * hk_61[k];

        t_55[k] = f_5 * gk_37[k]
                  + f_8 * hi0_36[k]
                  - f_9 * hi1_47[k]
                  + pb_y[k] * hk_62[k];

        t_56[k] = f_5 * gk_38[k]
                  + f_10 * hi0_37[k]
                  - f_11 * hi1_48[k]
                  + pb_y[k] * hk_63[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pb_y, fl0_7, fl1_7, gk_39, gk_40, gl_27, \
                         hi0_38, hi0_39, hi1_49, hi1_50, hk_64, hk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * gk_39[k]
                  + f_12 * hi0_38[k]
                  - f_13 * hi1_49[k]
                  + pb_y[k] * hk_64[k];

        t_58[k] = f_5 * gk_40[k]
                  + f_14 * hi0_39[k]
                  - f_15 * hi1_50[k]
                  + pb_y[k] * hk_65[k];

        t_59[k] = f_16 * fl0_7[k]
                  - f_17 * fl1_7[k]
                  + pa_y[k] * gl_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_z, pb_y, fl0_6, fl1_6, gk_42, gk_43, gl_21, \
                         hi0_41, hi0_42, hi1_52, hi1_53, hk_68, hk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_16 * fl0_6[k]
                  - f_17 * fl1_6[k]
                  + pa_z[k] * gl_21[k];

        t_61[k] = f_18 * gk_42[k]
                  + f_6 * hi0_41[k]
                  - f_7 * hi1_52[k]
                  + pb_y[k] * hk_68[k];

        t_62[k] = f_18 * gk_43[k]
                  + f_8 * hi0_42[k]
                  - f_9 * hi1_53[k]
                  + pb_y[k] * hk_69[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, gk_44, gk_45, gk_46, hi0_43, hi0_44, hi0_45, \
                         hi1_54, hi1_55, hi1_56, hk_70, hk_71, hk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_18 * gk_44[k]
                  + f_10 * hi0_43[k]
                  - f_11 * hi1_54[k]
                  + pb_y[k] * hk_70[k];

        t_64[k] = f_18 * gk_45[k]
                  + f_12 * hi0_44[k]
                  - f_13 * hi1_55[k]
                  + pb_y[k] * hk_71[k];

        t_65[k] = f_18 * gk_46[k]
                  + f_14 * hi0_45[k]
                  - f_15 * hi1_56[k]
                  + pb_y[k] * hk_72[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pb_z, fl0_8, fl1_8, gk_56, gl_28, gl_29, \
                         hi0_47, hi1_62, hk_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_y[k] * gl_28[k];

        t_67[k] = pa_y[k] * gl_29[k];

        t_68[k] = f_0 * gk_56[k]
                  + f_1 * hi0_47[k]
                  - f_2 * hi1_62[k]
                  + pb_z[k] * hk_80[k];
    }
}

auto
compute_prim_hl_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
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
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);
    const auto f_18 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_56 = buffer.data(gk + 56);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_4 = buffer.data(hi0 + 4);
    const auto *hi0_5 = buffer.data(hi0 + 5);
    const auto *hi0_6 = buffer.data(hi0 + 6);
    const auto *hi0_7 = buffer.data(hi0 + 7);
    const auto *hi0_8 = buffer.data(hi0 + 8);
    const auto *hi0_10 = buffer.data(hi0 + 10);
    const auto *hi0_11 = buffer.data(hi0 + 11);
    const auto *hi0_12 = buffer.data(hi0 + 12);
    const auto *hi0_13 = buffer.data(hi0 + 13);
    const auto *hi0_14 = buffer.data(hi0 + 14);
    const auto *hi0_16 = buffer.data(hi0 + 16);
    const auto *hi0_17 = buffer.data(hi0 + 17);
    const auto *hi0_18 = buffer.data(hi0 + 18);
    const auto *hi0_19 = buffer.data(hi0 + 19);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_25 = buffer.data(hi0 + 25);
    const auto *hi0_26 = buffer.data(hi0 + 26);
    const auto *hi0_27 = buffer.data(hi0 + 27);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_46 = buffer.data(hi0 + 46);
    const auto *hi0_47 = buffer.data(hi0 + 47);
    const auto *hi0_48 = buffer.data(hi0 + 48);
    const auto *hi0_49 = buffer.data(hi0 + 49);
    const auto *hi0_50 = buffer.data(hi0 + 50);
    const auto *hi0_52 = buffer.data(hi0 + 52);
    const auto *hi0_53 = buffer.data(hi0 + 53);
    const auto *hi0_54 = buffer.data(hi0 + 54);
    const auto *hi0_55 = buffer.data(hi0 + 55);
    const auto *hi0_56 = buffer.data(hi0 + 56);
    const auto *hi0_62 = buffer.data(hi0 + 62);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_20 = buffer.data(hi1 + 20);
    const auto *hi1_21 = buffer.data(hi1 + 21);
    const auto *hi1_22 = buffer.data(hi1 + 22);
    const auto *hi1_23 = buffer.data(hi1 + 23);
    const auto *hi1_24 = buffer.data(hi1 + 24);
    const auto *hi1_27 = buffer.data(hi1 + 27);
    const auto *hi1_28 = buffer.data(hi1 + 28);
    const auto *hi1_29 = buffer.data(hi1 + 29);
    const auto *hi1_30 = buffer.data(hi1 + 30);
    const auto *hi1_32 = buffer.data(hi1 + 32);
    const auto *hi1_34 = buffer.data(hi1 + 34);
    const auto *hi1_35 = buffer.data(hi1 + 35);
    const auto *hi1_36 = buffer.data(hi1 + 36);
    const auto *hi1_37 = buffer.data(hi1 + 37);
    const auto *hi1_38 = buffer.data(hi1 + 38);
    const auto *hi1_42 = buffer.data(hi1 + 42);
    const auto *hi1_43 = buffer.data(hi1 + 43);
    const auto *hi1_44 = buffer.data(hi1 + 44);
    const auto *hi1_45 = buffer.data(hi1 + 45);
    const auto *hi1_47 = buffer.data(hi1 + 47);
    const auto *hi1_74 = buffer.data(hi1 + 74);
    const auto *hi1_85 = buffer.data(hi1 + 85);
    const auto *hi1_86 = buffer.data(hi1 + 86);
    const auto *hi1_87 = buffer.data(hi1 + 87);
    const auto *hi1_88 = buffer.data(hi1 + 88);
    const auto *hi1_89 = buffer.data(hi1 + 89);
    const auto *hi1_93 = buffer.data(hi1 + 93);
    const auto *hi1_94 = buffer.data(hi1 + 94);
    const auto *hi1_95 = buffer.data(hi1 + 95);
    const auto *hi1_96 = buffer.data(hi1 + 96);
    const auto *hi1_97 = buffer.data(hi1 + 97);
    const auto *hi1_118 = buffer.data(hi1 + 118);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_25 = buffer.data(hk + 25);
    const auto *hk_26 = buffer.data(hk + 26);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_122 = buffer.data(hk + 122);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, fl0_0, fl1_0, gk_0, gl_0, gl_1, \
                         hi0_0, hi1_0, hk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = pa_y[k] * gl_0[k];

        t_2[k] = pa_z[k] * gl_0[k];

        t_3[k] = f_3 * fl0_0[k]
                 - f_4 * fl1_0[k]
                 + pa_y[k] * gl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, gk_4, gk_5, gk_6, hi0_4, hi0_5, hi0_6, hi1_20, \
                         hi1_21, hi1_22, hk_18, hk_19, hk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gk_4[k]
                 + f_6 * hi0_4[k]
                 - f_7 * hi1_20[k]
                 + pb_x[k] * hk_18[k];

        t_5[k] = f_5 * gk_5[k]
                 + f_8 * hi0_5[k]
                 - f_9 * hi1_21[k]
                 + pb_x[k] * hk_19[k];

        t_6[k] = f_5 * gk_6[k]
                 + f_10 * hi0_6[k]
                 - f_11 * hi1_22[k]
                 + pb_x[k] * hk_20[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, fl0_3, fl1_3, gk_7, gk_8, gl_9, hi0_7, \
                         hi0_8, hi1_23, hi1_24, hk_21, hk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * gk_7[k]
                 + f_12 * hi0_7[k]
                 - f_13 * hi1_23[k]
                 + pb_x[k] * hk_21[k];

        t_8[k] = f_5 * gk_8[k]
                 + f_14 * hi0_8[k]
                 - f_15 * hi1_24[k]
                 + pb_x[k] * hk_22[k];

        t_9[k] = f_16 * fl0_3[k]
                 - f_17 * fl1_3[k]
                 + pa_x[k] * gl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, fl0_0, fl1_0, gk_11, gk_12, gl_2, \
                         hi0_10, hi0_11, hi1_27, hi1_28, hk_25, hk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fl0_0[k]
                  - f_4 * fl1_0[k]
                  + pa_z[k] * gl_2[k];

        t_11[k] = f_5 * gk_11[k]
                  + f_6 * hi0_10[k]
                  - f_7 * hi1_27[k]
                  + pb_x[k] * hk_25[k];

        t_12[k] = f_5 * gk_12[k]
                  + f_8 * hi0_11[k]
                  - f_9 * hi1_28[k]
                  + pb_x[k] * hk_26[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, gk_13, gk_14, gk_15, hi0_12, hi0_13, hi0_14, \
                         hi1_29, hi1_30, hi1_32, hk_27, hk_28, hk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gk_13[k]
                  + f_10 * hi0_12[k]
                  - f_11 * hi1_29[k]
                  + pb_x[k] * hk_27[k];

        t_14[k] = f_5 * gk_14[k]
                  + f_12 * hi0_13[k]
                  - f_13 * hi1_30[k]
                  + pb_x[k] * hk_28[k];

        t_15[k] = f_5 * gk_15[k]
                  + f_14 * hi0_14[k]
                  - f_15 * hi1_32[k]
                  + pb_x[k] * hk_29[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_x, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gk_17, gl_3, gl_16, hi0_16, hi1_34, hk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_16 * fl0_4[k]
                  - f_17 * fl1_4[k]
                  + pa_x[k] * gl_16[k];

        t_17[k] = f_16 * fl0_1[k]
                  - f_17 * fl1_1[k]
                  + pa_y[k] * gl_3[k];

        t_18[k] = f_18 * gk_17[k]
                  + f_6 * hi0_16[k]
                  - f_7 * hi1_34[k]
                  + pb_x[k] * hk_32[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gk_18, gk_19, gk_20, hi0_17, hi0_18, hi0_19, \
                         hi1_35, hi1_36, hi1_37, hk_33, hk_34, hk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_18 * gk_18[k]
                  + f_8 * hi0_17[k]
                  - f_9 * hi1_35[k]
                  + pb_x[k] * hk_33[k];

        t_20[k] = f_18 * gk_19[k]
                  + f_10 * hi0_18[k]
                  - f_11 * hi1_36[k]
                  + pb_x[k] * hk_34[k];

        t_21[k] = f_18 * gk_20[k]
                  + f_12 * hi0_19[k]
                  - f_13 * hi1_37[k]
                  + pb_x[k] * hk_35[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, pb_x, fl0_5, fl1_5, gk_21, gl_4, \
                         gl_5, gl_17, hi0_20, hi1_38, hk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_18 * gk_21[k]
                  + f_14 * hi0_20[k]
                  - f_15 * hi1_38[k]
                  + pb_x[k] * hk_36[k];

        t_23[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_x[k] * gl_17[k];

        t_24[k] = pa_z[k] * gl_4[k];

        t_25[k] = pa_z[k] * gl_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, gl_6, gl_7, \
                         gl_8, gl_10, gl_11, gl_12, gl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * gl_6[k];

        t_27[k] = pa_z[k] * gl_7[k];

        t_28[k] = pa_z[k] * gl_8[k];

        t_29[k] = pa_y[k] * gl_10[k];

        t_30[k] = pa_y[k] * gl_11[k];

        t_31[k] = pa_y[k] * gl_12[k];

        t_32[k] = pa_y[k] * gl_13[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pa_z, pb_x, fl0_2, fl1_2, gk_23, gl_10, \
                         gl_14, gl_15, hi0_23, hi1_42, hk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * gl_14[k];

        t_34[k] = pa_y[k] * gl_15[k];

        t_35[k] = f_16 * fl0_2[k]
                  - f_17 * fl1_2[k]
                  + pa_z[k] * gl_10[k];

        t_36[k] = f_18 * gk_23[k]
                  + f_6 * hi0_23[k]
                  - f_7 * hi1_42[k]
                  + pb_x[k] * hk_48[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, gk_24, gk_25, gk_26, hi0_24, hi0_25, hi0_26, \
                         hi1_43, hi1_44, hi1_45, hk_49, hk_50, hk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_18 * gk_24[k]
                  + f_8 * hi0_24[k]
                  - f_9 * hi1_43[k]
                  + pb_x[k] * hk_49[k];

        t_38[k] = f_18 * gk_25[k]
                  + f_10 * hi0_25[k]
                  - f_11 * hi1_44[k]
                  + pb_x[k] * hk_50[k];

        t_39[k] = f_18 * gk_26[k]
                  + f_12 * hi0_26[k]
                  - f_13 * hi1_45[k]
                  + pb_x[k] * hk_51[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, fl0_8, fl1_8, gk_27, gl_18, \
                         gl_19, gl_21, hi0_27, hi1_47, hk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_18 * gk_27[k]
                  + f_14 * hi0_27[k]
                  - f_15 * hi1_47[k]
                  + pb_x[k] * hk_52[k];

        t_41[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_x[k] * gl_18[k];

        t_42[k] = pa_x[k] * gl_19[k];

        t_43[k] = pa_x[k] * gl_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, pa_x, gl_22, gl_23, gl_24, \
                         gl_25, gl_26, gl_27, gl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gl_22[k];

        t_45[k] = pa_x[k] * gl_23[k];

        t_46[k] = pa_x[k] * gl_24[k];

        t_47[k] = pa_x[k] * gl_25[k];

        t_48[k] = pa_x[k] * gl_26[k];

        t_49[k] = pa_x[k] * gl_27[k];

        t_50[k] = pa_x[k] * gl_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, fl0_5, fl1_5, gk_33, gl_19, gl_20, \
                         hi0_43, hi1_74, hk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * gk_33[k]
                  + f_1 * hi0_43[k]
                  - f_2 * hi1_74[k]
                  + pb_y[k] * hk_81[k];

        t_52[k] = pa_z[k] * gl_19[k];

        t_53[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_z[k] * gl_20[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, gk_36, gk_37, gk_38, hi0_46, hi0_47, hi0_48, \
                         hi1_85, hi1_86, hi1_87, hk_89, hk_90, hk_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * gk_36[k]
                  + f_6 * hi0_46[k]
                  - f_7 * hi1_85[k]
                  + pb_y[k] * hk_89[k];

        t_55[k] = f_5 * gk_37[k]
                  + f_8 * hi0_47[k]
                  - f_9 * hi1_86[k]
                  + pb_y[k] * hk_90[k];

        t_56[k] = f_5 * gk_38[k]
                  + f_10 * hi0_48[k]
                  - f_11 * hi1_87[k]
                  + pb_y[k] * hk_91[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pb_y, fl0_7, fl1_7, gk_39, gk_40, gl_27, \
                         hi0_49, hi0_50, hi1_88, hi1_89, hk_92, hk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * gk_39[k]
                  + f_12 * hi0_49[k]
                  - f_13 * hi1_88[k]
                  + pb_y[k] * hk_92[k];

        t_58[k] = f_5 * gk_40[k]
                  + f_14 * hi0_50[k]
                  - f_15 * hi1_89[k]
                  + pb_y[k] * hk_93[k];

        t_59[k] = f_16 * fl0_7[k]
                  - f_17 * fl1_7[k]
                  + pa_y[k] * gl_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_z, pb_y, fl0_6, fl1_6, gk_42, gk_43, gl_21, \
                         hi0_52, hi0_53, hi1_93, hi1_94, hk_96, hk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_16 * fl0_6[k]
                  - f_17 * fl1_6[k]
                  + pa_z[k] * gl_21[k];

        t_61[k] = f_18 * gk_42[k]
                  + f_6 * hi0_52[k]
                  - f_7 * hi1_93[k]
                  + pb_y[k] * hk_96[k];

        t_62[k] = f_18 * gk_43[k]
                  + f_8 * hi0_53[k]
                  - f_9 * hi1_94[k]
                  + pb_y[k] * hk_97[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, gk_44, gk_45, gk_46, hi0_54, hi0_55, hi0_56, \
                         hi1_95, hi1_96, hi1_97, hk_98, hk_99, hk_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_18 * gk_44[k]
                  + f_10 * hi0_54[k]
                  - f_11 * hi1_95[k]
                  + pb_y[k] * hk_98[k];

        t_64[k] = f_18 * gk_45[k]
                  + f_12 * hi0_55[k]
                  - f_13 * hi1_96[k]
                  + pb_y[k] * hk_99[k];

        t_65[k] = f_18 * gk_46[k]
                  + f_14 * hi0_56[k]
                  - f_15 * hi1_97[k]
                  + pb_y[k] * hk_100[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pb_z, fl0_8, fl1_8, gk_56, gl_28, gl_29, \
                         hi0_62, hi1_118, hk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_y[k] * gl_28[k];

        t_67[k] = pa_y[k] * gl_29[k];

        t_68[k] = f_0 * gk_56[k]
                  + f_1 * hi0_62[k]
                  - f_2 * hi1_118[k]
                  + pb_z[k] * hk_122[k];
    }
}

auto
compute_prim_hl_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fl0, const size_t fl1,
                                     const size_t gk, const size_t gl, const size_t hi0,
                                     const size_t hi1, const size_t hk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
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
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);
    const auto f_18 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fl0_0 = buffer.data(fl0 + 0);
    const auto *fl0_1 = buffer.data(fl0 + 1);
    const auto *fl0_2 = buffer.data(fl0 + 2);
    const auto *fl0_3 = buffer.data(fl0 + 3);
    const auto *fl0_4 = buffer.data(fl0 + 4);
    const auto *fl0_5 = buffer.data(fl0 + 5);
    const auto *fl0_6 = buffer.data(fl0 + 6);
    const auto *fl0_7 = buffer.data(fl0 + 7);
    const auto *fl0_8 = buffer.data(fl0 + 8);

    const auto *fl1_0 = buffer.data(fl1 + 0);
    const auto *fl1_1 = buffer.data(fl1 + 1);
    const auto *fl1_2 = buffer.data(fl1 + 2);
    const auto *fl1_3 = buffer.data(fl1 + 3);
    const auto *fl1_4 = buffer.data(fl1 + 4);
    const auto *fl1_5 = buffer.data(fl1 + 5);
    const auto *fl1_6 = buffer.data(fl1 + 6);
    const auto *fl1_7 = buffer.data(fl1 + 7);
    const auto *fl1_8 = buffer.data(fl1 + 8);

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
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_56 = buffer.data(gk + 56);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);

    const auto *hi0_0 = buffer.data(hi0 + 0);
    const auto *hi0_20 = buffer.data(hi0 + 20);
    const auto *hi0_21 = buffer.data(hi0 + 21);
    const auto *hi0_22 = buffer.data(hi0 + 22);
    const auto *hi0_23 = buffer.data(hi0 + 23);
    const auto *hi0_24 = buffer.data(hi0 + 24);
    const auto *hi0_27 = buffer.data(hi0 + 27);
    const auto *hi0_28 = buffer.data(hi0 + 28);
    const auto *hi0_29 = buffer.data(hi0 + 29);
    const auto *hi0_30 = buffer.data(hi0 + 30);
    const auto *hi0_32 = buffer.data(hi0 + 32);
    const auto *hi0_34 = buffer.data(hi0 + 34);
    const auto *hi0_35 = buffer.data(hi0 + 35);
    const auto *hi0_36 = buffer.data(hi0 + 36);
    const auto *hi0_37 = buffer.data(hi0 + 37);
    const auto *hi0_38 = buffer.data(hi0 + 38);
    const auto *hi0_42 = buffer.data(hi0 + 42);
    const auto *hi0_43 = buffer.data(hi0 + 43);
    const auto *hi0_44 = buffer.data(hi0 + 44);
    const auto *hi0_45 = buffer.data(hi0 + 45);
    const auto *hi0_47 = buffer.data(hi0 + 47);
    const auto *hi0_74 = buffer.data(hi0 + 74);
    const auto *hi0_85 = buffer.data(hi0 + 85);
    const auto *hi0_86 = buffer.data(hi0 + 86);
    const auto *hi0_87 = buffer.data(hi0 + 87);
    const auto *hi0_88 = buffer.data(hi0 + 88);
    const auto *hi0_89 = buffer.data(hi0 + 89);
    const auto *hi0_93 = buffer.data(hi0 + 93);
    const auto *hi0_94 = buffer.data(hi0 + 94);
    const auto *hi0_95 = buffer.data(hi0 + 95);
    const auto *hi0_96 = buffer.data(hi0 + 96);
    const auto *hi0_97 = buffer.data(hi0 + 97);
    const auto *hi0_118 = buffer.data(hi0 + 118);

    const auto *hi1_0 = buffer.data(hi1 + 0);
    const auto *hi1_18 = buffer.data(hi1 + 18);
    const auto *hi1_19 = buffer.data(hi1 + 19);
    const auto *hi1_20 = buffer.data(hi1 + 20);
    const auto *hi1_21 = buffer.data(hi1 + 21);
    const auto *hi1_22 = buffer.data(hi1 + 22);
    const auto *hi1_24 = buffer.data(hi1 + 24);
    const auto *hi1_25 = buffer.data(hi1 + 25);
    const auto *hi1_26 = buffer.data(hi1 + 26);
    const auto *hi1_27 = buffer.data(hi1 + 27);
    const auto *hi1_28 = buffer.data(hi1 + 28);
    const auto *hi1_30 = buffer.data(hi1 + 30);
    const auto *hi1_31 = buffer.data(hi1 + 31);
    const auto *hi1_32 = buffer.data(hi1 + 32);
    const auto *hi1_33 = buffer.data(hi1 + 33);
    const auto *hi1_34 = buffer.data(hi1 + 34);
    const auto *hi1_37 = buffer.data(hi1 + 37);
    const auto *hi1_38 = buffer.data(hi1 + 38);
    const auto *hi1_39 = buffer.data(hi1 + 39);
    const auto *hi1_40 = buffer.data(hi1 + 40);
    const auto *hi1_41 = buffer.data(hi1 + 41);
    const auto *hi1_66 = buffer.data(hi1 + 66);
    const auto *hi1_74 = buffer.data(hi1 + 74);
    const auto *hi1_75 = buffer.data(hi1 + 75);
    const auto *hi1_76 = buffer.data(hi1 + 76);
    const auto *hi1_77 = buffer.data(hi1 + 77);
    const auto *hi1_78 = buffer.data(hi1 + 78);
    const auto *hi1_80 = buffer.data(hi1 + 80);
    const auto *hi1_81 = buffer.data(hi1 + 81);
    const auto *hi1_82 = buffer.data(hi1 + 82);
    const auto *hi1_83 = buffer.data(hi1 + 83);
    const auto *hi1_84 = buffer.data(hi1 + 84);
    const auto *hi1_104 = buffer.data(hi1 + 104);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_65 = buffer.data(hk + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, fl0_0, fl1_0, gk_0, gl_0, gl_1, \
                         hi0_0, hi1_0, hk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_0[k]
                 + f_1 * hi0_0[k]
                 - f_2 * hi1_0[k]
                 + pb_x[k] * hk_0[k];

        t_1[k] = pa_y[k] * gl_0[k];

        t_2[k] = pa_z[k] * gl_0[k];

        t_3[k] = f_3 * fl0_0[k]
                 - f_4 * fl1_0[k]
                 + pa_y[k] * gl_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, gk_4, gk_5, gk_6, hi0_20, hi0_21, hi0_22, \
                         hi1_18, hi1_19, hi1_20, hk_4, hk_5, hk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gk_4[k]
                 + f_6 * hi0_20[k]
                 - f_7 * hi1_18[k]
                 + pb_x[k] * hk_4[k];

        t_5[k] = f_5 * gk_5[k]
                 + f_8 * hi0_21[k]
                 - f_9 * hi1_19[k]
                 + pb_x[k] * hk_5[k];

        t_6[k] = f_5 * gk_6[k]
                 + f_10 * hi0_22[k]
                 - f_11 * hi1_20[k]
                 + pb_x[k] * hk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, fl0_3, fl1_3, gk_7, gk_8, gl_9, hi0_23, \
                         hi0_24, hi1_21, hi1_22, hk_7, hk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * gk_7[k]
                 + f_12 * hi0_23[k]
                 - f_13 * hi1_21[k]
                 + pb_x[k] * hk_7[k];

        t_8[k] = f_5 * gk_8[k]
                 + f_14 * hi0_24[k]
                 - f_15 * hi1_22[k]
                 + pb_x[k] * hk_8[k];

        t_9[k] = f_16 * fl0_3[k]
                 - f_17 * fl1_3[k]
                 + pa_x[k] * gl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_x, fl0_0, fl1_0, gk_11, gk_12, gl_2, \
                         hi0_27, hi0_28, hi1_24, hi1_25, hk_11, hk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fl0_0[k]
                  - f_4 * fl1_0[k]
                  + pa_z[k] * gl_2[k];

        t_11[k] = f_5 * gk_11[k]
                  + f_6 * hi0_27[k]
                  - f_7 * hi1_24[k]
                  + pb_x[k] * hk_11[k];

        t_12[k] = f_5 * gk_12[k]
                  + f_8 * hi0_28[k]
                  - f_9 * hi1_25[k]
                  + pb_x[k] * hk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, gk_13, gk_14, gk_15, hi0_29, hi0_30, hi0_32, \
                         hi1_26, hi1_27, hi1_28, hk_13, hk_14, hk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gk_13[k]
                  + f_10 * hi0_29[k]
                  - f_11 * hi1_26[k]
                  + pb_x[k] * hk_13[k];

        t_14[k] = f_5 * gk_14[k]
                  + f_12 * hi0_30[k]
                  - f_13 * hi1_27[k]
                  + pb_x[k] * hk_14[k];

        t_15[k] = f_5 * gk_15[k]
                  + f_14 * hi0_32[k]
                  - f_15 * hi1_28[k]
                  + pb_x[k] * hk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_x, fl0_1, fl0_4, fl1_1, fl1_4, \
                         gk_17, gl_3, gl_16, hi0_34, hi1_30, hk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_16 * fl0_4[k]
                  - f_17 * fl1_4[k]
                  + pa_x[k] * gl_16[k];

        t_17[k] = f_16 * fl0_1[k]
                  - f_17 * fl1_1[k]
                  + pa_y[k] * gl_3[k];

        t_18[k] = f_18 * gk_17[k]
                  + f_6 * hi0_34[k]
                  - f_7 * hi1_30[k]
                  + pb_x[k] * hk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gk_18, gk_19, gk_20, hi0_35, hi0_36, hi0_37, \
                         hi1_31, hi1_32, hi1_33, hk_19, hk_20, hk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_18 * gk_18[k]
                  + f_8 * hi0_35[k]
                  - f_9 * hi1_31[k]
                  + pb_x[k] * hk_19[k];

        t_20[k] = f_18 * gk_19[k]
                  + f_10 * hi0_36[k]
                  - f_11 * hi1_32[k]
                  + pb_x[k] * hk_20[k];

        t_21[k] = f_18 * gk_20[k]
                  + f_12 * hi0_37[k]
                  - f_13 * hi1_33[k]
                  + pb_x[k] * hk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, pb_x, fl0_5, fl1_5, gk_21, gl_4, \
                         gl_5, gl_17, hi0_38, hi1_34, hk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_18 * gk_21[k]
                  + f_14 * hi0_38[k]
                  - f_15 * hi1_34[k]
                  + pb_x[k] * hk_22[k];

        t_23[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_x[k] * gl_17[k];

        t_24[k] = pa_z[k] * gl_4[k];

        t_25[k] = pa_z[k] * gl_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, gl_6, gl_7, \
                         gl_8, gl_10, gl_11, gl_12, gl_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * gl_6[k];

        t_27[k] = pa_z[k] * gl_7[k];

        t_28[k] = pa_z[k] * gl_8[k];

        t_29[k] = pa_y[k] * gl_10[k];

        t_30[k] = pa_y[k] * gl_11[k];

        t_31[k] = pa_y[k] * gl_12[k];

        t_32[k] = pa_y[k] * gl_13[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pa_z, pb_x, fl0_2, fl1_2, gk_23, gl_10, \
                         gl_14, gl_15, hi0_42, hi1_37, hk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * gl_14[k];

        t_34[k] = pa_y[k] * gl_15[k];

        t_35[k] = f_16 * fl0_2[k]
                  - f_17 * fl1_2[k]
                  + pa_z[k] * gl_10[k];

        t_36[k] = f_18 * gk_23[k]
                  + f_6 * hi0_42[k]
                  - f_7 * hi1_37[k]
                  + pb_x[k] * hk_34[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, gk_24, gk_25, gk_26, hi0_43, hi0_44, hi0_45, \
                         hi1_38, hi1_39, hi1_40, hk_35, hk_36, hk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_18 * gk_24[k]
                  + f_8 * hi0_43[k]
                  - f_9 * hi1_38[k]
                  + pb_x[k] * hk_35[k];

        t_38[k] = f_18 * gk_25[k]
                  + f_10 * hi0_44[k]
                  - f_11 * hi1_39[k]
                  + pb_x[k] * hk_36[k];

        t_39[k] = f_18 * gk_26[k]
                  + f_12 * hi0_45[k]
                  - f_13 * hi1_40[k]
                  + pb_x[k] * hk_37[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, fl0_8, fl1_8, gk_27, gl_18, \
                         gl_19, gl_21, hi0_47, hi1_41, hk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_18 * gk_27[k]
                  + f_14 * hi0_47[k]
                  - f_15 * hi1_41[k]
                  + pb_x[k] * hk_38[k];

        t_41[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_x[k] * gl_18[k];

        t_42[k] = pa_x[k] * gl_19[k];

        t_43[k] = pa_x[k] * gl_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, pa_x, gl_22, gl_23, gl_24, \
                         gl_25, gl_26, gl_27, gl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gl_22[k];

        t_45[k] = pa_x[k] * gl_23[k];

        t_46[k] = pa_x[k] * gl_24[k];

        t_47[k] = pa_x[k] * gl_25[k];

        t_48[k] = pa_x[k] * gl_26[k];

        t_49[k] = pa_x[k] * gl_27[k];

        t_50[k] = pa_x[k] * gl_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, fl0_5, fl1_5, gk_33, gl_19, gl_20, \
                         hi0_74, hi1_66, hk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * gk_33[k]
                  + f_1 * hi0_74[k]
                  - f_2 * hi1_66[k]
                  + pb_y[k] * hk_48[k];

        t_52[k] = pa_z[k] * gl_19[k];

        t_53[k] = f_3 * fl0_5[k]
                  - f_4 * fl1_5[k]
                  + pa_z[k] * gl_20[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, gk_36, gk_37, gk_38, hi0_85, hi0_86, hi0_87, \
                         hi1_74, hi1_75, hi1_76, hk_51, hk_52, hk_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * gk_36[k]
                  + f_6 * hi0_85[k]
                  - f_7 * hi1_74[k]
                  + pb_y[k] * hk_51[k];

        t_55[k] = f_5 * gk_37[k]
                  + f_8 * hi0_86[k]
                  - f_9 * hi1_75[k]
                  + pb_y[k] * hk_52[k];

        t_56[k] = f_5 * gk_38[k]
                  + f_10 * hi0_87[k]
                  - f_11 * hi1_76[k]
                  + pb_y[k] * hk_53[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pb_y, fl0_7, fl1_7, gk_39, gk_40, gl_27, \
                         hi0_88, hi0_89, hi1_77, hi1_78, hk_54, hk_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * gk_39[k]
                  + f_12 * hi0_88[k]
                  - f_13 * hi1_77[k]
                  + pb_y[k] * hk_54[k];

        t_58[k] = f_5 * gk_40[k]
                  + f_14 * hi0_89[k]
                  - f_15 * hi1_78[k]
                  + pb_y[k] * hk_55[k];

        t_59[k] = f_16 * fl0_7[k]
                  - f_17 * fl1_7[k]
                  + pa_y[k] * gl_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_z, pb_y, fl0_6, fl1_6, gk_42, gk_43, gl_21, \
                         hi0_93, hi0_94, hi1_80, hi1_81, hk_58, hk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_16 * fl0_6[k]
                  - f_17 * fl1_6[k]
                  + pa_z[k] * gl_21[k];

        t_61[k] = f_18 * gk_42[k]
                  + f_6 * hi0_93[k]
                  - f_7 * hi1_80[k]
                  + pb_y[k] * hk_58[k];

        t_62[k] = f_18 * gk_43[k]
                  + f_8 * hi0_94[k]
                  - f_9 * hi1_81[k]
                  + pb_y[k] * hk_59[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, gk_44, gk_45, gk_46, hi0_95, hi0_96, hi0_97, \
                         hi1_82, hi1_83, hi1_84, hk_60, hk_61, hk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_18 * gk_44[k]
                  + f_10 * hi0_95[k]
                  - f_11 * hi1_82[k]
                  + pb_y[k] * hk_60[k];

        t_64[k] = f_18 * gk_45[k]
                  + f_12 * hi0_96[k]
                  - f_13 * hi1_83[k]
                  + pb_y[k] * hk_61[k];

        t_65[k] = f_18 * gk_46[k]
                  + f_14 * hi0_97[k]
                  - f_15 * hi1_84[k]
                  + pb_y[k] * hk_62[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pb_z, fl0_8, fl1_8, gk_56, gl_28, gl_29, \
                         hi0_118, hi1_104, hk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * fl0_8[k]
                  - f_4 * fl1_8[k]
                  + pa_y[k] * gl_28[k];

        t_67[k] = pa_y[k] * gl_29[k];

        t_68[k] = f_0 * gk_56[k]
                  + f_1 * hi0_118[k]
                  - f_2 * hi1_104[k]
                  + pb_z[k] * hk_65[k];
    }
}

}  // namespace simdt2ceri
