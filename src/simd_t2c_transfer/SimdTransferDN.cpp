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


#include "SimdTransferDN.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dn(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pn, const size_t po, const size_t nmax) -> void
{
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *pn_0 = buffer.data(pn + 0);
    const auto *pn_1 = buffer.data(pn + 1);
    const auto *pn_2 = buffer.data(pn + 2);
    const auto *pn_3 = buffer.data(pn + 3);
    const auto *pn_4 = buffer.data(pn + 4);
    const auto *pn_5 = buffer.data(pn + 5);
    const auto *pn_6 = buffer.data(pn + 6);
    const auto *pn_7 = buffer.data(pn + 7);
    const auto *pn_8 = buffer.data(pn + 8);
    const auto *pn_9 = buffer.data(pn + 9);
    const auto *pn_10 = buffer.data(pn + 10);
    const auto *pn_11 = buffer.data(pn + 11);
    const auto *pn_12 = buffer.data(pn + 12);
    const auto *pn_13 = buffer.data(pn + 13);
    const auto *pn_14 = buffer.data(pn + 14);
    const auto *pn_15 = buffer.data(pn + 15);
    const auto *pn_16 = buffer.data(pn + 16);
    const auto *pn_17 = buffer.data(pn + 17);
    const auto *pn_18 = buffer.data(pn + 18);
    const auto *pn_19 = buffer.data(pn + 19);
    const auto *pn_20 = buffer.data(pn + 20);
    const auto *pn_21 = buffer.data(pn + 21);
    const auto *pn_22 = buffer.data(pn + 22);
    const auto *pn_23 = buffer.data(pn + 23);
    const auto *pn_24 = buffer.data(pn + 24);
    const auto *pn_25 = buffer.data(pn + 25);
    const auto *pn_26 = buffer.data(pn + 26);
    const auto *pn_27 = buffer.data(pn + 27);
    const auto *pn_28 = buffer.data(pn + 28);
    const auto *pn_29 = buffer.data(pn + 29);
    const auto *pn_30 = buffer.data(pn + 30);
    const auto *pn_31 = buffer.data(pn + 31);
    const auto *pn_32 = buffer.data(pn + 32);
    const auto *pn_33 = buffer.data(pn + 33);
    const auto *pn_34 = buffer.data(pn + 34);
    const auto *pn_35 = buffer.data(pn + 35);
    const auto *pn_36 = buffer.data(pn + 36);
    const auto *pn_37 = buffer.data(pn + 37);
    const auto *pn_38 = buffer.data(pn + 38);
    const auto *pn_39 = buffer.data(pn + 39);
    const auto *pn_40 = buffer.data(pn + 40);
    const auto *pn_41 = buffer.data(pn + 41);
    const auto *pn_42 = buffer.data(pn + 42);
    const auto *pn_43 = buffer.data(pn + 43);
    const auto *pn_44 = buffer.data(pn + 44);
    const auto *pn_45 = buffer.data(pn + 45);
    const auto *pn_46 = buffer.data(pn + 46);
    const auto *pn_47 = buffer.data(pn + 47);
    const auto *pn_48 = buffer.data(pn + 48);
    const auto *pn_49 = buffer.data(pn + 49);
    const auto *pn_50 = buffer.data(pn + 50);
    const auto *pn_51 = buffer.data(pn + 51);
    const auto *pn_52 = buffer.data(pn + 52);
    const auto *pn_53 = buffer.data(pn + 53);
    const auto *pn_54 = buffer.data(pn + 54);
    const auto *pn_55 = buffer.data(pn + 55);
    const auto *pn_56 = buffer.data(pn + 56);
    const auto *pn_57 = buffer.data(pn + 57);
    const auto *pn_58 = buffer.data(pn + 58);
    const auto *pn_59 = buffer.data(pn + 59);
    const auto *pn_60 = buffer.data(pn + 60);
    const auto *pn_61 = buffer.data(pn + 61);
    const auto *pn_62 = buffer.data(pn + 62);
    const auto *pn_63 = buffer.data(pn + 63);
    const auto *pn_64 = buffer.data(pn + 64);
    const auto *pn_65 = buffer.data(pn + 65);
    const auto *pn_66 = buffer.data(pn + 66);
    const auto *pn_67 = buffer.data(pn + 67);
    const auto *pn_68 = buffer.data(pn + 68);
    const auto *pn_69 = buffer.data(pn + 69);
    const auto *pn_70 = buffer.data(pn + 70);
    const auto *pn_71 = buffer.data(pn + 71);
    const auto *pn_72 = buffer.data(pn + 72);
    const auto *pn_73 = buffer.data(pn + 73);
    const auto *pn_74 = buffer.data(pn + 74);
    const auto *pn_75 = buffer.data(pn + 75);
    const auto *pn_76 = buffer.data(pn + 76);
    const auto *pn_77 = buffer.data(pn + 77);
    const auto *pn_78 = buffer.data(pn + 78);
    const auto *pn_79 = buffer.data(pn + 79);
    const auto *pn_80 = buffer.data(pn + 80);
    const auto *pn_81 = buffer.data(pn + 81);
    const auto *pn_82 = buffer.data(pn + 82);
    const auto *pn_83 = buffer.data(pn + 83);
    const auto *pn_84 = buffer.data(pn + 84);
    const auto *pn_85 = buffer.data(pn + 85);
    const auto *pn_86 = buffer.data(pn + 86);
    const auto *pn_87 = buffer.data(pn + 87);
    const auto *pn_88 = buffer.data(pn + 88);
    const auto *pn_89 = buffer.data(pn + 89);
    const auto *pn_90 = buffer.data(pn + 90);
    const auto *pn_91 = buffer.data(pn + 91);
    const auto *pn_92 = buffer.data(pn + 92);
    const auto *pn_93 = buffer.data(pn + 93);
    const auto *pn_94 = buffer.data(pn + 94);
    const auto *pn_95 = buffer.data(pn + 95);
    const auto *pn_96 = buffer.data(pn + 96);
    const auto *pn_97 = buffer.data(pn + 97);
    const auto *pn_98 = buffer.data(pn + 98);
    const auto *pn_99 = buffer.data(pn + 99);
    const auto *pn_100 = buffer.data(pn + 100);
    const auto *pn_101 = buffer.data(pn + 101);
    const auto *pn_102 = buffer.data(pn + 102);
    const auto *pn_103 = buffer.data(pn + 103);
    const auto *pn_104 = buffer.data(pn + 104);
    const auto *pn_105 = buffer.data(pn + 105);
    const auto *pn_106 = buffer.data(pn + 106);
    const auto *pn_107 = buffer.data(pn + 107);
    const auto *pn_108 = buffer.data(pn + 108);
    const auto *pn_109 = buffer.data(pn + 109);
    const auto *pn_110 = buffer.data(pn + 110);
    const auto *pn_111 = buffer.data(pn + 111);
    const auto *pn_112 = buffer.data(pn + 112);
    const auto *pn_113 = buffer.data(pn + 113);
    const auto *pn_114 = buffer.data(pn + 114);
    const auto *pn_115 = buffer.data(pn + 115);
    const auto *pn_116 = buffer.data(pn + 116);
    const auto *pn_117 = buffer.data(pn + 117);
    const auto *pn_118 = buffer.data(pn + 118);
    const auto *pn_119 = buffer.data(pn + 119);
    const auto *pn_120 = buffer.data(pn + 120);
    const auto *pn_121 = buffer.data(pn + 121);
    const auto *pn_122 = buffer.data(pn + 122);
    const auto *pn_123 = buffer.data(pn + 123);
    const auto *pn_124 = buffer.data(pn + 124);
    const auto *pn_125 = buffer.data(pn + 125);
    const auto *pn_126 = buffer.data(pn + 126);
    const auto *pn_127 = buffer.data(pn + 127);
    const auto *pn_128 = buffer.data(pn + 128);
    const auto *pn_129 = buffer.data(pn + 129);
    const auto *pn_130 = buffer.data(pn + 130);
    const auto *pn_131 = buffer.data(pn + 131);
    const auto *pn_132 = buffer.data(pn + 132);
    const auto *pn_133 = buffer.data(pn + 133);
    const auto *pn_134 = buffer.data(pn + 134);
    const auto *pn_135 = buffer.data(pn + 135);
    const auto *pn_136 = buffer.data(pn + 136);
    const auto *pn_137 = buffer.data(pn + 137);
    const auto *pn_138 = buffer.data(pn + 138);
    const auto *pn_139 = buffer.data(pn + 139);
    const auto *pn_140 = buffer.data(pn + 140);
    const auto *pn_141 = buffer.data(pn + 141);
    const auto *pn_142 = buffer.data(pn + 142);
    const auto *pn_143 = buffer.data(pn + 143);
    const auto *pn_144 = buffer.data(pn + 144);
    const auto *pn_145 = buffer.data(pn + 145);
    const auto *pn_146 = buffer.data(pn + 146);
    const auto *pn_147 = buffer.data(pn + 147);
    const auto *pn_148 = buffer.data(pn + 148);
    const auto *pn_149 = buffer.data(pn + 149);
    const auto *pn_150 = buffer.data(pn + 150);
    const auto *pn_151 = buffer.data(pn + 151);
    const auto *pn_152 = buffer.data(pn + 152);
    const auto *pn_153 = buffer.data(pn + 153);
    const auto *pn_154 = buffer.data(pn + 154);
    const auto *pn_155 = buffer.data(pn + 155);
    const auto *pn_156 = buffer.data(pn + 156);
    const auto *pn_157 = buffer.data(pn + 157);
    const auto *pn_158 = buffer.data(pn + 158);
    const auto *pn_159 = buffer.data(pn + 159);
    const auto *pn_160 = buffer.data(pn + 160);
    const auto *pn_161 = buffer.data(pn + 161);
    const auto *pn_162 = buffer.data(pn + 162);
    const auto *pn_163 = buffer.data(pn + 163);
    const auto *pn_164 = buffer.data(pn + 164);
    const auto *pn_165 = buffer.data(pn + 165);
    const auto *pn_166 = buffer.data(pn + 166);
    const auto *pn_167 = buffer.data(pn + 167);
    const auto *pn_168 = buffer.data(pn + 168);
    const auto *pn_169 = buffer.data(pn + 169);
    const auto *pn_170 = buffer.data(pn + 170);
    const auto *pn_171 = buffer.data(pn + 171);
    const auto *pn_172 = buffer.data(pn + 172);
    const auto *pn_173 = buffer.data(pn + 173);
    const auto *pn_174 = buffer.data(pn + 174);
    const auto *pn_175 = buffer.data(pn + 175);
    const auto *pn_176 = buffer.data(pn + 176);
    const auto *pn_177 = buffer.data(pn + 177);
    const auto *pn_178 = buffer.data(pn + 178);
    const auto *pn_179 = buffer.data(pn + 179);
    const auto *pn_180 = buffer.data(pn + 180);
    const auto *pn_181 = buffer.data(pn + 181);
    const auto *pn_182 = buffer.data(pn + 182);
    const auto *pn_183 = buffer.data(pn + 183);
    const auto *pn_184 = buffer.data(pn + 184);
    const auto *pn_185 = buffer.data(pn + 185);
    const auto *pn_186 = buffer.data(pn + 186);
    const auto *pn_187 = buffer.data(pn + 187);
    const auto *pn_188 = buffer.data(pn + 188);
    const auto *pn_189 = buffer.data(pn + 189);
    const auto *pn_190 = buffer.data(pn + 190);
    const auto *pn_191 = buffer.data(pn + 191);
    const auto *pn_192 = buffer.data(pn + 192);
    const auto *pn_193 = buffer.data(pn + 193);
    const auto *pn_194 = buffer.data(pn + 194);
    const auto *pn_195 = buffer.data(pn + 195);
    const auto *pn_196 = buffer.data(pn + 196);
    const auto *pn_197 = buffer.data(pn + 197);

    const auto *po_0 = buffer.data(po + 0);
    const auto *po_1 = buffer.data(po + 1);
    const auto *po_2 = buffer.data(po + 2);
    const auto *po_3 = buffer.data(po + 3);
    const auto *po_4 = buffer.data(po + 4);
    const auto *po_5 = buffer.data(po + 5);
    const auto *po_6 = buffer.data(po + 6);
    const auto *po_7 = buffer.data(po + 7);
    const auto *po_8 = buffer.data(po + 8);
    const auto *po_9 = buffer.data(po + 9);
    const auto *po_10 = buffer.data(po + 10);
    const auto *po_11 = buffer.data(po + 11);
    const auto *po_12 = buffer.data(po + 12);
    const auto *po_13 = buffer.data(po + 13);
    const auto *po_14 = buffer.data(po + 14);
    const auto *po_15 = buffer.data(po + 15);
    const auto *po_16 = buffer.data(po + 16);
    const auto *po_17 = buffer.data(po + 17);
    const auto *po_18 = buffer.data(po + 18);
    const auto *po_19 = buffer.data(po + 19);
    const auto *po_20 = buffer.data(po + 20);
    const auto *po_21 = buffer.data(po + 21);
    const auto *po_22 = buffer.data(po + 22);
    const auto *po_23 = buffer.data(po + 23);
    const auto *po_24 = buffer.data(po + 24);
    const auto *po_25 = buffer.data(po + 25);
    const auto *po_26 = buffer.data(po + 26);
    const auto *po_27 = buffer.data(po + 27);
    const auto *po_28 = buffer.data(po + 28);
    const auto *po_29 = buffer.data(po + 29);
    const auto *po_30 = buffer.data(po + 30);
    const auto *po_31 = buffer.data(po + 31);
    const auto *po_32 = buffer.data(po + 32);
    const auto *po_33 = buffer.data(po + 33);
    const auto *po_34 = buffer.data(po + 34);
    const auto *po_35 = buffer.data(po + 35);
    const auto *po_36 = buffer.data(po + 36);
    const auto *po_37 = buffer.data(po + 37);
    const auto *po_38 = buffer.data(po + 38);
    const auto *po_39 = buffer.data(po + 39);
    const auto *po_40 = buffer.data(po + 40);
    const auto *po_41 = buffer.data(po + 41);
    const auto *po_42 = buffer.data(po + 42);
    const auto *po_43 = buffer.data(po + 43);
    const auto *po_44 = buffer.data(po + 44);
    const auto *po_45 = buffer.data(po + 45);
    const auto *po_46 = buffer.data(po + 46);
    const auto *po_47 = buffer.data(po + 47);
    const auto *po_48 = buffer.data(po + 48);
    const auto *po_49 = buffer.data(po + 49);
    const auto *po_50 = buffer.data(po + 50);
    const auto *po_51 = buffer.data(po + 51);
    const auto *po_52 = buffer.data(po + 52);
    const auto *po_53 = buffer.data(po + 53);
    const auto *po_54 = buffer.data(po + 54);
    const auto *po_55 = buffer.data(po + 55);
    const auto *po_56 = buffer.data(po + 56);
    const auto *po_57 = buffer.data(po + 57);
    const auto *po_58 = buffer.data(po + 58);
    const auto *po_59 = buffer.data(po + 59);
    const auto *po_60 = buffer.data(po + 60);
    const auto *po_61 = buffer.data(po + 61);
    const auto *po_62 = buffer.data(po + 62);
    const auto *po_63 = buffer.data(po + 63);
    const auto *po_64 = buffer.data(po + 64);
    const auto *po_65 = buffer.data(po + 65);
    const auto *po_78 = buffer.data(po + 78);
    const auto *po_79 = buffer.data(po + 79);
    const auto *po_80 = buffer.data(po + 80);
    const auto *po_81 = buffer.data(po + 81);
    const auto *po_82 = buffer.data(po + 82);
    const auto *po_83 = buffer.data(po + 83);
    const auto *po_84 = buffer.data(po + 84);
    const auto *po_85 = buffer.data(po + 85);
    const auto *po_86 = buffer.data(po + 86);
    const auto *po_87 = buffer.data(po + 87);
    const auto *po_88 = buffer.data(po + 88);
    const auto *po_89 = buffer.data(po + 89);
    const auto *po_90 = buffer.data(po + 90);
    const auto *po_91 = buffer.data(po + 91);
    const auto *po_92 = buffer.data(po + 92);
    const auto *po_93 = buffer.data(po + 93);
    const auto *po_94 = buffer.data(po + 94);
    const auto *po_95 = buffer.data(po + 95);
    const auto *po_96 = buffer.data(po + 96);
    const auto *po_97 = buffer.data(po + 97);
    const auto *po_98 = buffer.data(po + 98);
    const auto *po_99 = buffer.data(po + 99);
    const auto *po_100 = buffer.data(po + 100);
    const auto *po_101 = buffer.data(po + 101);
    const auto *po_102 = buffer.data(po + 102);
    const auto *po_103 = buffer.data(po + 103);
    const auto *po_104 = buffer.data(po + 104);
    const auto *po_105 = buffer.data(po + 105);
    const auto *po_106 = buffer.data(po + 106);
    const auto *po_107 = buffer.data(po + 107);
    const auto *po_108 = buffer.data(po + 108);
    const auto *po_109 = buffer.data(po + 109);
    const auto *po_110 = buffer.data(po + 110);
    const auto *po_111 = buffer.data(po + 111);
    const auto *po_112 = buffer.data(po + 112);
    const auto *po_113 = buffer.data(po + 113);
    const auto *po_114 = buffer.data(po + 114);
    const auto *po_115 = buffer.data(po + 115);
    const auto *po_116 = buffer.data(po + 116);
    const auto *po_117 = buffer.data(po + 117);
    const auto *po_118 = buffer.data(po + 118);
    const auto *po_119 = buffer.data(po + 119);
    const auto *po_120 = buffer.data(po + 120);
    const auto *po_121 = buffer.data(po + 121);
    const auto *po_122 = buffer.data(po + 122);
    const auto *po_123 = buffer.data(po + 123);
    const auto *po_124 = buffer.data(po + 124);
    const auto *po_125 = buffer.data(po + 125);
    const auto *po_126 = buffer.data(po + 126);
    const auto *po_127 = buffer.data(po + 127);
    const auto *po_128 = buffer.data(po + 128);
    const auto *po_129 = buffer.data(po + 129);
    const auto *po_130 = buffer.data(po + 130);
    const auto *po_131 = buffer.data(po + 131);
    const auto *po_132 = buffer.data(po + 132);
    const auto *po_133 = buffer.data(po + 133);
    const auto *po_134 = buffer.data(po + 134);
    const auto *po_135 = buffer.data(po + 135);
    const auto *po_136 = buffer.data(po + 136);
    const auto *po_137 = buffer.data(po + 137);
    const auto *po_138 = buffer.data(po + 138);
    const auto *po_139 = buffer.data(po + 139);
    const auto *po_140 = buffer.data(po + 140);
    const auto *po_141 = buffer.data(po + 141);
    const auto *po_142 = buffer.data(po + 142);
    const auto *po_143 = buffer.data(po + 143);
    const auto *po_144 = buffer.data(po + 144);
    const auto *po_145 = buffer.data(po + 145);
    const auto *po_146 = buffer.data(po + 146);
    const auto *po_147 = buffer.data(po + 147);
    const auto *po_148 = buffer.data(po + 148);
    const auto *po_149 = buffer.data(po + 149);
    const auto *po_150 = buffer.data(po + 150);
    const auto *po_151 = buffer.data(po + 151);
    const auto *po_152 = buffer.data(po + 152);
    const auto *po_153 = buffer.data(po + 153);
    const auto *po_154 = buffer.data(po + 154);
    const auto *po_156 = buffer.data(po + 156);
    const auto *po_157 = buffer.data(po + 157);
    const auto *po_158 = buffer.data(po + 158);
    const auto *po_159 = buffer.data(po + 159);
    const auto *po_160 = buffer.data(po + 160);
    const auto *po_161 = buffer.data(po + 161);
    const auto *po_162 = buffer.data(po + 162);
    const auto *po_163 = buffer.data(po + 163);
    const auto *po_164 = buffer.data(po + 164);
    const auto *po_165 = buffer.data(po + 165);
    const auto *po_166 = buffer.data(po + 166);
    const auto *po_167 = buffer.data(po + 167);
    const auto *po_168 = buffer.data(po + 168);
    const auto *po_169 = buffer.data(po + 169);
    const auto *po_170 = buffer.data(po + 170);
    const auto *po_171 = buffer.data(po + 171);
    const auto *po_172 = buffer.data(po + 172);
    const auto *po_173 = buffer.data(po + 173);
    const auto *po_174 = buffer.data(po + 174);
    const auto *po_175 = buffer.data(po + 175);
    const auto *po_176 = buffer.data(po + 176);
    const auto *po_177 = buffer.data(po + 177);
    const auto *po_178 = buffer.data(po + 178);
    const auto *po_179 = buffer.data(po + 179);
    const auto *po_180 = buffer.data(po + 180);
    const auto *po_181 = buffer.data(po + 181);
    const auto *po_182 = buffer.data(po + 182);
    const auto *po_183 = buffer.data(po + 183);
    const auto *po_184 = buffer.data(po + 184);
    const auto *po_185 = buffer.data(po + 185);
    const auto *po_186 = buffer.data(po + 186);
    const auto *po_187 = buffer.data(po + 187);
    const auto *po_188 = buffer.data(po + 188);
    const auto *po_189 = buffer.data(po + 189);
    const auto *po_190 = buffer.data(po + 190);
    const auto *po_191 = buffer.data(po + 191);
    const auto *po_192 = buffer.data(po + 192);
    const auto *po_193 = buffer.data(po + 193);
    const auto *po_194 = buffer.data(po + 194);
    const auto *po_195 = buffer.data(po + 195);
    const auto *po_196 = buffer.data(po + 196);
    const auto *po_197 = buffer.data(po + 197);
    const auto *po_198 = buffer.data(po + 198);
    const auto *po_199 = buffer.data(po + 199);
    const auto *po_200 = buffer.data(po + 200);
    const auto *po_201 = buffer.data(po + 201);
    const auto *po_202 = buffer.data(po + 202);
    const auto *po_203 = buffer.data(po + 203);
    const auto *po_204 = buffer.data(po + 204);
    const auto *po_205 = buffer.data(po + 205);
    const auto *po_206 = buffer.data(po + 206);
    const auto *po_207 = buffer.data(po + 207);
    const auto *po_208 = buffer.data(po + 208);
    const auto *po_209 = buffer.data(po + 209);
    const auto *po_210 = buffer.data(po + 210);
    const auto *po_211 = buffer.data(po + 211);
    const auto *po_212 = buffer.data(po + 212);
    const auto *po_213 = buffer.data(po + 213);
    const auto *po_214 = buffer.data(po + 214);
    const auto *po_215 = buffer.data(po + 215);
    const auto *po_216 = buffer.data(po + 216);
    const auto *po_217 = buffer.data(po + 217);
    const auto *po_218 = buffer.data(po + 218);
    const auto *po_219 = buffer.data(po + 219);
    const auto *po_220 = buffer.data(po + 220);
    const auto *po_221 = buffer.data(po + 221);
    const auto *po_222 = buffer.data(po + 222);
    const auto *po_223 = buffer.data(po + 223);
    const auto *po_224 = buffer.data(po + 224);
    const auto *po_225 = buffer.data(po + 225);
    const auto *po_226 = buffer.data(po + 226);
    const auto *po_227 = buffer.data(po + 227);
    const auto *po_228 = buffer.data(po + 228);
    const auto *po_229 = buffer.data(po + 229);
    const auto *po_230 = buffer.data(po + 230);
    const auto *po_231 = buffer.data(po + 231);
    const auto *po_232 = buffer.data(po + 232);
    const auto *po_233 = buffer.data(po + 233);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pn_0, pn_1, pn_2, pn_3, pn_4, po_0, \
                         po_1, po_2, po_3, po_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * pn_0[k]
                 + po_0[k];

        t_1[k] = -ab_x[k] * pn_1[k]
                 + po_1[k];

        t_2[k] = -ab_x[k] * pn_2[k]
                 + po_2[k];

        t_3[k] = -ab_x[k] * pn_3[k]
                 + po_3[k];

        t_4[k] = -ab_x[k] * pn_4[k]
                 + po_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pn_5, pn_6, pn_7, pn_8, pn_9, po_5, \
                         po_6, po_7, po_8, po_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * pn_5[k]
                 + po_5[k];

        t_6[k] = -ab_x[k] * pn_6[k]
                 + po_6[k];

        t_7[k] = -ab_x[k] * pn_7[k]
                 + po_7[k];

        t_8[k] = -ab_x[k] * pn_8[k]
                 + po_8[k];

        t_9[k] = -ab_x[k] * pn_9[k]
                 + po_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pn_10, pn_11, pn_12, pn_13, \
                         pn_14, po_10, po_11, po_12, po_13, po_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * pn_10[k]
                  + po_10[k];

        t_11[k] = -ab_x[k] * pn_11[k]
                  + po_11[k];

        t_12[k] = -ab_x[k] * pn_12[k]
                  + po_12[k];

        t_13[k] = -ab_x[k] * pn_13[k]
                  + po_13[k];

        t_14[k] = -ab_x[k] * pn_14[k]
                  + po_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pn_15, pn_16, pn_17, pn_18, \
                         pn_19, po_15, po_16, po_17, po_18, po_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * pn_15[k]
                  + po_15[k];

        t_16[k] = -ab_x[k] * pn_16[k]
                  + po_16[k];

        t_17[k] = -ab_x[k] * pn_17[k]
                  + po_17[k];

        t_18[k] = -ab_x[k] * pn_18[k]
                  + po_18[k];

        t_19[k] = -ab_x[k] * pn_19[k]
                  + po_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pn_20, pn_21, pn_22, pn_23, \
                         pn_24, po_20, po_21, po_22, po_23, po_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * pn_20[k]
                  + po_20[k];

        t_21[k] = -ab_x[k] * pn_21[k]
                  + po_21[k];

        t_22[k] = -ab_x[k] * pn_22[k]
                  + po_22[k];

        t_23[k] = -ab_x[k] * pn_23[k]
                  + po_23[k];

        t_24[k] = -ab_x[k] * pn_24[k]
                  + po_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pn_25, pn_26, pn_27, pn_28, \
                         pn_29, po_25, po_26, po_27, po_28, po_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * pn_25[k]
                  + po_25[k];

        t_26[k] = -ab_x[k] * pn_26[k]
                  + po_26[k];

        t_27[k] = -ab_x[k] * pn_27[k]
                  + po_27[k];

        t_28[k] = -ab_x[k] * pn_28[k]
                  + po_28[k];

        t_29[k] = -ab_x[k] * pn_29[k]
                  + po_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pn_30, pn_31, pn_32, pn_33, \
                         pn_34, po_30, po_31, po_32, po_33, po_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * pn_30[k]
                  + po_30[k];

        t_31[k] = -ab_x[k] * pn_31[k]
                  + po_31[k];

        t_32[k] = -ab_x[k] * pn_32[k]
                  + po_32[k];

        t_33[k] = -ab_x[k] * pn_33[k]
                  + po_33[k];

        t_34[k] = -ab_x[k] * pn_34[k]
                  + po_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pn_35, pn_36, pn_37, pn_38, \
                         pn_39, po_35, po_36, po_37, po_38, po_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * pn_35[k]
                  + po_35[k];

        t_36[k] = -ab_x[k] * pn_36[k]
                  + po_36[k];

        t_37[k] = -ab_x[k] * pn_37[k]
                  + po_37[k];

        t_38[k] = -ab_x[k] * pn_38[k]
                  + po_38[k];

        t_39[k] = -ab_x[k] * pn_39[k]
                  + po_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pn_40, pn_41, pn_42, pn_43, \
                         pn_44, po_40, po_41, po_42, po_43, po_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * pn_40[k]
                  + po_40[k];

        t_41[k] = -ab_x[k] * pn_41[k]
                  + po_41[k];

        t_42[k] = -ab_x[k] * pn_42[k]
                  + po_42[k];

        t_43[k] = -ab_x[k] * pn_43[k]
                  + po_43[k];

        t_44[k] = -ab_x[k] * pn_44[k]
                  + po_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, pn_45, pn_46, pn_47, pn_48, \
                         pn_49, po_45, po_46, po_47, po_48, po_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * pn_45[k]
                  + po_45[k];

        t_46[k] = -ab_x[k] * pn_46[k]
                  + po_46[k];

        t_47[k] = -ab_x[k] * pn_47[k]
                  + po_47[k];

        t_48[k] = -ab_x[k] * pn_48[k]
                  + po_48[k];

        t_49[k] = -ab_x[k] * pn_49[k]
                  + po_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, pn_50, pn_51, pn_52, pn_53, \
                         pn_54, po_50, po_51, po_52, po_53, po_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * pn_50[k]
                  + po_50[k];

        t_51[k] = -ab_x[k] * pn_51[k]
                  + po_51[k];

        t_52[k] = -ab_x[k] * pn_52[k]
                  + po_52[k];

        t_53[k] = -ab_x[k] * pn_53[k]
                  + po_53[k];

        t_54[k] = -ab_x[k] * pn_54[k]
                  + po_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, pn_55, pn_56, pn_57, pn_58, \
                         pn_59, po_55, po_56, po_57, po_58, po_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * pn_55[k]
                  + po_55[k];

        t_56[k] = -ab_x[k] * pn_56[k]
                  + po_56[k];

        t_57[k] = -ab_x[k] * pn_57[k]
                  + po_57[k];

        t_58[k] = -ab_x[k] * pn_58[k]
                  + po_58[k];

        t_59[k] = -ab_x[k] * pn_59[k]
                  + po_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, pn_60, pn_61, pn_62, pn_63, \
                         pn_64, po_60, po_61, po_62, po_63, po_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * pn_60[k]
                  + po_60[k];

        t_61[k] = -ab_x[k] * pn_61[k]
                  + po_61[k];

        t_62[k] = -ab_x[k] * pn_62[k]
                  + po_62[k];

        t_63[k] = -ab_x[k] * pn_63[k]
                  + po_63[k];

        t_64[k] = -ab_x[k] * pn_64[k]
                  + po_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, pn_65, pn_66, pn_67, pn_68, \
                         pn_69, po_65, po_78, po_79, po_80, po_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * pn_65[k]
                  + po_65[k];

        t_66[k] = -ab_x[k] * pn_66[k]
                  + po_78[k];

        t_67[k] = -ab_x[k] * pn_67[k]
                  + po_79[k];

        t_68[k] = -ab_x[k] * pn_68[k]
                  + po_80[k];

        t_69[k] = -ab_x[k] * pn_69[k]
                  + po_81[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, pn_70, pn_71, pn_72, pn_73, \
                         pn_74, po_82, po_83, po_84, po_85, po_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * pn_70[k]
                  + po_82[k];

        t_71[k] = -ab_x[k] * pn_71[k]
                  + po_83[k];

        t_72[k] = -ab_x[k] * pn_72[k]
                  + po_84[k];

        t_73[k] = -ab_x[k] * pn_73[k]
                  + po_85[k];

        t_74[k] = -ab_x[k] * pn_74[k]
                  + po_86[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, pn_75, pn_76, pn_77, pn_78, \
                         pn_79, po_87, po_88, po_89, po_90, po_91 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * pn_75[k]
                  + po_87[k];

        t_76[k] = -ab_x[k] * pn_76[k]
                  + po_88[k];

        t_77[k] = -ab_x[k] * pn_77[k]
                  + po_89[k];

        t_78[k] = -ab_x[k] * pn_78[k]
                  + po_90[k];

        t_79[k] = -ab_x[k] * pn_79[k]
                  + po_91[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, pn_80, pn_81, pn_82, pn_83, \
                         pn_84, po_92, po_93, po_94, po_95, po_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * pn_80[k]
                  + po_92[k];

        t_81[k] = -ab_x[k] * pn_81[k]
                  + po_93[k];

        t_82[k] = -ab_x[k] * pn_82[k]
                  + po_94[k];

        t_83[k] = -ab_x[k] * pn_83[k]
                  + po_95[k];

        t_84[k] = -ab_x[k] * pn_84[k]
                  + po_96[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, pn_85, pn_86, pn_87, pn_88, \
                         pn_89, po_97, po_98, po_99, po_100, po_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * pn_85[k]
                  + po_97[k];

        t_86[k] = -ab_x[k] * pn_86[k]
                  + po_98[k];

        t_87[k] = -ab_x[k] * pn_87[k]
                  + po_99[k];

        t_88[k] = -ab_x[k] * pn_88[k]
                  + po_100[k];

        t_89[k] = -ab_x[k] * pn_89[k]
                  + po_101[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, pn_90, pn_91, pn_92, pn_93, \
                         pn_94, po_102, po_103, po_104, po_105, \
                         po_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * pn_90[k]
                  + po_102[k];

        t_91[k] = -ab_x[k] * pn_91[k]
                  + po_103[k];

        t_92[k] = -ab_x[k] * pn_92[k]
                  + po_104[k];

        t_93[k] = -ab_x[k] * pn_93[k]
                  + po_105[k];

        t_94[k] = -ab_x[k] * pn_94[k]
                  + po_106[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, pn_95, pn_96, pn_97, pn_98, \
                         pn_99, po_107, po_108, po_109, po_110, \
                         po_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * pn_95[k]
                  + po_107[k];

        t_96[k] = -ab_x[k] * pn_96[k]
                  + po_108[k];

        t_97[k] = -ab_x[k] * pn_97[k]
                  + po_109[k];

        t_98[k] = -ab_x[k] * pn_98[k]
                  + po_110[k];

        t_99[k] = -ab_x[k] * pn_99[k]
                  + po_111[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, pn_100, pn_101, pn_102, \
                         pn_103, pn_104, po_112, po_113, po_114, po_115, \
                         po_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * pn_100[k]
                   + po_112[k];

        t_101[k] = -ab_x[k] * pn_101[k]
                   + po_113[k];

        t_102[k] = -ab_x[k] * pn_102[k]
                   + po_114[k];

        t_103[k] = -ab_x[k] * pn_103[k]
                   + po_115[k];

        t_104[k] = -ab_x[k] * pn_104[k]
                   + po_116[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, pn_105, pn_106, pn_107, \
                         pn_108, pn_109, po_117, po_118, po_119, po_120, \
                         po_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * pn_105[k]
                   + po_117[k];

        t_106[k] = -ab_x[k] * pn_106[k]
                   + po_118[k];

        t_107[k] = -ab_x[k] * pn_107[k]
                   + po_119[k];

        t_108[k] = -ab_x[k] * pn_108[k]
                   + po_120[k];

        t_109[k] = -ab_x[k] * pn_109[k]
                   + po_121[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, pn_110, pn_111, pn_112, \
                         pn_113, pn_114, po_122, po_123, po_124, po_125, \
                         po_126 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * pn_110[k]
                   + po_122[k];

        t_111[k] = -ab_x[k] * pn_111[k]
                   + po_123[k];

        t_112[k] = -ab_x[k] * pn_112[k]
                   + po_124[k];

        t_113[k] = -ab_x[k] * pn_113[k]
                   + po_125[k];

        t_114[k] = -ab_x[k] * pn_114[k]
                   + po_126[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, pn_115, pn_116, pn_117, \
                         pn_118, pn_119, po_127, po_128, po_129, po_130, \
                         po_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * pn_115[k]
                   + po_127[k];

        t_116[k] = -ab_x[k] * pn_116[k]
                   + po_128[k];

        t_117[k] = -ab_x[k] * pn_117[k]
                   + po_129[k];

        t_118[k] = -ab_x[k] * pn_118[k]
                   + po_130[k];

        t_119[k] = -ab_x[k] * pn_119[k]
                   + po_131[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, pn_120, pn_121, pn_122, \
                         pn_123, pn_124, po_132, po_133, po_134, po_135, \
                         po_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * pn_120[k]
                   + po_132[k];

        t_121[k] = -ab_x[k] * pn_121[k]
                   + po_133[k];

        t_122[k] = -ab_x[k] * pn_122[k]
                   + po_134[k];

        t_123[k] = -ab_x[k] * pn_123[k]
                   + po_135[k];

        t_124[k] = -ab_x[k] * pn_124[k]
                   + po_136[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, pn_125, pn_126, pn_127, \
                         pn_128, pn_129, po_137, po_138, po_139, po_140, \
                         po_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * pn_125[k]
                   + po_137[k];

        t_126[k] = -ab_x[k] * pn_126[k]
                   + po_138[k];

        t_127[k] = -ab_x[k] * pn_127[k]
                   + po_139[k];

        t_128[k] = -ab_x[k] * pn_128[k]
                   + po_140[k];

        t_129[k] = -ab_x[k] * pn_129[k]
                   + po_141[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, pn_130, pn_131, pn_132, \
                         pn_133, pn_134, po_142, po_143, po_156, po_157, \
                         po_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_x[k] * pn_130[k]
                   + po_142[k];

        t_131[k] = -ab_x[k] * pn_131[k]
                   + po_143[k];

        t_132[k] = -ab_x[k] * pn_132[k]
                   + po_156[k];

        t_133[k] = -ab_x[k] * pn_133[k]
                   + po_157[k];

        t_134[k] = -ab_x[k] * pn_134[k]
                   + po_158[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, pn_135, pn_136, pn_137, \
                         pn_138, pn_139, po_159, po_160, po_161, po_162, \
                         po_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_x[k] * pn_135[k]
                   + po_159[k];

        t_136[k] = -ab_x[k] * pn_136[k]
                   + po_160[k];

        t_137[k] = -ab_x[k] * pn_137[k]
                   + po_161[k];

        t_138[k] = -ab_x[k] * pn_138[k]
                   + po_162[k];

        t_139[k] = -ab_x[k] * pn_139[k]
                   + po_163[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, pn_140, pn_141, pn_142, \
                         pn_143, pn_144, po_164, po_165, po_166, po_167, \
                         po_168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_x[k] * pn_140[k]
                   + po_164[k];

        t_141[k] = -ab_x[k] * pn_141[k]
                   + po_165[k];

        t_142[k] = -ab_x[k] * pn_142[k]
                   + po_166[k];

        t_143[k] = -ab_x[k] * pn_143[k]
                   + po_167[k];

        t_144[k] = -ab_x[k] * pn_144[k]
                   + po_168[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, pn_145, pn_146, pn_147, \
                         pn_148, pn_149, po_169, po_170, po_171, po_172, \
                         po_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_x[k] * pn_145[k]
                   + po_169[k];

        t_146[k] = -ab_x[k] * pn_146[k]
                   + po_170[k];

        t_147[k] = -ab_x[k] * pn_147[k]
                   + po_171[k];

        t_148[k] = -ab_x[k] * pn_148[k]
                   + po_172[k];

        t_149[k] = -ab_x[k] * pn_149[k]
                   + po_173[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, pn_150, pn_151, pn_152, \
                         pn_153, pn_154, po_174, po_175, po_176, po_177, \
                         po_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_x[k] * pn_150[k]
                   + po_174[k];

        t_151[k] = -ab_x[k] * pn_151[k]
                   + po_175[k];

        t_152[k] = -ab_x[k] * pn_152[k]
                   + po_176[k];

        t_153[k] = -ab_x[k] * pn_153[k]
                   + po_177[k];

        t_154[k] = -ab_x[k] * pn_154[k]
                   + po_178[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, pn_155, pn_156, pn_157, \
                         pn_158, pn_159, po_179, po_180, po_181, po_182, \
                         po_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_x[k] * pn_155[k]
                   + po_179[k];

        t_156[k] = -ab_x[k] * pn_156[k]
                   + po_180[k];

        t_157[k] = -ab_x[k] * pn_157[k]
                   + po_181[k];

        t_158[k] = -ab_x[k] * pn_158[k]
                   + po_182[k];

        t_159[k] = -ab_x[k] * pn_159[k]
                   + po_183[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, pn_160, pn_161, pn_162, \
                         pn_163, pn_164, po_184, po_185, po_186, po_187, \
                         po_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_x[k] * pn_160[k]
                   + po_184[k];

        t_161[k] = -ab_x[k] * pn_161[k]
                   + po_185[k];

        t_162[k] = -ab_x[k] * pn_162[k]
                   + po_186[k];

        t_163[k] = -ab_x[k] * pn_163[k]
                   + po_187[k];

        t_164[k] = -ab_x[k] * pn_164[k]
                   + po_188[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, pn_165, pn_166, pn_167, \
                         pn_168, pn_169, po_189, po_190, po_191, po_192, \
                         po_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = -ab_x[k] * pn_165[k]
                   + po_189[k];

        t_166[k] = -ab_x[k] * pn_166[k]
                   + po_190[k];

        t_167[k] = -ab_x[k] * pn_167[k]
                   + po_191[k];

        t_168[k] = -ab_x[k] * pn_168[k]
                   + po_192[k];

        t_169[k] = -ab_x[k] * pn_169[k]
                   + po_193[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, pn_170, pn_171, pn_172, \
                         pn_173, pn_174, po_194, po_195, po_196, po_197, \
                         po_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = -ab_x[k] * pn_170[k]
                   + po_194[k];

        t_171[k] = -ab_x[k] * pn_171[k]
                   + po_195[k];

        t_172[k] = -ab_x[k] * pn_172[k]
                   + po_196[k];

        t_173[k] = -ab_x[k] * pn_173[k]
                   + po_197[k];

        t_174[k] = -ab_x[k] * pn_174[k]
                   + po_198[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, pn_175, pn_176, pn_177, \
                         pn_178, pn_179, po_199, po_200, po_201, po_202, \
                         po_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = -ab_x[k] * pn_175[k]
                   + po_199[k];

        t_176[k] = -ab_x[k] * pn_176[k]
                   + po_200[k];

        t_177[k] = -ab_x[k] * pn_177[k]
                   + po_201[k];

        t_178[k] = -ab_x[k] * pn_178[k]
                   + po_202[k];

        t_179[k] = -ab_x[k] * pn_179[k]
                   + po_203[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, pn_180, pn_181, pn_182, \
                         pn_183, pn_184, po_204, po_205, po_206, po_207, \
                         po_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = -ab_x[k] * pn_180[k]
                   + po_204[k];

        t_181[k] = -ab_x[k] * pn_181[k]
                   + po_205[k];

        t_182[k] = -ab_x[k] * pn_182[k]
                   + po_206[k];

        t_183[k] = -ab_x[k] * pn_183[k]
                   + po_207[k];

        t_184[k] = -ab_x[k] * pn_184[k]
                   + po_208[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, pn_185, pn_186, pn_187, \
                         pn_188, pn_189, po_209, po_210, po_211, po_212, \
                         po_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = -ab_x[k] * pn_185[k]
                   + po_209[k];

        t_186[k] = -ab_x[k] * pn_186[k]
                   + po_210[k];

        t_187[k] = -ab_x[k] * pn_187[k]
                   + po_211[k];

        t_188[k] = -ab_x[k] * pn_188[k]
                   + po_212[k];

        t_189[k] = -ab_x[k] * pn_189[k]
                   + po_213[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, pn_190, pn_191, pn_192, \
                         pn_193, pn_194, po_214, po_215, po_216, po_217, \
                         po_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = -ab_x[k] * pn_190[k]
                   + po_214[k];

        t_191[k] = -ab_x[k] * pn_191[k]
                   + po_215[k];

        t_192[k] = -ab_x[k] * pn_192[k]
                   + po_216[k];

        t_193[k] = -ab_x[k] * pn_193[k]
                   + po_217[k];

        t_194[k] = -ab_x[k] * pn_194[k]
                   + po_218[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, ab_x, ab_y, pn_66, pn_195, pn_196, \
                         pn_197, po_79, po_219, po_220, po_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = -ab_x[k] * pn_195[k]
                   + po_219[k];

        t_196[k] = -ab_x[k] * pn_196[k]
                   + po_220[k];

        t_197[k] = -ab_x[k] * pn_197[k]
                   + po_221[k];

        t_198[k] = -ab_y[k] * pn_66[k]
                   + po_79[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, ab_y, pn_67, pn_68, pn_69, pn_70, \
                         pn_71, po_81, po_82, po_84, po_85, po_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_199[k] = -ab_y[k] * pn_67[k]
                   + po_81[k];

        t_200[k] = -ab_y[k] * pn_68[k]
                   + po_82[k];

        t_201[k] = -ab_y[k] * pn_69[k]
                   + po_84[k];

        t_202[k] = -ab_y[k] * pn_70[k]
                   + po_85[k];

        t_203[k] = -ab_y[k] * pn_71[k]
                   + po_86[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, ab_y, pn_72, pn_73, pn_74, pn_75, \
                         pn_76, po_88, po_89, po_90, po_91, po_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_204[k] = -ab_y[k] * pn_72[k]
                   + po_88[k];

        t_205[k] = -ab_y[k] * pn_73[k]
                   + po_89[k];

        t_206[k] = -ab_y[k] * pn_74[k]
                   + po_90[k];

        t_207[k] = -ab_y[k] * pn_75[k]
                   + po_91[k];

        t_208[k] = -ab_y[k] * pn_76[k]
                   + po_93[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, ab_y, pn_77, pn_78, pn_79, pn_80, \
                         pn_81, po_94, po_95, po_96, po_97, po_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_209[k] = -ab_y[k] * pn_77[k]
                   + po_94[k];

        t_210[k] = -ab_y[k] * pn_78[k]
                   + po_95[k];

        t_211[k] = -ab_y[k] * pn_79[k]
                   + po_96[k];

        t_212[k] = -ab_y[k] * pn_80[k]
                   + po_97[k];

        t_213[k] = -ab_y[k] * pn_81[k]
                   + po_99[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, ab_y, pn_82, pn_83, pn_84, pn_85, \
                         pn_86, po_100, po_101, po_102, po_103, \
                         po_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_214[k] = -ab_y[k] * pn_82[k]
                   + po_100[k];

        t_215[k] = -ab_y[k] * pn_83[k]
                   + po_101[k];

        t_216[k] = -ab_y[k] * pn_84[k]
                   + po_102[k];

        t_217[k] = -ab_y[k] * pn_85[k]
                   + po_103[k];

        t_218[k] = -ab_y[k] * pn_86[k]
                   + po_104[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, ab_y, pn_87, pn_88, pn_89, pn_90, \
                         pn_91, po_106, po_107, po_108, po_109, \
                         po_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_219[k] = -ab_y[k] * pn_87[k]
                   + po_106[k];

        t_220[k] = -ab_y[k] * pn_88[k]
                   + po_107[k];

        t_221[k] = -ab_y[k] * pn_89[k]
                   + po_108[k];

        t_222[k] = -ab_y[k] * pn_90[k]
                   + po_109[k];

        t_223[k] = -ab_y[k] * pn_91[k]
                   + po_110[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, ab_y, pn_92, pn_93, pn_94, pn_95, \
                         pn_96, po_111, po_112, po_114, po_115, \
                         po_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_224[k] = -ab_y[k] * pn_92[k]
                   + po_111[k];

        t_225[k] = -ab_y[k] * pn_93[k]
                   + po_112[k];

        t_226[k] = -ab_y[k] * pn_94[k]
                   + po_114[k];

        t_227[k] = -ab_y[k] * pn_95[k]
                   + po_115[k];

        t_228[k] = -ab_y[k] * pn_96[k]
                   + po_116[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, ab_y, pn_97, pn_98, pn_99, pn_100, \
                         pn_101, po_117, po_118, po_119, po_120, \
                         po_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_229[k] = -ab_y[k] * pn_97[k]
                   + po_117[k];

        t_230[k] = -ab_y[k] * pn_98[k]
                   + po_118[k];

        t_231[k] = -ab_y[k] * pn_99[k]
                   + po_119[k];

        t_232[k] = -ab_y[k] * pn_100[k]
                   + po_120[k];

        t_233[k] = -ab_y[k] * pn_101[k]
                   + po_121[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_y, pn_102, pn_103, pn_104, \
                         pn_105, pn_106, po_123, po_124, po_125, po_126, \
                         po_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_234[k] = -ab_y[k] * pn_102[k]
                   + po_123[k];

        t_235[k] = -ab_y[k] * pn_103[k]
                   + po_124[k];

        t_236[k] = -ab_y[k] * pn_104[k]
                   + po_125[k];

        t_237[k] = -ab_y[k] * pn_105[k]
                   + po_126[k];

        t_238[k] = -ab_y[k] * pn_106[k]
                   + po_127[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_y, pn_107, pn_108, pn_109, \
                         pn_110, pn_111, po_128, po_129, po_130, po_131, \
                         po_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_239[k] = -ab_y[k] * pn_107[k]
                   + po_128[k];

        t_240[k] = -ab_y[k] * pn_108[k]
                   + po_129[k];

        t_241[k] = -ab_y[k] * pn_109[k]
                   + po_130[k];

        t_242[k] = -ab_y[k] * pn_110[k]
                   + po_131[k];

        t_243[k] = -ab_y[k] * pn_111[k]
                   + po_133[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, ab_y, pn_112, pn_113, pn_114, \
                         pn_115, pn_116, po_134, po_135, po_136, po_137, \
                         po_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_244[k] = -ab_y[k] * pn_112[k]
                   + po_134[k];

        t_245[k] = -ab_y[k] * pn_113[k]
                   + po_135[k];

        t_246[k] = -ab_y[k] * pn_114[k]
                   + po_136[k];

        t_247[k] = -ab_y[k] * pn_115[k]
                   + po_137[k];

        t_248[k] = -ab_y[k] * pn_116[k]
                   + po_138[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, ab_y, pn_117, pn_118, pn_119, \
                         pn_120, pn_121, po_139, po_140, po_141, po_142, \
                         po_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_249[k] = -ab_y[k] * pn_117[k]
                   + po_139[k];

        t_250[k] = -ab_y[k] * pn_118[k]
                   + po_140[k];

        t_251[k] = -ab_y[k] * pn_119[k]
                   + po_141[k];

        t_252[k] = -ab_y[k] * pn_120[k]
                   + po_142[k];

        t_253[k] = -ab_y[k] * pn_121[k]
                   + po_144[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, ab_y, pn_122, pn_123, pn_124, \
                         pn_125, pn_126, po_145, po_146, po_147, po_148, \
                         po_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_254[k] = -ab_y[k] * pn_122[k]
                   + po_145[k];

        t_255[k] = -ab_y[k] * pn_123[k]
                   + po_146[k];

        t_256[k] = -ab_y[k] * pn_124[k]
                   + po_147[k];

        t_257[k] = -ab_y[k] * pn_125[k]
                   + po_148[k];

        t_258[k] = -ab_y[k] * pn_126[k]
                   + po_149[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, ab_y, pn_127, pn_128, pn_129, \
                         pn_130, pn_131, po_150, po_151, po_152, po_153, \
                         po_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_259[k] = -ab_y[k] * pn_127[k]
                   + po_150[k];

        t_260[k] = -ab_y[k] * pn_128[k]
                   + po_151[k];

        t_261[k] = -ab_y[k] * pn_129[k]
                   + po_152[k];

        t_262[k] = -ab_y[k] * pn_130[k]
                   + po_153[k];

        t_263[k] = -ab_y[k] * pn_131[k]
                   + po_154[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, ab_y, pn_132, pn_133, pn_134, \
                         pn_135, pn_136, po_157, po_159, po_160, po_162, \
                         po_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_264[k] = -ab_y[k] * pn_132[k]
                   + po_157[k];

        t_265[k] = -ab_y[k] * pn_133[k]
                   + po_159[k];

        t_266[k] = -ab_y[k] * pn_134[k]
                   + po_160[k];

        t_267[k] = -ab_y[k] * pn_135[k]
                   + po_162[k];

        t_268[k] = -ab_y[k] * pn_136[k]
                   + po_163[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, ab_y, pn_137, pn_138, pn_139, \
                         pn_140, pn_141, po_164, po_166, po_167, po_168, \
                         po_169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_269[k] = -ab_y[k] * pn_137[k]
                   + po_164[k];

        t_270[k] = -ab_y[k] * pn_138[k]
                   + po_166[k];

        t_271[k] = -ab_y[k] * pn_139[k]
                   + po_167[k];

        t_272[k] = -ab_y[k] * pn_140[k]
                   + po_168[k];

        t_273[k] = -ab_y[k] * pn_141[k]
                   + po_169[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, ab_y, pn_142, pn_143, pn_144, \
                         pn_145, pn_146, po_171, po_172, po_173, po_174, \
                         po_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_274[k] = -ab_y[k] * pn_142[k]
                   + po_171[k];

        t_275[k] = -ab_y[k] * pn_143[k]
                   + po_172[k];

        t_276[k] = -ab_y[k] * pn_144[k]
                   + po_173[k];

        t_277[k] = -ab_y[k] * pn_145[k]
                   + po_174[k];

        t_278[k] = -ab_y[k] * pn_146[k]
                   + po_175[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, t_283, ab_y, pn_147, pn_148, pn_149, \
                         pn_150, pn_151, po_177, po_178, po_179, po_180, \
                         po_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_279[k] = -ab_y[k] * pn_147[k]
                   + po_177[k];

        t_280[k] = -ab_y[k] * pn_148[k]
                   + po_178[k];

        t_281[k] = -ab_y[k] * pn_149[k]
                   + po_179[k];

        t_282[k] = -ab_y[k] * pn_150[k]
                   + po_180[k];

        t_283[k] = -ab_y[k] * pn_151[k]
                   + po_181[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, ab_y, pn_152, pn_153, pn_154, \
                         pn_155, pn_156, po_182, po_184, po_185, po_186, \
                         po_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_284[k] = -ab_y[k] * pn_152[k]
                   + po_182[k];

        t_285[k] = -ab_y[k] * pn_153[k]
                   + po_184[k];

        t_286[k] = -ab_y[k] * pn_154[k]
                   + po_185[k];

        t_287[k] = -ab_y[k] * pn_155[k]
                   + po_186[k];

        t_288[k] = -ab_y[k] * pn_156[k]
                   + po_187[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, ab_y, pn_157, pn_158, pn_159, \
                         pn_160, pn_161, po_188, po_189, po_190, po_192, \
                         po_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_289[k] = -ab_y[k] * pn_157[k]
                   + po_188[k];

        t_290[k] = -ab_y[k] * pn_158[k]
                   + po_189[k];

        t_291[k] = -ab_y[k] * pn_159[k]
                   + po_190[k];

        t_292[k] = -ab_y[k] * pn_160[k]
                   + po_192[k];

        t_293[k] = -ab_y[k] * pn_161[k]
                   + po_193[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_y, pn_162, pn_163, pn_164, \
                         pn_165, pn_166, po_194, po_195, po_196, po_197, \
                         po_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_294[k] = -ab_y[k] * pn_162[k]
                   + po_194[k];

        t_295[k] = -ab_y[k] * pn_163[k]
                   + po_195[k];

        t_296[k] = -ab_y[k] * pn_164[k]
                   + po_196[k];

        t_297[k] = -ab_y[k] * pn_165[k]
                   + po_197[k];

        t_298[k] = -ab_y[k] * pn_166[k]
                   + po_198[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_y, pn_167, pn_168, pn_169, \
                         pn_170, pn_171, po_199, po_201, po_202, po_203, \
                         po_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_299[k] = -ab_y[k] * pn_167[k]
                   + po_199[k];

        t_300[k] = -ab_y[k] * pn_168[k]
                   + po_201[k];

        t_301[k] = -ab_y[k] * pn_169[k]
                   + po_202[k];

        t_302[k] = -ab_y[k] * pn_170[k]
                   + po_203[k];

        t_303[k] = -ab_y[k] * pn_171[k]
                   + po_204[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_y, pn_172, pn_173, pn_174, \
                         pn_175, pn_176, po_205, po_206, po_207, po_208, \
                         po_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_304[k] = -ab_y[k] * pn_172[k]
                   + po_205[k];

        t_305[k] = -ab_y[k] * pn_173[k]
                   + po_206[k];

        t_306[k] = -ab_y[k] * pn_174[k]
                   + po_207[k];

        t_307[k] = -ab_y[k] * pn_175[k]
                   + po_208[k];

        t_308[k] = -ab_y[k] * pn_176[k]
                   + po_209[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_y, pn_177, pn_178, pn_179, \
                         pn_180, pn_181, po_211, po_212, po_213, po_214, \
                         po_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_309[k] = -ab_y[k] * pn_177[k]
                   + po_211[k];

        t_310[k] = -ab_y[k] * pn_178[k]
                   + po_212[k];

        t_311[k] = -ab_y[k] * pn_179[k]
                   + po_213[k];

        t_312[k] = -ab_y[k] * pn_180[k]
                   + po_214[k];

        t_313[k] = -ab_y[k] * pn_181[k]
                   + po_215[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, ab_y, pn_182, pn_183, pn_184, \
                         pn_185, pn_186, po_216, po_217, po_218, po_219, \
                         po_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_314[k] = -ab_y[k] * pn_182[k]
                   + po_216[k];

        t_315[k] = -ab_y[k] * pn_183[k]
                   + po_217[k];

        t_316[k] = -ab_y[k] * pn_184[k]
                   + po_218[k];

        t_317[k] = -ab_y[k] * pn_185[k]
                   + po_219[k];

        t_318[k] = -ab_y[k] * pn_186[k]
                   + po_220[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, ab_y, pn_187, pn_188, pn_189, \
                         pn_190, pn_191, po_222, po_223, po_224, po_225, \
                         po_226 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_319[k] = -ab_y[k] * pn_187[k]
                   + po_222[k];

        t_320[k] = -ab_y[k] * pn_188[k]
                   + po_223[k];

        t_321[k] = -ab_y[k] * pn_189[k]
                   + po_224[k];

        t_322[k] = -ab_y[k] * pn_190[k]
                   + po_225[k];

        t_323[k] = -ab_y[k] * pn_191[k]
                   + po_226[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, ab_y, pn_192, pn_193, pn_194, \
                         pn_195, pn_196, po_227, po_228, po_229, po_230, \
                         po_231 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_324[k] = -ab_y[k] * pn_192[k]
                   + po_227[k];

        t_325[k] = -ab_y[k] * pn_193[k]
                   + po_228[k];

        t_326[k] = -ab_y[k] * pn_194[k]
                   + po_229[k];

        t_327[k] = -ab_y[k] * pn_195[k]
                   + po_230[k];

        t_328[k] = -ab_y[k] * pn_196[k]
                   + po_231[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, ab_y, ab_z, pn_132, pn_133, pn_134, \
                         pn_197, po_158, po_160, po_161, po_232 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_329[k] = -ab_y[k] * pn_197[k]
                   + po_232[k];

        t_330[k] = -ab_z[k] * pn_132[k]
                   + po_158[k];

        t_331[k] = -ab_z[k] * pn_133[k]
                   + po_160[k];

        t_332[k] = -ab_z[k] * pn_134[k]
                   + po_161[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, ab_z, pn_135, pn_136, pn_137, \
                         pn_138, pn_139, po_163, po_164, po_165, po_167, \
                         po_168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_333[k] = -ab_z[k] * pn_135[k]
                   + po_163[k];

        t_334[k] = -ab_z[k] * pn_136[k]
                   + po_164[k];

        t_335[k] = -ab_z[k] * pn_137[k]
                   + po_165[k];

        t_336[k] = -ab_z[k] * pn_138[k]
                   + po_167[k];

        t_337[k] = -ab_z[k] * pn_139[k]
                   + po_168[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, ab_z, pn_140, pn_141, pn_142, \
                         pn_143, pn_144, po_169, po_170, po_172, po_173, \
                         po_174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_338[k] = -ab_z[k] * pn_140[k]
                   + po_169[k];

        t_339[k] = -ab_z[k] * pn_141[k]
                   + po_170[k];

        t_340[k] = -ab_z[k] * pn_142[k]
                   + po_172[k];

        t_341[k] = -ab_z[k] * pn_143[k]
                   + po_173[k];

        t_342[k] = -ab_z[k] * pn_144[k]
                   + po_174[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, ab_z, pn_145, pn_146, pn_147, \
                         pn_148, pn_149, po_175, po_176, po_178, po_179, \
                         po_180 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_343[k] = -ab_z[k] * pn_145[k]
                   + po_175[k];

        t_344[k] = -ab_z[k] * pn_146[k]
                   + po_176[k];

        t_345[k] = -ab_z[k] * pn_147[k]
                   + po_178[k];

        t_346[k] = -ab_z[k] * pn_148[k]
                   + po_179[k];

        t_347[k] = -ab_z[k] * pn_149[k]
                   + po_180[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, ab_z, pn_150, pn_151, pn_152, \
                         pn_153, pn_154, po_181, po_182, po_183, po_185, \
                         po_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_348[k] = -ab_z[k] * pn_150[k]
                   + po_181[k];

        t_349[k] = -ab_z[k] * pn_151[k]
                   + po_182[k];

        t_350[k] = -ab_z[k] * pn_152[k]
                   + po_183[k];

        t_351[k] = -ab_z[k] * pn_153[k]
                   + po_185[k];

        t_352[k] = -ab_z[k] * pn_154[k]
                   + po_186[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, ab_z, pn_155, pn_156, pn_157, \
                         pn_158, pn_159, po_187, po_188, po_189, po_190, \
                         po_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_353[k] = -ab_z[k] * pn_155[k]
                   + po_187[k];

        t_354[k] = -ab_z[k] * pn_156[k]
                   + po_188[k];

        t_355[k] = -ab_z[k] * pn_157[k]
                   + po_189[k];

        t_356[k] = -ab_z[k] * pn_158[k]
                   + po_190[k];

        t_357[k] = -ab_z[k] * pn_159[k]
                   + po_191[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, ab_z, pn_160, pn_161, pn_162, \
                         pn_163, pn_164, po_193, po_194, po_195, po_196, \
                         po_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_358[k] = -ab_z[k] * pn_160[k]
                   + po_193[k];

        t_359[k] = -ab_z[k] * pn_161[k]
                   + po_194[k];

        t_360[k] = -ab_z[k] * pn_162[k]
                   + po_195[k];

        t_361[k] = -ab_z[k] * pn_163[k]
                   + po_196[k];

        t_362[k] = -ab_z[k] * pn_164[k]
                   + po_197[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, ab_z, pn_165, pn_166, pn_167, \
                         pn_168, pn_169, po_198, po_199, po_200, po_202, \
                         po_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_363[k] = -ab_z[k] * pn_165[k]
                   + po_198[k];

        t_364[k] = -ab_z[k] * pn_166[k]
                   + po_199[k];

        t_365[k] = -ab_z[k] * pn_167[k]
                   + po_200[k];

        t_366[k] = -ab_z[k] * pn_168[k]
                   + po_202[k];

        t_367[k] = -ab_z[k] * pn_169[k]
                   + po_203[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, ab_z, pn_170, pn_171, pn_172, \
                         pn_173, pn_174, po_204, po_205, po_206, po_207, \
                         po_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_368[k] = -ab_z[k] * pn_170[k]
                   + po_204[k];

        t_369[k] = -ab_z[k] * pn_171[k]
                   + po_205[k];

        t_370[k] = -ab_z[k] * pn_172[k]
                   + po_206[k];

        t_371[k] = -ab_z[k] * pn_173[k]
                   + po_207[k];

        t_372[k] = -ab_z[k] * pn_174[k]
                   + po_208[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, ab_z, pn_175, pn_176, pn_177, \
                         pn_178, pn_179, po_209, po_210, po_212, po_213, \
                         po_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_373[k] = -ab_z[k] * pn_175[k]
                   + po_209[k];

        t_374[k] = -ab_z[k] * pn_176[k]
                   + po_210[k];

        t_375[k] = -ab_z[k] * pn_177[k]
                   + po_212[k];

        t_376[k] = -ab_z[k] * pn_178[k]
                   + po_213[k];

        t_377[k] = -ab_z[k] * pn_179[k]
                   + po_214[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, ab_z, pn_180, pn_181, pn_182, \
                         pn_183, pn_184, po_215, po_216, po_217, po_218, \
                         po_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_378[k] = -ab_z[k] * pn_180[k]
                   + po_215[k];

        t_379[k] = -ab_z[k] * pn_181[k]
                   + po_216[k];

        t_380[k] = -ab_z[k] * pn_182[k]
                   + po_217[k];

        t_381[k] = -ab_z[k] * pn_183[k]
                   + po_218[k];

        t_382[k] = -ab_z[k] * pn_184[k]
                   + po_219[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, ab_z, pn_185, pn_186, pn_187, \
                         pn_188, pn_189, po_220, po_221, po_223, po_224, \
                         po_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_383[k] = -ab_z[k] * pn_185[k]
                   + po_220[k];

        t_384[k] = -ab_z[k] * pn_186[k]
                   + po_221[k];

        t_385[k] = -ab_z[k] * pn_187[k]
                   + po_223[k];

        t_386[k] = -ab_z[k] * pn_188[k]
                   + po_224[k];

        t_387[k] = -ab_z[k] * pn_189[k]
                   + po_225[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, ab_z, pn_190, pn_191, pn_192, \
                         pn_193, pn_194, po_226, po_227, po_228, po_229, \
                         po_230 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_388[k] = -ab_z[k] * pn_190[k]
                   + po_226[k];

        t_389[k] = -ab_z[k] * pn_191[k]
                   + po_227[k];

        t_390[k] = -ab_z[k] * pn_192[k]
                   + po_228[k];

        t_391[k] = -ab_z[k] * pn_193[k]
                   + po_229[k];

        t_392[k] = -ab_z[k] * pn_194[k]
                   + po_230[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, ab_z, pn_195, pn_196, pn_197, po_231, po_232, \
                         po_233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_393[k] = -ab_z[k] * pn_195[k]
                   + po_231[k];

        t_394[k] = -ab_z[k] * pn_196[k]
                   + po_232[k];

        t_395[k] = -ab_z[k] * pn_197[k]
                   + po_233[k];
    }
}

}  // namespace simdtrf
