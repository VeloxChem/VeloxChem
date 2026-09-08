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


#include "SimdTransferKF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_kf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t kd, const size_t ld, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_169 = buffer.data(kd + 169);
    const auto *kd_170 = buffer.data(kd + 170);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_211 = buffer.data(kd + 211);
    const auto *kd_212 = buffer.data(kd + 212);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

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
    const auto *ld_177 = buffer.data(ld + 177);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_181 = buffer.data(ld + 181);
    const auto *ld_182 = buffer.data(ld + 182);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_187 = buffer.data(ld + 187);
    const auto *ld_188 = buffer.data(ld + 188);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_193 = buffer.data(ld + 193);
    const auto *ld_194 = buffer.data(ld + 194);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_199 = buffer.data(ld + 199);
    const auto *ld_200 = buffer.data(ld + 200);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_205 = buffer.data(ld + 205);
    const auto *ld_206 = buffer.data(ld + 206);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_209 = buffer.data(ld + 209);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_211 = buffer.data(ld + 211);
    const auto *ld_212 = buffer.data(ld + 212);
    const auto *ld_213 = buffer.data(ld + 213);
    const auto *ld_214 = buffer.data(ld + 214);
    const auto *ld_215 = buffer.data(ld + 215);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_220 = buffer.data(ld + 220);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_269 = buffer.data(ld + 269);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, kd_0, kd_1, kd_2, kd_3, kd_4, ld_0, \
                         ld_1, ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * kd_0[k]
                 + ld_0[k];

        t_1[k] = ab_x[k] * kd_1[k]
                 + ld_1[k];

        t_2[k] = ab_x[k] * kd_2[k]
                 + ld_2[k];

        t_3[k] = ab_x[k] * kd_3[k]
                 + ld_3[k];

        t_4[k] = ab_x[k] * kd_4[k]
                 + ld_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, kd_3, kd_4, kd_5, ld_5, \
                         ld_9, ld_10, ld_11, ld_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_x[k] * kd_5[k]
                 + ld_5[k];

        t_6[k] = ab_y[k] * kd_3[k]
                 + ld_9[k];

        t_7[k] = ab_y[k] * kd_4[k]
                 + ld_10[k];

        t_8[k] = ab_y[k] * kd_5[k]
                 + ld_11[k];

        t_9[k] = ab_z[k] * kd_5[k]
                 + ld_17[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, kd_6, kd_7, kd_8, kd_9, kd_10, \
                         ld_6, ld_7, ld_8, ld_9, ld_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = ab_x[k] * kd_6[k]
                  + ld_6[k];

        t_11[k] = ab_x[k] * kd_7[k]
                  + ld_7[k];

        t_12[k] = ab_x[k] * kd_8[k]
                  + ld_8[k];

        t_13[k] = ab_x[k] * kd_9[k]
                  + ld_9[k];

        t_14[k] = ab_x[k] * kd_10[k]
                  + ld_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, kd_9, kd_10, kd_11, \
                         ld_11, ld_21, ld_22, ld_23, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = ab_x[k] * kd_11[k]
                  + ld_11[k];

        t_16[k] = ab_y[k] * kd_9[k]
                  + ld_21[k];

        t_17[k] = ab_y[k] * kd_10[k]
                  + ld_22[k];

        t_18[k] = ab_y[k] * kd_11[k]
                  + ld_23[k];

        t_19[k] = ab_z[k] * kd_11[k]
                  + ld_29[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, kd_12, kd_13, kd_14, kd_15, \
                         kd_16, ld_12, ld_13, ld_14, ld_15, ld_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = ab_x[k] * kd_12[k]
                  + ld_12[k];

        t_21[k] = ab_x[k] * kd_13[k]
                  + ld_13[k];

        t_22[k] = ab_x[k] * kd_14[k]
                  + ld_14[k];

        t_23[k] = ab_x[k] * kd_15[k]
                  + ld_15[k];

        t_24[k] = ab_x[k] * kd_16[k]
                  + ld_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, kd_15, kd_16, kd_17, \
                         ld_17, ld_27, ld_28, ld_29, ld_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = ab_x[k] * kd_17[k]
                  + ld_17[k];

        t_26[k] = ab_y[k] * kd_15[k]
                  + ld_27[k];

        t_27[k] = ab_y[k] * kd_16[k]
                  + ld_28[k];

        t_28[k] = ab_y[k] * kd_17[k]
                  + ld_29[k];

        t_29[k] = ab_z[k] * kd_17[k]
                  + ld_35[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, kd_18, kd_19, kd_20, kd_21, \
                         kd_22, ld_18, ld_19, ld_20, ld_21, ld_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = ab_x[k] * kd_18[k]
                  + ld_18[k];

        t_31[k] = ab_x[k] * kd_19[k]
                  + ld_19[k];

        t_32[k] = ab_x[k] * kd_20[k]
                  + ld_20[k];

        t_33[k] = ab_x[k] * kd_21[k]
                  + ld_21[k];

        t_34[k] = ab_x[k] * kd_22[k]
                  + ld_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, kd_21, kd_22, kd_23, \
                         ld_23, ld_39, ld_40, ld_41, ld_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = ab_x[k] * kd_23[k]
                  + ld_23[k];

        t_36[k] = ab_y[k] * kd_21[k]
                  + ld_39[k];

        t_37[k] = ab_y[k] * kd_22[k]
                  + ld_40[k];

        t_38[k] = ab_y[k] * kd_23[k]
                  + ld_41[k];

        t_39[k] = ab_z[k] * kd_23[k]
                  + ld_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, kd_24, kd_25, kd_26, kd_27, \
                         kd_28, ld_24, ld_25, ld_26, ld_27, ld_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = ab_x[k] * kd_24[k]
                  + ld_24[k];

        t_41[k] = ab_x[k] * kd_25[k]
                  + ld_25[k];

        t_42[k] = ab_x[k] * kd_26[k]
                  + ld_26[k];

        t_43[k] = ab_x[k] * kd_27[k]
                  + ld_27[k];

        t_44[k] = ab_x[k] * kd_28[k]
                  + ld_28[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, kd_27, kd_28, kd_29, \
                         ld_29, ld_45, ld_46, ld_47, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_x[k] * kd_29[k]
                  + ld_29[k];

        t_46[k] = ab_y[k] * kd_27[k]
                  + ld_45[k];

        t_47[k] = ab_y[k] * kd_28[k]
                  + ld_46[k];

        t_48[k] = ab_y[k] * kd_29[k]
                  + ld_47[k];

        t_49[k] = ab_z[k] * kd_29[k]
                  + ld_53[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, kd_30, kd_31, kd_32, kd_33, \
                         kd_34, ld_30, ld_31, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = ab_x[k] * kd_30[k]
                  + ld_30[k];

        t_51[k] = ab_x[k] * kd_31[k]
                  + ld_31[k];

        t_52[k] = ab_x[k] * kd_32[k]
                  + ld_32[k];

        t_53[k] = ab_x[k] * kd_33[k]
                  + ld_33[k];

        t_54[k] = ab_x[k] * kd_34[k]
                  + ld_34[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, kd_33, kd_34, kd_35, \
                         ld_35, ld_51, ld_52, ld_53, ld_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = ab_x[k] * kd_35[k]
                  + ld_35[k];

        t_56[k] = ab_y[k] * kd_33[k]
                  + ld_51[k];

        t_57[k] = ab_y[k] * kd_34[k]
                  + ld_52[k];

        t_58[k] = ab_y[k] * kd_35[k]
                  + ld_53[k];

        t_59[k] = ab_z[k] * kd_35[k]
                  + ld_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, kd_36, kd_37, kd_38, kd_39, \
                         kd_40, ld_36, ld_37, ld_38, ld_39, ld_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * kd_36[k]
                  + ld_36[k];

        t_61[k] = ab_x[k] * kd_37[k]
                  + ld_37[k];

        t_62[k] = ab_x[k] * kd_38[k]
                  + ld_38[k];

        t_63[k] = ab_x[k] * kd_39[k]
                  + ld_39[k];

        t_64[k] = ab_x[k] * kd_40[k]
                  + ld_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, kd_39, kd_40, kd_41, \
                         ld_41, ld_63, ld_64, ld_65, ld_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_x[k] * kd_41[k]
                  + ld_41[k];

        t_66[k] = ab_y[k] * kd_39[k]
                  + ld_63[k];

        t_67[k] = ab_y[k] * kd_40[k]
                  + ld_64[k];

        t_68[k] = ab_y[k] * kd_41[k]
                  + ld_65[k];

        t_69[k] = ab_z[k] * kd_41[k]
                  + ld_71[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, kd_42, kd_43, kd_44, kd_45, \
                         kd_46, ld_42, ld_43, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = ab_x[k] * kd_42[k]
                  + ld_42[k];

        t_71[k] = ab_x[k] * kd_43[k]
                  + ld_43[k];

        t_72[k] = ab_x[k] * kd_44[k]
                  + ld_44[k];

        t_73[k] = ab_x[k] * kd_45[k]
                  + ld_45[k];

        t_74[k] = ab_x[k] * kd_46[k]
                  + ld_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, kd_45, kd_46, kd_47, \
                         ld_47, ld_69, ld_70, ld_71, ld_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = ab_x[k] * kd_47[k]
                  + ld_47[k];

        t_76[k] = ab_y[k] * kd_45[k]
                  + ld_69[k];

        t_77[k] = ab_y[k] * kd_46[k]
                  + ld_70[k];

        t_78[k] = ab_y[k] * kd_47[k]
                  + ld_71[k];

        t_79[k] = ab_z[k] * kd_47[k]
                  + ld_77[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, kd_48, kd_49, kd_50, kd_51, \
                         kd_52, ld_48, ld_49, ld_50, ld_51, ld_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = ab_x[k] * kd_48[k]
                  + ld_48[k];

        t_81[k] = ab_x[k] * kd_49[k]
                  + ld_49[k];

        t_82[k] = ab_x[k] * kd_50[k]
                  + ld_50[k];

        t_83[k] = ab_x[k] * kd_51[k]
                  + ld_51[k];

        t_84[k] = ab_x[k] * kd_52[k]
                  + ld_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, kd_51, kd_52, kd_53, \
                         ld_53, ld_75, ld_76, ld_77, ld_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * kd_53[k]
                  + ld_53[k];

        t_86[k] = ab_y[k] * kd_51[k]
                  + ld_75[k];

        t_87[k] = ab_y[k] * kd_52[k]
                  + ld_76[k];

        t_88[k] = ab_y[k] * kd_53[k]
                  + ld_77[k];

        t_89[k] = ab_z[k] * kd_53[k]
                  + ld_83[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, kd_54, kd_55, kd_56, kd_57, \
                         kd_58, ld_54, ld_55, ld_56, ld_57, ld_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * kd_54[k]
                  + ld_54[k];

        t_91[k] = ab_x[k] * kd_55[k]
                  + ld_55[k];

        t_92[k] = ab_x[k] * kd_56[k]
                  + ld_56[k];

        t_93[k] = ab_x[k] * kd_57[k]
                  + ld_57[k];

        t_94[k] = ab_x[k] * kd_58[k]
                  + ld_58[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, kd_57, kd_58, kd_59, \
                         ld_59, ld_81, ld_82, ld_83, ld_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_x[k] * kd_59[k]
                  + ld_59[k];

        t_96[k] = ab_y[k] * kd_57[k]
                  + ld_81[k];

        t_97[k] = ab_y[k] * kd_58[k]
                  + ld_82[k];

        t_98[k] = ab_y[k] * kd_59[k]
                  + ld_83[k];

        t_99[k] = ab_z[k] * kd_59[k]
                  + ld_89[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, kd_60, kd_61, kd_62, kd_63, \
                         kd_64, ld_60, ld_61, ld_62, ld_63, ld_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = ab_x[k] * kd_60[k]
                   + ld_60[k];

        t_101[k] = ab_x[k] * kd_61[k]
                   + ld_61[k];

        t_102[k] = ab_x[k] * kd_62[k]
                   + ld_62[k];

        t_103[k] = ab_x[k] * kd_63[k]
                   + ld_63[k];

        t_104[k] = ab_x[k] * kd_64[k]
                   + ld_64[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, kd_63, kd_64, \
                         kd_65, ld_65, ld_93, ld_94, ld_95, ld_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = ab_x[k] * kd_65[k]
                   + ld_65[k];

        t_106[k] = ab_y[k] * kd_63[k]
                   + ld_93[k];

        t_107[k] = ab_y[k] * kd_64[k]
                   + ld_94[k];

        t_108[k] = ab_y[k] * kd_65[k]
                   + ld_95[k];

        t_109[k] = ab_z[k] * kd_65[k]
                   + ld_101[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, kd_66, kd_67, kd_68, kd_69, \
                         kd_70, ld_66, ld_67, ld_68, ld_69, ld_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = ab_x[k] * kd_66[k]
                   + ld_66[k];

        t_111[k] = ab_x[k] * kd_67[k]
                   + ld_67[k];

        t_112[k] = ab_x[k] * kd_68[k]
                   + ld_68[k];

        t_113[k] = ab_x[k] * kd_69[k]
                   + ld_69[k];

        t_114[k] = ab_x[k] * kd_70[k]
                   + ld_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, kd_69, kd_70, \
                         kd_71, ld_71, ld_99, ld_100, ld_101, ld_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = ab_x[k] * kd_71[k]
                   + ld_71[k];

        t_116[k] = ab_y[k] * kd_69[k]
                   + ld_99[k];

        t_117[k] = ab_y[k] * kd_70[k]
                   + ld_100[k];

        t_118[k] = ab_y[k] * kd_71[k]
                   + ld_101[k];

        t_119[k] = ab_z[k] * kd_71[k]
                   + ld_107[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, kd_72, kd_73, kd_74, kd_75, \
                         kd_76, ld_72, ld_73, ld_74, ld_75, ld_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = ab_x[k] * kd_72[k]
                   + ld_72[k];

        t_121[k] = ab_x[k] * kd_73[k]
                   + ld_73[k];

        t_122[k] = ab_x[k] * kd_74[k]
                   + ld_74[k];

        t_123[k] = ab_x[k] * kd_75[k]
                   + ld_75[k];

        t_124[k] = ab_x[k] * kd_76[k]
                   + ld_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, kd_75, kd_76, \
                         kd_77, ld_77, ld_105, ld_106, ld_107, ld_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = ab_x[k] * kd_77[k]
                   + ld_77[k];

        t_126[k] = ab_y[k] * kd_75[k]
                   + ld_105[k];

        t_127[k] = ab_y[k] * kd_76[k]
                   + ld_106[k];

        t_128[k] = ab_y[k] * kd_77[k]
                   + ld_107[k];

        t_129[k] = ab_z[k] * kd_77[k]
                   + ld_113[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, kd_78, kd_79, kd_80, kd_81, \
                         kd_82, ld_78, ld_79, ld_80, ld_81, ld_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = ab_x[k] * kd_78[k]
                   + ld_78[k];

        t_131[k] = ab_x[k] * kd_79[k]
                   + ld_79[k];

        t_132[k] = ab_x[k] * kd_80[k]
                   + ld_80[k];

        t_133[k] = ab_x[k] * kd_81[k]
                   + ld_81[k];

        t_134[k] = ab_x[k] * kd_82[k]
                   + ld_82[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, kd_81, kd_82, \
                         kd_83, ld_83, ld_111, ld_112, ld_113, ld_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_x[k] * kd_83[k]
                   + ld_83[k];

        t_136[k] = ab_y[k] * kd_81[k]
                   + ld_111[k];

        t_137[k] = ab_y[k] * kd_82[k]
                   + ld_112[k];

        t_138[k] = ab_y[k] * kd_83[k]
                   + ld_113[k];

        t_139[k] = ab_z[k] * kd_83[k]
                   + ld_119[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, kd_84, kd_85, kd_86, kd_87, \
                         kd_88, ld_84, ld_85, ld_86, ld_87, ld_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = ab_x[k] * kd_84[k]
                   + ld_84[k];

        t_141[k] = ab_x[k] * kd_85[k]
                   + ld_85[k];

        t_142[k] = ab_x[k] * kd_86[k]
                   + ld_86[k];

        t_143[k] = ab_x[k] * kd_87[k]
                   + ld_87[k];

        t_144[k] = ab_x[k] * kd_88[k]
                   + ld_88[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, kd_87, kd_88, \
                         kd_89, ld_89, ld_117, ld_118, ld_119, ld_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = ab_x[k] * kd_89[k]
                   + ld_89[k];

        t_146[k] = ab_y[k] * kd_87[k]
                   + ld_117[k];

        t_147[k] = ab_y[k] * kd_88[k]
                   + ld_118[k];

        t_148[k] = ab_y[k] * kd_89[k]
                   + ld_119[k];

        t_149[k] = ab_z[k] * kd_89[k]
                   + ld_125[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, kd_90, kd_91, kd_92, kd_93, \
                         kd_94, ld_90, ld_91, ld_92, ld_93, ld_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = ab_x[k] * kd_90[k]
                   + ld_90[k];

        t_151[k] = ab_x[k] * kd_91[k]
                   + ld_91[k];

        t_152[k] = ab_x[k] * kd_92[k]
                   + ld_92[k];

        t_153[k] = ab_x[k] * kd_93[k]
                   + ld_93[k];

        t_154[k] = ab_x[k] * kd_94[k]
                   + ld_94[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, kd_93, kd_94, \
                         kd_95, ld_95, ld_129, ld_130, ld_131, ld_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = ab_x[k] * kd_95[k]
                   + ld_95[k];

        t_156[k] = ab_y[k] * kd_93[k]
                   + ld_129[k];

        t_157[k] = ab_y[k] * kd_94[k]
                   + ld_130[k];

        t_158[k] = ab_y[k] * kd_95[k]
                   + ld_131[k];

        t_159[k] = ab_z[k] * kd_95[k]
                   + ld_137[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, kd_96, kd_97, kd_98, kd_99, \
                         kd_100, ld_96, ld_97, ld_98, ld_99, ld_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = ab_x[k] * kd_96[k]
                   + ld_96[k];

        t_161[k] = ab_x[k] * kd_97[k]
                   + ld_97[k];

        t_162[k] = ab_x[k] * kd_98[k]
                   + ld_98[k];

        t_163[k] = ab_x[k] * kd_99[k]
                   + ld_99[k];

        t_164[k] = ab_x[k] * kd_100[k]
                   + ld_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, kd_99, kd_100, \
                         kd_101, ld_101, ld_135, ld_136, ld_137, \
                         ld_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = ab_x[k] * kd_101[k]
                   + ld_101[k];

        t_166[k] = ab_y[k] * kd_99[k]
                   + ld_135[k];

        t_167[k] = ab_y[k] * kd_100[k]
                   + ld_136[k];

        t_168[k] = ab_y[k] * kd_101[k]
                   + ld_137[k];

        t_169[k] = ab_z[k] * kd_101[k]
                   + ld_143[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, kd_102, kd_103, kd_104, \
                         kd_105, kd_106, ld_102, ld_103, ld_104, ld_105, \
                         ld_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = ab_x[k] * kd_102[k]
                   + ld_102[k];

        t_171[k] = ab_x[k] * kd_103[k]
                   + ld_103[k];

        t_172[k] = ab_x[k] * kd_104[k]
                   + ld_104[k];

        t_173[k] = ab_x[k] * kd_105[k]
                   + ld_105[k];

        t_174[k] = ab_x[k] * kd_106[k]
                   + ld_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, kd_105, kd_106, \
                         kd_107, ld_107, ld_141, ld_142, ld_143, \
                         ld_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = ab_x[k] * kd_107[k]
                   + ld_107[k];

        t_176[k] = ab_y[k] * kd_105[k]
                   + ld_141[k];

        t_177[k] = ab_y[k] * kd_106[k]
                   + ld_142[k];

        t_178[k] = ab_y[k] * kd_107[k]
                   + ld_143[k];

        t_179[k] = ab_z[k] * kd_107[k]
                   + ld_149[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, kd_108, kd_109, kd_110, \
                         kd_111, kd_112, ld_108, ld_109, ld_110, ld_111, \
                         ld_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = ab_x[k] * kd_108[k]
                   + ld_108[k];

        t_181[k] = ab_x[k] * kd_109[k]
                   + ld_109[k];

        t_182[k] = ab_x[k] * kd_110[k]
                   + ld_110[k];

        t_183[k] = ab_x[k] * kd_111[k]
                   + ld_111[k];

        t_184[k] = ab_x[k] * kd_112[k]
                   + ld_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, kd_111, kd_112, \
                         kd_113, ld_113, ld_147, ld_148, ld_149, \
                         ld_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = ab_x[k] * kd_113[k]
                   + ld_113[k];

        t_186[k] = ab_y[k] * kd_111[k]
                   + ld_147[k];

        t_187[k] = ab_y[k] * kd_112[k]
                   + ld_148[k];

        t_188[k] = ab_y[k] * kd_113[k]
                   + ld_149[k];

        t_189[k] = ab_z[k] * kd_113[k]
                   + ld_155[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, kd_114, kd_115, kd_116, \
                         kd_117, kd_118, ld_114, ld_115, ld_116, ld_117, \
                         ld_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = ab_x[k] * kd_114[k]
                   + ld_114[k];

        t_191[k] = ab_x[k] * kd_115[k]
                   + ld_115[k];

        t_192[k] = ab_x[k] * kd_116[k]
                   + ld_116[k];

        t_193[k] = ab_x[k] * kd_117[k]
                   + ld_117[k];

        t_194[k] = ab_x[k] * kd_118[k]
                   + ld_118[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, kd_117, kd_118, \
                         kd_119, ld_119, ld_153, ld_154, ld_155, \
                         ld_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = ab_x[k] * kd_119[k]
                   + ld_119[k];

        t_196[k] = ab_y[k] * kd_117[k]
                   + ld_153[k];

        t_197[k] = ab_y[k] * kd_118[k]
                   + ld_154[k];

        t_198[k] = ab_y[k] * kd_119[k]
                   + ld_155[k];

        t_199[k] = ab_z[k] * kd_119[k]
                   + ld_161[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, kd_120, kd_121, kd_122, \
                         kd_123, kd_124, ld_120, ld_121, ld_122, ld_123, \
                         ld_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = ab_x[k] * kd_120[k]
                   + ld_120[k];

        t_201[k] = ab_x[k] * kd_121[k]
                   + ld_121[k];

        t_202[k] = ab_x[k] * kd_122[k]
                   + ld_122[k];

        t_203[k] = ab_x[k] * kd_123[k]
                   + ld_123[k];

        t_204[k] = ab_x[k] * kd_124[k]
                   + ld_124[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, kd_123, kd_124, \
                         kd_125, ld_125, ld_159, ld_160, ld_161, \
                         ld_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = ab_x[k] * kd_125[k]
                   + ld_125[k];

        t_206[k] = ab_y[k] * kd_123[k]
                   + ld_159[k];

        t_207[k] = ab_y[k] * kd_124[k]
                   + ld_160[k];

        t_208[k] = ab_y[k] * kd_125[k]
                   + ld_161[k];

        t_209[k] = ab_z[k] * kd_125[k]
                   + ld_167[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, kd_126, kd_127, kd_128, \
                         kd_129, kd_130, ld_126, ld_127, ld_128, ld_129, \
                         ld_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = ab_x[k] * kd_126[k]
                   + ld_126[k];

        t_211[k] = ab_x[k] * kd_127[k]
                   + ld_127[k];

        t_212[k] = ab_x[k] * kd_128[k]
                   + ld_128[k];

        t_213[k] = ab_x[k] * kd_129[k]
                   + ld_129[k];

        t_214[k] = ab_x[k] * kd_130[k]
                   + ld_130[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ab_y, ab_z, kd_129, kd_130, \
                         kd_131, ld_131, ld_171, ld_172, ld_173, \
                         ld_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = ab_x[k] * kd_131[k]
                   + ld_131[k];

        t_216[k] = ab_y[k] * kd_129[k]
                   + ld_171[k];

        t_217[k] = ab_y[k] * kd_130[k]
                   + ld_172[k];

        t_218[k] = ab_y[k] * kd_131[k]
                   + ld_173[k];

        t_219[k] = ab_z[k] * kd_131[k]
                   + ld_179[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, kd_132, kd_133, kd_134, \
                         kd_135, kd_136, ld_132, ld_133, ld_134, ld_135, \
                         ld_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = ab_x[k] * kd_132[k]
                   + ld_132[k];

        t_221[k] = ab_x[k] * kd_133[k]
                   + ld_133[k];

        t_222[k] = ab_x[k] * kd_134[k]
                   + ld_134[k];

        t_223[k] = ab_x[k] * kd_135[k]
                   + ld_135[k];

        t_224[k] = ab_x[k] * kd_136[k]
                   + ld_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, ab_y, ab_z, kd_135, kd_136, \
                         kd_137, ld_137, ld_177, ld_178, ld_179, \
                         ld_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = ab_x[k] * kd_137[k]
                   + ld_137[k];

        t_226[k] = ab_y[k] * kd_135[k]
                   + ld_177[k];

        t_227[k] = ab_y[k] * kd_136[k]
                   + ld_178[k];

        t_228[k] = ab_y[k] * kd_137[k]
                   + ld_179[k];

        t_229[k] = ab_z[k] * kd_137[k]
                   + ld_185[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, kd_138, kd_139, kd_140, \
                         kd_141, kd_142, ld_138, ld_139, ld_140, ld_141, \
                         ld_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = ab_x[k] * kd_138[k]
                   + ld_138[k];

        t_231[k] = ab_x[k] * kd_139[k]
                   + ld_139[k];

        t_232[k] = ab_x[k] * kd_140[k]
                   + ld_140[k];

        t_233[k] = ab_x[k] * kd_141[k]
                   + ld_141[k];

        t_234[k] = ab_x[k] * kd_142[k]
                   + ld_142[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, ab_y, ab_z, kd_141, kd_142, \
                         kd_143, ld_143, ld_183, ld_184, ld_185, \
                         ld_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = ab_x[k] * kd_143[k]
                   + ld_143[k];

        t_236[k] = ab_y[k] * kd_141[k]
                   + ld_183[k];

        t_237[k] = ab_y[k] * kd_142[k]
                   + ld_184[k];

        t_238[k] = ab_y[k] * kd_143[k]
                   + ld_185[k];

        t_239[k] = ab_z[k] * kd_143[k]
                   + ld_191[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, kd_144, kd_145, kd_146, \
                         kd_147, kd_148, ld_144, ld_145, ld_146, ld_147, \
                         ld_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = ab_x[k] * kd_144[k]
                   + ld_144[k];

        t_241[k] = ab_x[k] * kd_145[k]
                   + ld_145[k];

        t_242[k] = ab_x[k] * kd_146[k]
                   + ld_146[k];

        t_243[k] = ab_x[k] * kd_147[k]
                   + ld_147[k];

        t_244[k] = ab_x[k] * kd_148[k]
                   + ld_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, ab_y, ab_z, kd_147, kd_148, \
                         kd_149, ld_149, ld_189, ld_190, ld_191, \
                         ld_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = ab_x[k] * kd_149[k]
                   + ld_149[k];

        t_246[k] = ab_y[k] * kd_147[k]
                   + ld_189[k];

        t_247[k] = ab_y[k] * kd_148[k]
                   + ld_190[k];

        t_248[k] = ab_y[k] * kd_149[k]
                   + ld_191[k];

        t_249[k] = ab_z[k] * kd_149[k]
                   + ld_197[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, kd_150, kd_151, kd_152, \
                         kd_153, kd_154, ld_150, ld_151, ld_152, ld_153, \
                         ld_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = ab_x[k] * kd_150[k]
                   + ld_150[k];

        t_251[k] = ab_x[k] * kd_151[k]
                   + ld_151[k];

        t_252[k] = ab_x[k] * kd_152[k]
                   + ld_152[k];

        t_253[k] = ab_x[k] * kd_153[k]
                   + ld_153[k];

        t_254[k] = ab_x[k] * kd_154[k]
                   + ld_154[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, ab_y, ab_z, kd_153, kd_154, \
                         kd_155, ld_155, ld_195, ld_196, ld_197, \
                         ld_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = ab_x[k] * kd_155[k]
                   + ld_155[k];

        t_256[k] = ab_y[k] * kd_153[k]
                   + ld_195[k];

        t_257[k] = ab_y[k] * kd_154[k]
                   + ld_196[k];

        t_258[k] = ab_y[k] * kd_155[k]
                   + ld_197[k];

        t_259[k] = ab_z[k] * kd_155[k]
                   + ld_203[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, kd_156, kd_157, kd_158, \
                         kd_159, kd_160, ld_156, ld_157, ld_158, ld_159, \
                         ld_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = ab_x[k] * kd_156[k]
                   + ld_156[k];

        t_261[k] = ab_x[k] * kd_157[k]
                   + ld_157[k];

        t_262[k] = ab_x[k] * kd_158[k]
                   + ld_158[k];

        t_263[k] = ab_x[k] * kd_159[k]
                   + ld_159[k];

        t_264[k] = ab_x[k] * kd_160[k]
                   + ld_160[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, kd_159, kd_160, \
                         kd_161, ld_161, ld_201, ld_202, ld_203, \
                         ld_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = ab_x[k] * kd_161[k]
                   + ld_161[k];

        t_266[k] = ab_y[k] * kd_159[k]
                   + ld_201[k];

        t_267[k] = ab_y[k] * kd_160[k]
                   + ld_202[k];

        t_268[k] = ab_y[k] * kd_161[k]
                   + ld_203[k];

        t_269[k] = ab_z[k] * kd_161[k]
                   + ld_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, kd_162, kd_163, kd_164, \
                         kd_165, kd_166, ld_162, ld_163, ld_164, ld_165, \
                         ld_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = ab_x[k] * kd_162[k]
                   + ld_162[k];

        t_271[k] = ab_x[k] * kd_163[k]
                   + ld_163[k];

        t_272[k] = ab_x[k] * kd_164[k]
                   + ld_164[k];

        t_273[k] = ab_x[k] * kd_165[k]
                   + ld_165[k];

        t_274[k] = ab_x[k] * kd_166[k]
                   + ld_166[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, ab_y, ab_z, kd_165, kd_166, \
                         kd_167, ld_167, ld_207, ld_208, ld_209, \
                         ld_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = ab_x[k] * kd_167[k]
                   + ld_167[k];

        t_276[k] = ab_y[k] * kd_165[k]
                   + ld_207[k];

        t_277[k] = ab_y[k] * kd_166[k]
                   + ld_208[k];

        t_278[k] = ab_y[k] * kd_167[k]
                   + ld_209[k];

        t_279[k] = ab_z[k] * kd_167[k]
                   + ld_215[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, kd_168, kd_169, kd_170, \
                         kd_171, kd_172, ld_168, ld_169, ld_170, ld_171, \
                         ld_172 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_280[k] = ab_x[k] * kd_168[k]
                   + ld_168[k];

        t_281[k] = ab_x[k] * kd_169[k]
                   + ld_169[k];

        t_282[k] = ab_x[k] * kd_170[k]
                   + ld_170[k];

        t_283[k] = ab_x[k] * kd_171[k]
                   + ld_171[k];

        t_284[k] = ab_x[k] * kd_172[k]
                   + ld_172[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, ab_y, ab_z, kd_171, kd_172, \
                         kd_173, ld_173, ld_219, ld_220, ld_221, \
                         ld_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_285[k] = ab_x[k] * kd_173[k]
                   + ld_173[k];

        t_286[k] = ab_y[k] * kd_171[k]
                   + ld_219[k];

        t_287[k] = ab_y[k] * kd_172[k]
                   + ld_220[k];

        t_288[k] = ab_y[k] * kd_173[k]
                   + ld_221[k];

        t_289[k] = ab_z[k] * kd_173[k]
                   + ld_227[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, kd_174, kd_175, kd_176, \
                         kd_177, kd_178, ld_174, ld_175, ld_176, ld_177, \
                         ld_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_290[k] = ab_x[k] * kd_174[k]
                   + ld_174[k];

        t_291[k] = ab_x[k] * kd_175[k]
                   + ld_175[k];

        t_292[k] = ab_x[k] * kd_176[k]
                   + ld_176[k];

        t_293[k] = ab_x[k] * kd_177[k]
                   + ld_177[k];

        t_294[k] = ab_x[k] * kd_178[k]
                   + ld_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, ab_y, ab_z, kd_177, kd_178, \
                         kd_179, ld_179, ld_225, ld_226, ld_227, \
                         ld_233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_295[k] = ab_x[k] * kd_179[k]
                   + ld_179[k];

        t_296[k] = ab_y[k] * kd_177[k]
                   + ld_225[k];

        t_297[k] = ab_y[k] * kd_178[k]
                   + ld_226[k];

        t_298[k] = ab_y[k] * kd_179[k]
                   + ld_227[k];

        t_299[k] = ab_z[k] * kd_179[k]
                   + ld_233[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, kd_180, kd_181, kd_182, \
                         kd_183, kd_184, ld_180, ld_181, ld_182, ld_183, \
                         ld_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_300[k] = ab_x[k] * kd_180[k]
                   + ld_180[k];

        t_301[k] = ab_x[k] * kd_181[k]
                   + ld_181[k];

        t_302[k] = ab_x[k] * kd_182[k]
                   + ld_182[k];

        t_303[k] = ab_x[k] * kd_183[k]
                   + ld_183[k];

        t_304[k] = ab_x[k] * kd_184[k]
                   + ld_184[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, ab_y, ab_z, kd_183, kd_184, \
                         kd_185, ld_185, ld_231, ld_232, ld_233, \
                         ld_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_305[k] = ab_x[k] * kd_185[k]
                   + ld_185[k];

        t_306[k] = ab_y[k] * kd_183[k]
                   + ld_231[k];

        t_307[k] = ab_y[k] * kd_184[k]
                   + ld_232[k];

        t_308[k] = ab_y[k] * kd_185[k]
                   + ld_233[k];

        t_309[k] = ab_z[k] * kd_185[k]
                   + ld_239[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, kd_186, kd_187, kd_188, \
                         kd_189, kd_190, ld_186, ld_187, ld_188, ld_189, \
                         ld_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_310[k] = ab_x[k] * kd_186[k]
                   + ld_186[k];

        t_311[k] = ab_x[k] * kd_187[k]
                   + ld_187[k];

        t_312[k] = ab_x[k] * kd_188[k]
                   + ld_188[k];

        t_313[k] = ab_x[k] * kd_189[k]
                   + ld_189[k];

        t_314[k] = ab_x[k] * kd_190[k]
                   + ld_190[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, ab_y, ab_z, kd_189, kd_190, \
                         kd_191, ld_191, ld_237, ld_238, ld_239, \
                         ld_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_315[k] = ab_x[k] * kd_191[k]
                   + ld_191[k];

        t_316[k] = ab_y[k] * kd_189[k]
                   + ld_237[k];

        t_317[k] = ab_y[k] * kd_190[k]
                   + ld_238[k];

        t_318[k] = ab_y[k] * kd_191[k]
                   + ld_239[k];

        t_319[k] = ab_z[k] * kd_191[k]
                   + ld_245[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, kd_192, kd_193, kd_194, \
                         kd_195, kd_196, ld_192, ld_193, ld_194, ld_195, \
                         ld_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_320[k] = ab_x[k] * kd_192[k]
                   + ld_192[k];

        t_321[k] = ab_x[k] * kd_193[k]
                   + ld_193[k];

        t_322[k] = ab_x[k] * kd_194[k]
                   + ld_194[k];

        t_323[k] = ab_x[k] * kd_195[k]
                   + ld_195[k];

        t_324[k] = ab_x[k] * kd_196[k]
                   + ld_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, ab_y, ab_z, kd_195, kd_196, \
                         kd_197, ld_197, ld_243, ld_244, ld_245, \
                         ld_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_325[k] = ab_x[k] * kd_197[k]
                   + ld_197[k];

        t_326[k] = ab_y[k] * kd_195[k]
                   + ld_243[k];

        t_327[k] = ab_y[k] * kd_196[k]
                   + ld_244[k];

        t_328[k] = ab_y[k] * kd_197[k]
                   + ld_245[k];

        t_329[k] = ab_z[k] * kd_197[k]
                   + ld_251[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, kd_198, kd_199, kd_200, \
                         kd_201, kd_202, ld_198, ld_199, ld_200, ld_201, \
                         ld_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_330[k] = ab_x[k] * kd_198[k]
                   + ld_198[k];

        t_331[k] = ab_x[k] * kd_199[k]
                   + ld_199[k];

        t_332[k] = ab_x[k] * kd_200[k]
                   + ld_200[k];

        t_333[k] = ab_x[k] * kd_201[k]
                   + ld_201[k];

        t_334[k] = ab_x[k] * kd_202[k]
                   + ld_202[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, ab_y, ab_z, kd_201, kd_202, \
                         kd_203, ld_203, ld_249, ld_250, ld_251, \
                         ld_257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_335[k] = ab_x[k] * kd_203[k]
                   + ld_203[k];

        t_336[k] = ab_y[k] * kd_201[k]
                   + ld_249[k];

        t_337[k] = ab_y[k] * kd_202[k]
                   + ld_250[k];

        t_338[k] = ab_y[k] * kd_203[k]
                   + ld_251[k];

        t_339[k] = ab_z[k] * kd_203[k]
                   + ld_257[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, kd_204, kd_205, kd_206, \
                         kd_207, kd_208, ld_204, ld_205, ld_206, ld_207, \
                         ld_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_340[k] = ab_x[k] * kd_204[k]
                   + ld_204[k];

        t_341[k] = ab_x[k] * kd_205[k]
                   + ld_205[k];

        t_342[k] = ab_x[k] * kd_206[k]
                   + ld_206[k];

        t_343[k] = ab_x[k] * kd_207[k]
                   + ld_207[k];

        t_344[k] = ab_x[k] * kd_208[k]
                   + ld_208[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, ab_y, ab_z, kd_207, kd_208, \
                         kd_209, ld_209, ld_255, ld_256, ld_257, \
                         ld_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_345[k] = ab_x[k] * kd_209[k]
                   + ld_209[k];

        t_346[k] = ab_y[k] * kd_207[k]
                   + ld_255[k];

        t_347[k] = ab_y[k] * kd_208[k]
                   + ld_256[k];

        t_348[k] = ab_y[k] * kd_209[k]
                   + ld_257[k];

        t_349[k] = ab_z[k] * kd_209[k]
                   + ld_263[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, kd_210, kd_211, kd_212, \
                         kd_213, kd_214, ld_210, ld_211, ld_212, ld_213, \
                         ld_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_350[k] = ab_x[k] * kd_210[k]
                   + ld_210[k];

        t_351[k] = ab_x[k] * kd_211[k]
                   + ld_211[k];

        t_352[k] = ab_x[k] * kd_212[k]
                   + ld_212[k];

        t_353[k] = ab_x[k] * kd_213[k]
                   + ld_213[k];

        t_354[k] = ab_x[k] * kd_214[k]
                   + ld_214[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, ab_y, ab_z, kd_213, kd_214, \
                         kd_215, ld_215, ld_261, ld_262, ld_263, \
                         ld_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_355[k] = ab_x[k] * kd_215[k]
                   + ld_215[k];

        t_356[k] = ab_y[k] * kd_213[k]
                   + ld_261[k];

        t_357[k] = ab_y[k] * kd_214[k]
                   + ld_262[k];

        t_358[k] = ab_y[k] * kd_215[k]
                   + ld_263[k];

        t_359[k] = ab_z[k] * kd_215[k]
                   + ld_269[k];
    }
}

}  // namespace simdtrf
