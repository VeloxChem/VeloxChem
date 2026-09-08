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


#include "SimdElectronRepulsionVrrRecFK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_fk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 2.0 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_13 = buffer.data(fh0 + 13);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_15 = buffer.data(fh0 + 15);
    const auto *fh0_16 = buffer.data(fh0 + 16);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_18 = buffer.data(fh0 + 18);
    const auto *fh0_19 = buffer.data(fh0 + 19);
    const auto *fh0_20 = buffer.data(fh0 + 20);
    const auto *fh0_21 = buffer.data(fh0 + 21);
    const auto *fh0_22 = buffer.data(fh0 + 22);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_36 = buffer.data(fh0 + 36);
    const auto *fh0_37 = buffer.data(fh0 + 37);
    const auto *fh0_38 = buffer.data(fh0 + 38);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_13 = buffer.data(fh1 + 13);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_15 = buffer.data(fh1 + 15);
    const auto *fh1_16 = buffer.data(fh1 + 16);
    const auto *fh1_17 = buffer.data(fh1 + 17);
    const auto *fh1_18 = buffer.data(fh1 + 18);
    const auto *fh1_19 = buffer.data(fh1 + 19);
    const auto *fh1_20 = buffer.data(fh1 + 20);
    const auto *fh1_21 = buffer.data(fh1 + 21);
    const auto *fh1_22 = buffer.data(fh1 + 22);
    const auto *fh1_23 = buffer.data(fh1 + 23);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_36 = buffer.data(fh1 + 36);
    const auto *fh1_37 = buffer.data(fh1 + 37);
    const auto *fh1_38 = buffer.data(fh1 + 38);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_112 = buffer.data(fi + 112);
    const auto *fi_113 = buffer.data(fi + 113);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_116 = buffer.data(fi + 116);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_119 = buffer.data(fi + 119);
    const auto *fi_120 = buffer.data(fi + 120);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_123 = buffer.data(fi + 123);
    const auto *fi_124 = buffer.data(fi + 124);
    const auto *fi_125 = buffer.data(fi + 125);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_127 = buffer.data(fi + 127);
    const auto *fi_128 = buffer.data(fi + 128);
    const auto *fi_129 = buffer.data(fi + 129);
    const auto *fi_130 = buffer.data(fi + 130);
    const auto *fi_131 = buffer.data(fi + 131);
    const auto *fi_132 = buffer.data(fi + 132);
    const auto *fi_133 = buffer.data(fi + 133);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_139 = buffer.data(fi + 139);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_141 = buffer.data(fi + 141);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_144 = buffer.data(fi + 144);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_147 = buffer.data(fi + 147);
    const auto *fi_148 = buffer.data(fi + 148);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_151 = buffer.data(fi + 151);
    const auto *fi_152 = buffer.data(fi + 152);
    const auto *fi_153 = buffer.data(fi + 153);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_155 = buffer.data(fi + 155);
    const auto *fi_156 = buffer.data(fi + 156);
    const auto *fi_157 = buffer.data(fi + 157);
    const auto *fi_158 = buffer.data(fi + 158);
    const auto *fi_159 = buffer.data(fi + 159);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, di_0, fh0_0, fh1_0, \
                         fi_0, fi_1, fi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pb_y[k] * fi_0[k];

        t_2[k] = pb_z[k] * fi_0[k];

        t_3[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_4[k] = pb_y[k] * fi_2[k];

        t_5[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, fh0_1, fh0_2, fh0_3, fh1_1, \
                         fh1_2, fh1_3, fi_3, fi_4, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_7[k] = pb_z[k] * fi_3[k];

        t_8[k] = pb_y[k] * fi_4[k];

        t_9[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_10[k] = f_7 * fh0_3[k]
                  - f_8 * fh1_3[k]
                  + pb_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, fh0_4, fh0_5, fh1_4, \
                         fh1_5, fi_5, fi_6, fi_7, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * fi_5[k];

        t_12[k] = f_3 * fh0_4[k]
                  - f_4 * fh1_4[k]
                  + pb_y[k] * fi_6[k];

        t_13[k] = pb_y[k] * fi_7[k];

        t_14[k] = f_7 * fh0_4[k]
                  - f_8 * fh1_4[k]
                  + pb_z[k] * fi_7[k];

        t_15[k] = f_9 * fh0_5[k]
                  - f_10 * fh1_5[k]
                  + pb_y[k] * fi_8[k];

        t_16[k] = pb_z[k] * fi_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, fh0_6, fh0_7, fh1_6, fh1_7, fi_9, \
                         fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * fh0_6[k]
                  - f_6 * fh1_6[k]
                  + pb_y[k] * fi_9[k];

        t_18[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_19[k] = pb_y[k] * fi_11[k];

        t_20[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, di_14, di_16, di_17, di_18, \
                         fi_12, fi_14, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * di_14[k]
                  + pb_x[k] * fi_14[k];

        t_22[k] = pb_z[k] * fi_12[k];

        t_23[k] = f_0 * di_16[k]
                  + pb_x[k] * fi_15[k];

        t_24[k] = f_0 * di_17[k]
                  + pb_x[k] * fi_16[k];

        t_25[k] = f_0 * di_18[k]
                  + pb_x[k] * fi_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, di_20, fh0_8, fh1_8, fi_13, \
                         fi_14, fi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * fi_13[k];

        t_27[k] = f_0 * di_20[k]
                  + pb_x[k] * fi_19[k];

        t_28[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_14[k];

        t_29[k] = pb_z[k] * fi_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, fh0_9, fh0_10, fh0_11, fh1_9, fh1_10, fh1_11, \
                         fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_9[k]
                  + pb_y[k] * fi_15[k];

        t_31[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_10[k]
                  + pb_y[k] * fi_16[k];

        t_32[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_11[k]
                  + pb_y[k] * fi_17[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, di_0, dk_0, \
                         fh0_12, fh1_12, fi_18, fi_19, fi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_12[k]
                  + pb_y[k] * fi_18[k];

        t_34[k] = pb_y[k] * fi_19[k];

        t_35[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_12[k]
                  + pb_z[k] * fi_19[k];

        t_36[k] = pa_y[k] * dk_0[k];

        t_37[k] = f_11 * di_0[k]
                  + pb_y[k] * fi_20[k];

        t_38[k] = pb_z[k] * fi_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, di_1, di_3, dk_1, dk_2, \
                         dk_3, fi_21, fi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * di_1[k]
                  + pa_y[k] * dk_1[k];

        t_40[k] = pb_z[k] * fi_21[k];

        t_41[k] = pa_y[k] * dk_2[k];

        t_42[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_3[k];

        t_43[k] = pb_z[k] * fi_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, di_4, di_5, di_7, \
                         dk_4, dk_5, dk_6, fi_23, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * di_4[k]
                  + pb_y[k] * fi_23[k];

        t_45[k] = pa_y[k] * dk_4[k];

        t_46[k] = f_13 * di_5[k]
                  + pa_y[k] * dk_5[k];

        t_47[k] = pb_z[k] * fi_24[k];

        t_48[k] = f_12 * di_7[k]
                  + pa_y[k] * dk_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, di_8, di_9, di_11, \
                         dk_7, dk_8, dk_9, fi_25, fi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * di_8[k]
                  + pb_y[k] * fi_25[k];

        t_50[k] = pa_y[k] * dk_7[k];

        t_51[k] = f_14 * di_9[k]
                  + pa_y[k] * dk_8[k];

        t_52[k] = pb_z[k] * fi_26[k];

        t_53[k] = f_0 * di_11[k]
                  + pa_y[k] * dk_9[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, di_12, di_13, di_28, dk_10, \
                         dk_11, fi_27, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * di_12[k]
                  + pa_y[k] * dk_10[k];

        t_55[k] = f_11 * di_13[k]
                  + pb_y[k] * fi_27[k];

        t_56[k] = pa_y[k] * dk_11[k];

        t_57[k] = f_12 * di_28[k]
                  + pb_x[k] * fi_29[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, di_29, di_30, di_31, di_32, \
                         fi_28, fi_30, fi_31, fi_32, fi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * fi_28[k];

        t_59[k] = f_12 * di_29[k]
                  + pb_x[k] * fi_30[k];

        t_60[k] = f_12 * di_30[k]
                  + pb_x[k] * fi_31[k];

        t_61[k] = f_12 * di_31[k]
                  + pb_x[k] * fi_32[k];

        t_62[k] = f_12 * di_32[k]
                  + pb_x[k] * fi_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, di_14, di_16, di_17, dk_13, \
                         dk_14, dk_15, dk_16, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * dk_13[k];

        t_64[k] = f_15 * di_14[k]
                  + pa_y[k] * dk_14[k];

        t_65[k] = pb_z[k] * fi_29[k];

        t_66[k] = f_14 * di_16[k]
                  + pa_y[k] * dk_15[k];

        t_67[k] = f_13 * di_17[k]
                  + pa_y[k] * dk_16[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, di_18, di_19, di_20, \
                         dk_0, dk_17, dk_18, dk_19, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * di_18[k]
                  + pa_y[k] * dk_17[k];

        t_69[k] = f_12 * di_19[k]
                  + pa_y[k] * dk_18[k];

        t_70[k] = f_11 * di_20[k]
                  + pb_y[k] * fi_34[k];

        t_71[k] = pa_y[k] * dk_19[k];

        t_72[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, di_0, di_2, \
                         dk_1, dk_2, dk_3, fi_35, fi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * fi_35[k];

        t_74[k] = f_11 * di_0[k]
                  + pb_z[k] * fi_35[k];

        t_75[k] = pa_z[k] * dk_1[k];

        t_76[k] = pb_y[k] * fi_36[k];

        t_77[k] = f_12 * di_2[k]
                  + pa_z[k] * dk_2[k];

        t_78[k] = pa_z[k] * dk_3[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, di_3, di_4, di_5, \
                         dk_4, dk_5, fi_37, fi_38, fi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * di_3[k]
                  + pb_z[k] * fi_37[k];

        t_80[k] = pb_y[k] * fi_38[k];

        t_81[k] = f_0 * di_4[k]
                  + pa_z[k] * dk_4[k];

        t_82[k] = pa_z[k] * dk_5[k];

        t_83[k] = f_11 * di_5[k]
                  + pb_z[k] * fi_39[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, di_6, di_8, di_9, \
                         dk_6, dk_7, dk_8, fi_40, fi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * di_6[k]
                  + pa_z[k] * dk_6[k];

        t_85[k] = pb_y[k] * fi_40[k];

        t_86[k] = f_13 * di_8[k]
                  + pa_z[k] * dk_7[k];

        t_87[k] = pa_z[k] * dk_8[k];

        t_88[k] = f_11 * di_9[k]
                  + pb_z[k] * fi_41[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, di_10, di_11, di_13, dk_9, \
                         dk_10, dk_11, dk_12, fi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * di_10[k]
                  + pa_z[k] * dk_9[k];

        t_90[k] = f_0 * di_11[k]
                  + pa_z[k] * dk_10[k];

        t_91[k] = pb_y[k] * fi_42[k];

        t_92[k] = f_14 * di_13[k]
                  + pa_z[k] * dk_11[k];

        t_93[k] = pa_z[k] * dk_12[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, di_41, di_42, di_43, di_44, \
                         fi_43, fi_45, fi_46, fi_47, fi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_12 * di_41[k]
                  + pb_x[k] * fi_45[k];

        t_95[k] = f_12 * di_42[k]
                  + pb_x[k] * fi_46[k];

        t_96[k] = f_12 * di_43[k]
                  + pb_x[k] * fi_47[k];

        t_97[k] = f_12 * di_44[k]
                  + pb_x[k] * fi_48[k];

        t_98[k] = pb_y[k] * fi_43[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, di_14, di_15, di_45, \
                         dk_14, dk_15, fi_44, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_12 * di_45[k]
                  + pb_x[k] * fi_49[k];

        t_100[k] = pa_z[k] * dk_14[k];

        t_101[k] = f_11 * di_14[k]
                   + pb_z[k] * fi_44[k];

        t_102[k] = f_12 * di_15[k]
                   + pa_z[k] * dk_15[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, di_16, di_17, di_18, \
                         di_20, dk_16, dk_17, dk_18, dk_19, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * di_16[k]
                   + pa_z[k] * dk_16[k];

        t_104[k] = f_13 * di_17[k]
                   + pa_z[k] * dk_17[k];

        t_105[k] = f_14 * di_18[k]
                   + pa_z[k] * dk_18[k];

        t_106[k] = pb_y[k] * fi_49[k];

        t_107[k] = f_15 * di_20[k]
                   + pa_z[k] * dk_19[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_x, pb_y, pb_z, di_21, di_46, \
                         di_48, dk_33, dk_35, fi_50, fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * di_46[k]
                   + pa_x[k] * dk_33[k];

        t_109[k] = f_12 * di_21[k]
                   + pb_y[k] * fi_50[k];

        t_110[k] = pb_z[k] * fi_50[k];

        t_111[k] = f_14 * di_48[k]
                   + pa_x[k] * dk_35[k];

        t_112[k] = pb_z[k] * fi_51[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pb_y, pb_z, di_23, di_49, di_50, \
                         dk_36, dk_37, fi_52, fi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * di_49[k]
                   + pa_x[k] * dk_36[k];

        t_114[k] = f_13 * di_50[k]
                   + pa_x[k] * dk_37[k];

        t_115[k] = pb_z[k] * fi_52[k];

        t_116[k] = f_12 * di_23[k]
                   + pb_y[k] * fi_53[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pb_z, di_52, di_53, di_55, dk_38, \
                         dk_39, dk_40, fi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_13 * di_52[k]
                   + pa_x[k] * dk_38[k];

        t_118[k] = f_0 * di_53[k]
                   + pa_x[k] * dk_39[k];

        t_119[k] = pb_z[k] * fi_54[k];

        t_120[k] = f_0 * di_55[k]
                   + pa_x[k] * dk_40[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pb_y, pb_z, di_25, di_56, di_57, \
                         dk_41, dk_42, fi_55, fi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * di_25[k]
                   + pb_y[k] * fi_55[k];

        t_122[k] = f_0 * di_56[k]
                   + pa_x[k] * dk_41[k];

        t_123[k] = f_12 * di_57[k]
                   + pa_x[k] * dk_42[k];

        t_124[k] = pb_z[k] * fi_56[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_x, pb_y, di_27, di_58, di_59, di_60, \
                         dk_43, dk_44, dk_45, fi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * di_58[k]
                   + pa_x[k] * dk_43[k];

        t_126[k] = f_12 * di_59[k]
                   + pa_x[k] * dk_44[k];

        t_127[k] = f_12 * di_27[k]
                   + pb_y[k] * fi_57[k];

        t_128[k] = f_12 * di_60[k]
                   + pa_x[k] * dk_45[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, di_61, di_63, di_64, \
                         di_65, fi_58, fi_59, fi_60, fi_61, fi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * di_61[k]
                   + pb_x[k] * fi_59[k];

        t_130[k] = pb_z[k] * fi_58[k];

        t_131[k] = f_11 * di_63[k]
                   + pb_x[k] * fi_60[k];

        t_132[k] = f_11 * di_64[k]
                   + pb_x[k] * fi_61[k];

        t_133[k] = f_11 * di_65[k]
                   + pb_x[k] * fi_62[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_x, pb_z, di_66, di_67, \
                         dk_46, dk_47, fi_59, fi_63, fi_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * di_66[k]
                   + pb_x[k] * fi_63[k];

        t_135[k] = f_11 * di_67[k]
                   + pb_x[k] * fi_64[k];

        t_136[k] = pa_x[k] * dk_46[k];

        t_137[k] = pb_z[k] * fi_59[k];

        t_138[k] = pa_x[k] * dk_47[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, pa_x, pa_y, dk_26, dk_48, \
                         dk_49, dk_50, dk_51, dk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_x[k] * dk_48[k];

        t_140[k] = pa_x[k] * dk_49[k];

        t_141[k] = pa_x[k] * dk_50[k];

        t_142[k] = pa_x[k] * dk_51[k];

        t_143[k] = pa_x[k] * dk_52[k];

        t_144[k] = pa_y[k] * dk_26[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, di_34, \
                         dk_20, dk_21, dk_22, dk_27, dk_28, fi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * dk_20[k];

        t_146[k] = pa_y[k] * dk_27[k];

        t_147[k] = pa_z[k] * dk_21[k];

        t_148[k] = f_11 * di_34[k]
                   + pb_y[k] * fi_65[k];

        t_149[k] = pa_y[k] * dk_28[k];

        t_150[k] = pa_z[k] * dk_22[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, di_22, di_36, \
                         dk_23, dk_29, fi_66, fi_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * di_22[k]
                   + pb_z[k] * fi_66[k];

        t_152[k] = f_11 * di_36[k]
                   + pb_y[k] * fi_67[k];

        t_153[k] = pa_y[k] * dk_29[k];

        t_154[k] = pa_z[k] * dk_23[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_x, pa_y, pb_y, pb_z, di_24, di_38, \
                         di_74, dk_30, dk_53, fi_68, fi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * di_24[k]
                   + pb_z[k] * fi_68[k];

        t_156[k] = f_0 * di_74[k]
                   + pa_x[k] * dk_53[k];

        t_157[k] = f_11 * di_38[k]
                   + pb_y[k] * fi_69[k];

        t_158[k] = pa_y[k] * dk_30[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pa_z, pb_z, di_26, di_76, di_77, \
                         dk_24, dk_54, dk_55, fi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * dk_24[k];

        t_160[k] = f_11 * di_26[k]
                   + pb_z[k] * fi_70[k];

        t_161[k] = f_12 * di_76[k]
                   + pa_x[k] * dk_54[k];

        t_162[k] = f_12 * di_77[k]
                   + pa_x[k] * dk_55[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, di_40, di_79, \
                         dk_25, dk_31, fi_71, fi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * di_40[k]
                   + pb_y[k] * fi_71[k];

        t_164[k] = pa_y[k] * dk_31[k];

        t_165[k] = pa_z[k] * dk_25[k];

        t_166[k] = f_11 * di_79[k]
                   + pb_x[k] * fi_72[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, di_80, di_81, di_82, \
                         di_83, dk_32, fi_73, fi_74, fi_75, fi_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_11 * di_80[k]
                   + pb_x[k] * fi_73[k];

        t_168[k] = f_11 * di_81[k]
                   + pb_x[k] * fi_74[k];

        t_169[k] = f_11 * di_82[k]
                   + pb_x[k] * fi_75[k];

        t_170[k] = f_11 * di_83[k]
                   + pb_x[k] * fi_76[k];

        t_171[k] = pa_y[k] * dk_32[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, t_177, t_178, pa_x, dk_56, dk_57, \
                         dk_58, dk_59, dk_60, dk_61, dk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_x[k] * dk_56[k];

        t_173[k] = pa_x[k] * dk_57[k];

        t_174[k] = pa_x[k] * dk_58[k];

        t_175[k] = pa_x[k] * dk_59[k];

        t_176[k] = pa_x[k] * dk_60[k];

        t_177[k] = pa_x[k] * dk_61[k];

        t_178[k] = pa_x[k] * dk_62[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pa_x, pb_y, pb_z, di_33, di_85, \
                         di_88, dk_63, dk_64, dk_66, fi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_x[k] * dk_63[k];

        t_180[k] = f_15 * di_85[k]
                   + pa_x[k] * dk_64[k];

        t_181[k] = pb_y[k] * fi_77[k];

        t_182[k] = f_12 * di_33[k]
                   + pb_z[k] * fi_77[k];

        t_183[k] = f_14 * di_88[k]
                   + pa_x[k] * dk_66[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pa_x, pb_y, pb_z, di_35, di_89, \
                         di_90, dk_67, dk_68, fi_78, fi_79, fi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * fi_78[k];

        t_185[k] = f_14 * di_89[k]
                   + pa_x[k] * dk_67[k];

        t_186[k] = f_13 * di_90[k]
                   + pa_x[k] * dk_68[k];

        t_187[k] = f_12 * di_35[k]
                   + pb_z[k] * fi_79[k];

        t_188[k] = pb_y[k] * fi_80[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pb_z, di_37, di_92, di_93, di_94, \
                         dk_69, dk_70, dk_71, fi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * di_92[k]
                   + pa_x[k] * dk_69[k];

        t_190[k] = f_0 * di_93[k]
                   + pa_x[k] * dk_70[k];

        t_191[k] = f_12 * di_37[k]
                   + pb_z[k] * fi_81[k];

        t_192[k] = f_0 * di_94[k]
                   + pa_x[k] * dk_71[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_x, pb_y, pb_z, di_39, di_96, di_97, \
                         dk_72, dk_73, fi_82, fi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_y[k] * fi_82[k];

        t_194[k] = f_0 * di_96[k]
                   + pa_x[k] * dk_72[k];

        t_195[k] = f_12 * di_97[k]
                   + pa_x[k] * dk_73[k];

        t_196[k] = f_12 * di_39[k]
                   + pb_z[k] * fi_83[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_x, pb_y, di_98, di_99, di_100, dk_74, \
                         dk_75, dk_76, fi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * di_98[k]
                   + pa_x[k] * dk_74[k];

        t_198[k] = f_12 * di_99[k]
                   + pa_x[k] * dk_75[k];

        t_199[k] = pb_y[k] * fi_84[k];

        t_200[k] = f_12 * di_100[k]
                   + pa_x[k] * dk_76[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pb_x, di_101, di_102, di_103, \
                         di_104, di_105, fi_86, fi_87, fi_88, fi_89, \
                         fi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_11 * di_101[k]
                   + pb_x[k] * fi_86[k];

        t_202[k] = f_11 * di_102[k]
                   + pb_x[k] * fi_87[k];

        t_203[k] = f_11 * di_103[k]
                   + pb_x[k] * fi_88[k];

        t_204[k] = f_11 * di_104[k]
                   + pb_x[k] * fi_89[k];

        t_205[k] = f_11 * di_105[k]
                   + pb_x[k] * fi_90[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, t_211, pa_x, pb_x, pb_y, di_107, \
                         dk_77, dk_78, dk_79, dk_80, fi_85, fi_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_y[k] * fi_85[k];

        t_207[k] = f_11 * di_107[k]
                   + pb_x[k] * fi_91[k];

        t_208[k] = pa_x[k] * dk_77[k];

        t_209[k] = pa_x[k] * dk_78[k];

        t_210[k] = pa_x[k] * dk_79[k];

        t_211[k] = pa_x[k] * dk_80[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pb_x, pb_y, dk_81, dk_82, \
                         dk_83, fh0_13, fh1_13, fi_91, fi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pa_x[k] * dk_81[k];

        t_213[k] = pa_x[k] * dk_82[k];

        t_214[k] = pb_y[k] * fi_91[k];

        t_215[k] = pa_x[k] * dk_83[k];

        t_216[k] = f_1 * fh0_13[k]
                   - f_2 * fh1_13[k]
                   + pb_x[k] * fi_92[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_x, pb_y, pb_z, di_46, fh0_14, fh1_14, \
                         fi_92, fi_93, fi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_0 * di_46[k]
                   + pb_y[k] * fi_92[k];

        t_218[k] = pb_z[k] * fi_92[k];

        t_219[k] = f_9 * fh0_14[k]
                   - f_10 * fh1_14[k]
                   + pb_x[k] * fi_94[k];

        t_220[k] = pb_z[k] * fi_93[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_x, pb_y, pb_z, di_49, fh0_15, fh0_16, \
                         fh1_15, fh1_16, fi_94, fi_95, fi_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * fh0_15[k]
                   - f_10 * fh1_15[k]
                   + pb_x[k] * fi_95[k];

        t_222[k] = f_7 * fh0_16[k]
                   - f_8 * fh1_16[k]
                   + pb_x[k] * fi_96[k];

        t_223[k] = pb_z[k] * fi_94[k];

        t_224[k] = f_0 * di_49[k]
                   + pb_y[k] * fi_95[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_x, pb_z, fh0_17, fh0_18, fh0_19, \
                         fh1_17, fh1_18, fh1_19, fi_96, fi_97, fi_98, \
                         fi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_7 * fh0_17[k]
                   - f_8 * fh1_17[k]
                   + pb_x[k] * fi_97[k];

        t_226[k] = f_5 * fh0_18[k]
                   - f_6 * fh1_18[k]
                   + pb_x[k] * fi_98[k];

        t_227[k] = pb_z[k] * fi_96[k];

        t_228[k] = f_5 * fh0_19[k]
                   - f_6 * fh1_19[k]
                   + pb_x[k] * fi_99[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, di_52, fh0_20, fh0_21, \
                         fh1_20, fh1_21, fi_97, fi_98, fi_100, fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * di_52[k]
                   + pb_y[k] * fi_97[k];

        t_230[k] = f_5 * fh0_20[k]
                   - f_6 * fh1_20[k]
                   + pb_x[k] * fi_100[k];

        t_231[k] = f_3 * fh0_21[k]
                   - f_4 * fh1_21[k]
                   + pb_x[k] * fi_101[k];

        t_232[k] = pb_z[k] * fi_98[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_x, pb_y, di_56, fh0_23, fh0_24, fh1_23, \
                         fh1_24, fi_100, fi_102, fi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * fh0_23[k]
                   - f_4 * fh1_23[k]
                   + pb_x[k] * fi_102[k];

        t_234[k] = f_3 * fh0_24[k]
                   - f_4 * fh1_24[k]
                   + pb_x[k] * fi_103[k];

        t_235[k] = f_0 * di_56[k]
                   + pb_y[k] * fi_100[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, pb_x, fh0_25, fh1_25, \
                         fi_104, fi_105, fi_106, fi_107, fi_108, \
                         fi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_3 * fh0_25[k]
                   - f_4 * fh1_25[k]
                   + pb_x[k] * fi_104[k];

        t_237[k] = pb_x[k] * fi_105[k];

        t_238[k] = pb_x[k] * fi_106[k];

        t_239[k] = pb_x[k] * fi_107[k];

        t_240[k] = pb_x[k] * fi_108[k];

        t_241[k] = pb_x[k] * fi_109[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, pb_y, pb_z, di_61, fh0_21, \
                         fh1_21, fi_105, fi_106, fi_110, fi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = pb_x[k] * fi_110[k];

        t_243[k] = pb_x[k] * fi_111[k];

        t_244[k] = f_0 * di_61[k]
                   + f_1 * fh0_21[k]
                   - f_2 * fh1_21[k]
                   + pb_y[k] * fi_105[k];

        t_245[k] = pb_z[k] * fi_105[k];

        t_246[k] = f_3 * fh0_21[k]
                   - f_4 * fh1_21[k]
                   + pb_z[k] * fi_106[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pb_z, fh0_22, fh0_23, fh0_24, fh1_22, fh1_23, \
                         fh1_24, fi_107, fi_108, fi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_5 * fh0_22[k]
                   - f_6 * fh1_22[k]
                   + pb_z[k] * fi_107[k];

        t_248[k] = f_7 * fh0_23[k]
                   - f_8 * fh1_23[k]
                   + pb_z[k] * fi_108[k];

        t_249[k] = f_9 * fh0_24[k]
                   - f_10 * fh1_24[k]
                   + pb_z[k] * fi_109[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, pa_z, pb_y, pb_z, di_46, di_67, \
                         dk_33, dk_34, fh0_25, fh1_25, fi_111, fi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * di_67[k]
                   + pb_y[k] * fi_111[k];

        t_251[k] = f_1 * fh0_25[k]
                   - f_2 * fh1_25[k]
                   + pb_z[k] * fi_111[k];

        t_252[k] = pa_z[k] * dk_33[k];

        t_253[k] = pa_z[k] * dk_34[k];

        t_254[k] = f_11 * di_46[k]
                   + pb_z[k] * fi_112[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pa_z, pb_y, pb_z, di_47, di_48, \
                         di_68, dk_35, dk_36, dk_37, fi_113, fi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pa_z[k] * dk_35[k];

        t_256[k] = f_12 * di_68[k]
                   + pb_y[k] * fi_113[k];

        t_257[k] = f_12 * di_47[k]
                   + pa_z[k] * dk_36[k];

        t_258[k] = pa_z[k] * dk_37[k];

        t_259[k] = f_11 * di_48[k]
                   + pb_z[k] * fi_114[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_z, pb_y, pb_z, di_49, di_50, di_70, \
                         dk_38, dk_39, fi_115, fi_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_12 * di_70[k]
                   + pb_y[k] * fi_115[k];

        t_261[k] = f_0 * di_49[k]
                   + pa_z[k] * dk_38[k];

        t_262[k] = pa_z[k] * dk_39[k];

        t_263[k] = f_11 * di_50[k]
                   + pb_z[k] * fi_116[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_z, pb_y, di_51, di_52, di_72, dk_40, \
                         dk_41, dk_42, fi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_12 * di_51[k]
                   + pa_z[k] * dk_40[k];

        t_265[k] = f_12 * di_72[k]
                   + pb_y[k] * fi_117[k];

        t_266[k] = f_13 * di_52[k]
                   + pa_z[k] * dk_41[k];

        t_267[k] = pa_z[k] * dk_42[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_z, pb_y, pb_z, di_53, di_54, di_55, \
                         di_75, dk_43, dk_44, fi_118, fi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_11 * di_53[k]
                   + pb_z[k] * fi_118[k];

        t_269[k] = f_12 * di_54[k]
                   + pa_z[k] * dk_43[k];

        t_270[k] = f_0 * di_55[k]
                   + pa_z[k] * dk_44[k];

        t_271[k] = f_12 * di_75[k]
                   + pb_y[k] * fi_119[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, t_277, pa_z, pb_x, di_56, dk_45, \
                         fi_120, fi_121, fi_122, fi_123, fi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_14 * di_56[k]
                   + pa_z[k] * dk_45[k];

        t_273[k] = pb_x[k] * fi_120[k];

        t_274[k] = pb_x[k] * fi_121[k];

        t_275[k] = pb_x[k] * fi_122[k];

        t_276[k] = pb_x[k] * fi_123[k];

        t_277[k] = pb_x[k] * fi_124[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, pa_z, pb_x, pb_z, di_61, di_62, \
                         dk_46, dk_47, fi_120, fi_125, fi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pb_x[k] * fi_125[k];

        t_279[k] = pb_x[k] * fi_126[k];

        t_280[k] = pa_z[k] * dk_46[k];

        t_281[k] = f_11 * di_61[k]
                   + pb_z[k] * fi_120[k];

        t_282[k] = f_12 * di_62[k]
                   + pa_z[k] * dk_47[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_z, pb_y, di_63, di_64, di_65, di_84, \
                         dk_48, dk_49, dk_50, fi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_0 * di_63[k]
                   + pa_z[k] * dk_48[k];

        t_284[k] = f_13 * di_64[k]
                   + pa_z[k] * dk_49[k];

        t_285[k] = f_14 * di_65[k]
                   + pa_z[k] * dk_50[k];

        t_286[k] = f_12 * di_84[k]
                   + pb_y[k] * fi_126[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_y, pa_z, pb_y, di_67, di_85, \
                         di_86, dk_52, dk_64, dk_65, dk_66, fi_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_15 * di_67[k]
                   + pa_z[k] * dk_52[k];

        t_288[k] = pa_y[k] * dk_64[k];

        t_289[k] = f_11 * di_85[k]
                   + pb_y[k] * fi_127[k];

        t_290[k] = pa_y[k] * dk_65[k];

        t_291[k] = f_12 * di_86[k]
                   + pa_y[k] * dk_66[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_y, pb_z, di_69, di_87, di_88, \
                         dk_67, dk_68, fi_128, fi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_11 * di_87[k]
                   + pb_y[k] * fi_128[k];

        t_293[k] = pa_y[k] * dk_67[k];

        t_294[k] = f_0 * di_88[k]
                   + pa_y[k] * dk_68[k];

        t_295[k] = f_12 * di_69[k]
                   + pb_z[k] * fi_129[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_y, pb_z, di_71, di_89, di_90, \
                         dk_69, dk_70, fi_130, fi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_11 * di_89[k]
                   + pb_y[k] * fi_130[k];

        t_297[k] = pa_y[k] * dk_69[k];

        t_298[k] = f_13 * di_90[k]
                   + pa_y[k] * dk_70[k];

        t_299[k] = f_12 * di_71[k]
                   + pb_z[k] * fi_131[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_y, pb_y, di_91, di_92, di_93, dk_71, \
                         dk_72, dk_73, fi_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_12 * di_91[k]
                   + pa_y[k] * dk_71[k];

        t_301[k] = f_11 * di_92[k]
                   + pb_y[k] * fi_132[k];

        t_302[k] = pa_y[k] * dk_72[k];

        t_303[k] = f_14 * di_93[k]
                   + pa_y[k] * dk_73[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_y, pb_y, pb_z, di_73, di_94, di_95, \
                         di_96, dk_74, dk_75, fi_133, fi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_12 * di_73[k]
                   + pb_z[k] * fi_133[k];

        t_305[k] = f_0 * di_94[k]
                   + pa_y[k] * dk_74[k];

        t_306[k] = f_12 * di_95[k]
                   + pa_y[k] * dk_75[k];

        t_307[k] = f_11 * di_96[k]
                   + pb_y[k] * fi_134[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, t_314, pa_y, pb_x, dk_76, \
                         fi_135, fi_136, fi_137, fi_138, fi_139, \
                         fi_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_y[k] * dk_76[k];

        t_309[k] = pb_x[k] * fi_135[k];

        t_310[k] = pb_x[k] * fi_136[k];

        t_311[k] = pb_x[k] * fi_137[k];

        t_312[k] = pb_x[k] * fi_138[k];

        t_313[k] = pb_x[k] * fi_139[k];

        t_314[k] = pb_x[k] * fi_140[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_y, pb_x, pb_z, di_78, di_101, di_103, \
                         dk_77, dk_79, fi_135, fi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_x[k] * fi_141[k];

        t_316[k] = f_15 * di_101[k]
                   + pa_y[k] * dk_77[k];

        t_317[k] = f_12 * di_78[k]
                   + pb_z[k] * fi_135[k];

        t_318[k] = f_14 * di_103[k]
                   + pa_y[k] * dk_79[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, pa_y, pb_y, di_104, di_105, \
                         di_106, di_107, dk_80, dk_81, dk_82, dk_83, \
                         fi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_13 * di_104[k]
                   + pa_y[k] * dk_80[k];

        t_320[k] = f_0 * di_105[k]
                   + pa_y[k] * dk_81[k];

        t_321[k] = f_12 * di_106[k]
                   + pa_y[k] * dk_82[k];

        t_322[k] = f_11 * di_107[k]
                   + pb_y[k] * fi_141[k];

        t_323[k] = pa_y[k] * dk_83[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pb_x, pb_y, pb_z, di_85, fh0_26, \
                         fh0_27, fh1_26, fh1_27, fi_142, fi_143, \
                         fi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * fh0_26[k]
                   - f_2 * fh1_26[k]
                   + pb_x[k] * fi_142[k];

        t_325[k] = pb_y[k] * fi_142[k];

        t_326[k] = f_0 * di_85[k]
                   + pb_z[k] * fi_142[k];

        t_327[k] = f_9 * fh0_27[k]
                   - f_10 * fh1_27[k]
                   + pb_x[k] * fi_144[k];

        t_328[k] = pb_y[k] * fi_143[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pb_x, pb_y, pb_z, di_88, fh0_28, fh0_29, \
                         fh1_28, fh1_29, fi_144, fi_145, fi_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_9 * fh0_28[k]
                   - f_10 * fh1_28[k]
                   + pb_x[k] * fi_145[k];

        t_330[k] = f_7 * fh0_29[k]
                   - f_8 * fh1_29[k]
                   + pb_x[k] * fi_146[k];

        t_331[k] = f_0 * di_88[k]
                   + pb_z[k] * fi_144[k];

        t_332[k] = pb_y[k] * fi_145[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pb_x, pb_z, di_90, fh0_30, fh0_31, fh1_30, \
                         fh1_31, fi_146, fi_147, fi_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_7 * fh0_30[k]
                   - f_8 * fh1_30[k]
                   + pb_x[k] * fi_147[k];

        t_334[k] = f_5 * fh0_31[k]
                   - f_6 * fh1_31[k]
                   + pb_x[k] * fi_148[k];

        t_335[k] = f_0 * di_90[k]
                   + pb_z[k] * fi_146[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pb_x, pb_y, fh0_32, fh0_33, fh0_34, \
                         fh1_32, fh1_33, fh1_34, fi_147, fi_149, fi_150, \
                         fi_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_5 * fh0_32[k]
                   - f_6 * fh1_32[k]
                   + pb_x[k] * fi_149[k];

        t_337[k] = pb_y[k] * fi_147[k];

        t_338[k] = f_5 * fh0_33[k]
                   - f_6 * fh1_33[k]
                   + pb_x[k] * fi_150[k];

        t_339[k] = f_3 * fh0_34[k]
                   - f_4 * fh1_34[k]
                   + pb_x[k] * fi_151[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pb_x, pb_y, pb_z, di_93, fh0_35, fh0_36, \
                         fh1_35, fh1_36, fi_148, fi_150, fi_152, \
                         fi_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_0 * di_93[k]
                   + pb_z[k] * fi_148[k];

        t_341[k] = f_3 * fh0_35[k]
                   - f_4 * fh1_35[k]
                   + pb_x[k] * fi_152[k];

        t_342[k] = f_3 * fh0_36[k]
                   - f_4 * fh1_36[k]
                   + pb_x[k] * fi_153[k];

        t_343[k] = pb_y[k] * fi_150[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, t_349, pb_x, fh0_38, fh1_38, \
                         fi_154, fi_155, fi_156, fi_157, fi_158, \
                         fi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_3 * fh0_38[k]
                   - f_4 * fh1_38[k]
                   + pb_x[k] * fi_154[k];

        t_345[k] = pb_x[k] * fi_155[k];

        t_346[k] = pb_x[k] * fi_156[k];

        t_347[k] = pb_x[k] * fi_157[k];

        t_348[k] = pb_x[k] * fi_158[k];

        t_349[k] = pb_x[k] * fi_159[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_x, pb_y, pb_z, di_101, fh0_34, fh1_34, \
                         fi_155, fi_160, fi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_x[k] * fi_160[k];

        t_351[k] = pb_x[k] * fi_161[k];

        t_352[k] = f_1 * fh0_34[k]
                   - f_2 * fh1_34[k]
                   + pb_y[k] * fi_155[k];

        t_353[k] = f_0 * di_101[k]
                   + pb_z[k] * fi_155[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pb_y, fh0_35, fh0_36, fh0_37, fh1_35, fh1_36, \
                         fh1_37, fi_157, fi_158, fi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_9 * fh0_35[k]
                   - f_10 * fh1_35[k]
                   + pb_y[k] * fi_157[k];

        t_355[k] = f_7 * fh0_36[k]
                   - f_8 * fh1_36[k]
                   + pb_y[k] * fi_158[k];

        t_356[k] = f_5 * fh0_37[k]
                   - f_6 * fh1_37[k]
                   + pb_y[k] * fi_159[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pb_y, pb_z, di_107, fh0_38, fh1_38, fi_160, \
                         fi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_3 * fh0_38[k]
                   - f_4 * fh1_38[k]
                   + pb_y[k] * fi_160[k];

        t_358[k] = pb_y[k] * fi_161[k];

        t_359[k] = f_0 * di_107[k]
                   + f_1 * fh0_38[k]
                   - f_2 * fh1_38[k]
                   + pb_z[k] * fi_161[k];
    }
}

auto
compute_prim_fk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 2.0 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_95 = buffer.data(dk + 95);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_13 = buffer.data(fh0 + 13);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_15 = buffer.data(fh0 + 15);
    const auto *fh0_16 = buffer.data(fh0 + 16);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_18 = buffer.data(fh0 + 18);
    const auto *fh0_19 = buffer.data(fh0 + 19);
    const auto *fh0_20 = buffer.data(fh0 + 20);
    const auto *fh0_21 = buffer.data(fh0 + 21);
    const auto *fh0_22 = buffer.data(fh0 + 22);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_36 = buffer.data(fh0 + 36);
    const auto *fh0_37 = buffer.data(fh0 + 37);
    const auto *fh0_38 = buffer.data(fh0 + 38);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_13 = buffer.data(fh1 + 13);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_15 = buffer.data(fh1 + 15);
    const auto *fh1_16 = buffer.data(fh1 + 16);
    const auto *fh1_17 = buffer.data(fh1 + 17);
    const auto *fh1_18 = buffer.data(fh1 + 18);
    const auto *fh1_19 = buffer.data(fh1 + 19);
    const auto *fh1_20 = buffer.data(fh1 + 20);
    const auto *fh1_21 = buffer.data(fh1 + 21);
    const auto *fh1_22 = buffer.data(fh1 + 22);
    const auto *fh1_23 = buffer.data(fh1 + 23);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_36 = buffer.data(fh1 + 36);
    const auto *fh1_37 = buffer.data(fh1 + 37);
    const auto *fh1_38 = buffer.data(fh1 + 38);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, di_0, fh0_0, fh1_0, fi_0, \
                         fi_1, fi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pb_y[k] * fi_0[k];

        t_2[k] = pb_z[k] * fi_0[k];

        t_3[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_4[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, fh0_1, fh0_2, fh0_3, fh1_1, fh1_2, \
                         fh1_3, fi_3, fi_4, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_6[k] = pb_y[k] * fi_4[k];

        t_7[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_8[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, fh0_4, fh0_5, fh1_4, fh1_5, fi_6, \
                         fi_7, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * fh0_4[k]
                 - f_4 * fh1_4[k]
                 + pb_y[k] * fi_6[k];

        t_10[k] = pb_y[k] * fi_7[k];

        t_11[k] = f_7 * fh0_4[k]
                  - f_8 * fh1_4[k]
                  + pb_z[k] * fi_7[k];

        t_12[k] = f_9 * fh0_5[k]
                  - f_10 * fh1_5[k]
                  + pb_y[k] * fi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, fh0_6, fh0_7, fh1_6, fh1_7, fi_9, \
                         fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fh0_6[k]
                  - f_6 * fh1_6[k]
                  + pb_y[k] * fi_9[k];

        t_14[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_15[k] = pb_y[k] * fi_11[k];

        t_16[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_x, pb_y, di_12, di_18, fh0_8, fh0_9, \
                         fh1_8, fh1_9, fi_12, fi_13, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * di_12[k]
                  + pb_x[k] * fi_12[k];

        t_18[k] = f_0 * di_18[k]
                  + pb_x[k] * fi_17[k];

        t_19[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_20[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_9[k]
                  + pb_y[k] * fi_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, fh0_10, fh0_11, fh0_12, fh1_10, fh1_11, \
                         fh1_12, fi_14, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_10[k]
                  + pb_y[k] * fi_14[k];

        t_22[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_11[k]
                  + pb_y[k] * fi_15[k];

        t_23[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_12[k]
                  + pb_y[k] * fi_16[k];

        t_24[k] = pb_y[k] * fi_17[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_y, pb_z, di_0, di_1, dk_0, dk_3, \
                         fh0_12, fh1_12, fi_17, fi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_12[k]
                  + pb_z[k] * fi_17[k];

        t_26[k] = pa_y[k] * dk_0[k];

        t_27[k] = f_11 * di_0[k]
                  + pb_y[k] * fi_18[k];

        t_28[k] = f_12 * di_1[k]
                  + pa_y[k] * dk_3[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_y, di_3, di_5, di_8, dk_4, \
                         dk_5, dk_7, dk_8, dk_11, dk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * dk_4[k];

        t_30[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_5[k];

        t_31[k] = pa_y[k] * dk_7[k];

        t_32[k] = f_13 * di_5[k]
                  + pa_y[k] * dk_8[k];

        t_33[k] = pa_y[k] * dk_11[k];

        t_34[k] = f_14 * di_8[k]
                  + pa_y[k] * dk_12[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pb_x, di_12, di_14, di_15, di_20, \
                         dk_16, dk_17, dk_18, dk_19, fi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_y[k] * dk_16[k];

        t_36[k] = f_12 * di_20[k]
                  + pb_x[k] * fi_19[k];

        t_37[k] = f_15 * di_12[k]
                  + pa_y[k] * dk_17[k];

        t_38[k] = f_14 * di_14[k]
                  + pa_y[k] * dk_18[k];

        t_39[k] = f_13 * di_15[k]
                  + pa_y[k] * dk_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, di_16, di_17, di_18, \
                         dk_0, dk_20, dk_21, dk_23, fi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * di_16[k]
                  + pa_y[k] * dk_20[k];

        t_41[k] = f_12 * di_17[k]
                  + pa_y[k] * dk_21[k];

        t_42[k] = f_11 * di_18[k]
                  + pb_y[k] * fi_20[k];

        t_43[k] = pa_y[k] * dk_23[k];

        t_44[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_z, pb_z, di_0, di_2, di_4, dk_3, \
                         dk_4, dk_5, dk_7, fi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_11 * di_0[k]
                  + pb_z[k] * fi_21[k];

        t_46[k] = pa_z[k] * dk_3[k];

        t_47[k] = f_12 * di_2[k]
                  + pa_z[k] * dk_4[k];

        t_48[k] = pa_z[k] * dk_5[k];

        t_49[k] = f_0 * di_4[k]
                  + pa_z[k] * dk_7[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_x, di_7, di_11, di_22, dk_8, \
                         dk_11, dk_12, dk_16, fi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_z[k] * dk_8[k];

        t_51[k] = f_13 * di_7[k]
                  + pa_z[k] * dk_11[k];

        t_52[k] = pa_z[k] * dk_12[k];

        t_53[k] = f_14 * di_11[k]
                  + pa_z[k] * dk_16[k];

        t_54[k] = f_12 * di_22[k]
                  + pb_x[k] * fi_23[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_z, pb_z, di_12, di_13, di_14, di_15, \
                         dk_17, dk_18, dk_19, dk_20, fi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * dk_17[k];

        t_56[k] = f_11 * di_12[k]
                  + pb_z[k] * fi_22[k];

        t_57[k] = f_12 * di_13[k]
                  + pa_z[k] * dk_18[k];

        t_58[k] = f_0 * di_14[k]
                  + pa_z[k] * dk_19[k];

        t_59[k] = f_13 * di_15[k]
                  + pa_z[k] * dk_20[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pa_z, pb_y, di_16, di_18, di_19, di_23, \
                         dk_21, dk_23, dk_37, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_14 * di_16[k]
                  + pa_z[k] * dk_21[k];

        t_61[k] = f_15 * di_18[k]
                  + pa_z[k] * dk_23[k];

        t_62[k] = f_15 * di_23[k]
                  + pa_x[k] * dk_37[k];

        t_63[k] = f_12 * di_19[k]
                  + pb_y[k] * fi_24[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_x, di_25, di_26, di_27, di_28, \
                         di_29, dk_38, dk_39, dk_40, dk_41, dk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_14 * di_25[k]
                  + pa_x[k] * dk_38[k];

        t_65[k] = f_14 * di_26[k]
                  + pa_x[k] * dk_39[k];

        t_66[k] = f_13 * di_27[k]
                  + pa_x[k] * dk_40[k];

        t_67[k] = f_13 * di_28[k]
                  + pa_x[k] * dk_41[k];

        t_68[k] = f_0 * di_29[k]
                  + pa_x[k] * dk_42[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_x, pb_x, di_31, di_32, di_35, di_36, \
                         dk_44, dk_45, dk_48, dk_54, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_0 * di_31[k]
                  + pa_x[k] * dk_44[k];

        t_70[k] = f_12 * di_32[k]
                  + pa_x[k] * dk_45[k];

        t_71[k] = f_12 * di_35[k]
                  + pa_x[k] * dk_48[k];

        t_72[k] = f_11 * di_36[k]
                  + pb_x[k] * fi_25[k];

        t_73[k] = pa_x[k] * dk_54[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, pa_x, pa_y, dk_31, dk_56, \
                         dk_57, dk_58, dk_59, dk_60, dk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_x[k] * dk_56[k];

        t_75[k] = pa_x[k] * dk_57[k];

        t_76[k] = pa_x[k] * dk_58[k];

        t_77[k] = pa_x[k] * dk_59[k];

        t_78[k] = pa_x[k] * dk_60[k];

        t_79[k] = pa_x[k] * dk_61[k];

        t_80[k] = pa_y[k] * dk_31[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, t_87, pa_y, pa_z, dk_25, dk_26, \
                         dk_27, dk_28, dk_32, dk_33, dk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * dk_25[k];

        t_82[k] = pa_y[k] * dk_32[k];

        t_83[k] = pa_z[k] * dk_26[k];

        t_84[k] = pa_y[k] * dk_33[k];

        t_85[k] = pa_z[k] * dk_27[k];

        t_86[k] = pa_y[k] * dk_34[k];

        t_87[k] = pa_z[k] * dk_28[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, t_94, pa_x, pa_y, dk_35, dk_63, \
                         dk_64, dk_65, dk_66, dk_67, dk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pa_y[k] * dk_35[k];

        t_89[k] = pa_x[k] * dk_63[k];

        t_90[k] = pa_x[k] * dk_64[k];

        t_91[k] = pa_x[k] * dk_65[k];

        t_92[k] = pa_x[k] * dk_66[k];

        t_93[k] = pa_x[k] * dk_67[k];

        t_94[k] = pa_x[k] * dk_68[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_x, pb_z, di_21, di_47, di_49, di_50, \
                         dk_70, dk_72, dk_73, fi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * di_47[k]
                  + pa_x[k] * dk_70[k];

        t_96[k] = f_12 * di_21[k]
                  + pb_z[k] * fi_26[k];

        t_97[k] = f_14 * di_49[k]
                  + pa_x[k] * dk_72[k];

        t_98[k] = f_14 * di_50[k]
                  + pa_x[k] * dk_73[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_x, di_51, di_52, di_53, di_55, \
                         di_56, dk_74, dk_75, dk_76, dk_78, dk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * di_51[k]
                  + pa_x[k] * dk_74[k];

        t_100[k] = f_13 * di_52[k]
                   + pa_x[k] * dk_75[k];

        t_101[k] = f_0 * di_53[k]
                   + pa_x[k] * dk_76[k];

        t_102[k] = f_0 * di_55[k]
                   + pa_x[k] * dk_78[k];

        t_103[k] = f_12 * di_56[k]
                   + pa_x[k] * dk_79[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, t_109, pa_x, pb_x, di_59, di_65, \
                         dk_82, dk_88, dk_89, dk_90, dk_91, fi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_12 * di_59[k]
                   + pa_x[k] * dk_82[k];

        t_105[k] = f_11 * di_65[k]
                   + pb_x[k] * fi_27[k];

        t_106[k] = pa_x[k] * dk_88[k];

        t_107[k] = pa_x[k] * dk_89[k];

        t_108[k] = pa_x[k] * dk_90[k];

        t_109[k] = pa_x[k] * dk_91[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_x, pb_x, pb_y, di_23, dk_92, \
                         dk_93, dk_95, fh0_13, fh1_13, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_x[k] * dk_92[k];

        t_111[k] = pa_x[k] * dk_93[k];

        t_112[k] = pa_x[k] * dk_95[k];

        t_113[k] = f_1 * fh0_13[k]
                   - f_2 * fh1_13[k]
                   + pb_x[k] * fi_28[k];

        t_114[k] = f_0 * di_23[k]
                   + pb_y[k] * fi_28[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_x, fh0_14, fh0_15, fh0_16, fh1_14, fh1_15, \
                         fh1_16, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_9 * fh0_14[k]
                   - f_10 * fh1_14[k]
                   + pb_x[k] * fi_29[k];

        t_116[k] = f_9 * fh0_15[k]
                   - f_10 * fh1_15[k]
                   + pb_x[k] * fi_30[k];

        t_117[k] = f_7 * fh0_16[k]
                   - f_8 * fh1_16[k]
                   + pb_x[k] * fi_31[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, fh0_17, fh0_18, fh0_19, fh1_17, fh1_18, \
                         fh1_19, fi_32, fi_33, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_7 * fh0_17[k]
                   - f_8 * fh1_17[k]
                   + pb_x[k] * fi_32[k];

        t_119[k] = f_5 * fh0_18[k]
                   - f_6 * fh1_18[k]
                   + pb_x[k] * fi_33[k];

        t_120[k] = f_5 * fh0_19[k]
                   - f_6 * fh1_19[k]
                   + pb_x[k] * fi_34[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_x, fh0_20, fh0_21, fh0_23, fh1_20, fh1_21, \
                         fh1_23, fi_35, fi_36, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * fh0_20[k]
                   - f_6 * fh1_20[k]
                   + pb_x[k] * fi_35[k];

        t_122[k] = f_3 * fh0_21[k]
                   - f_4 * fh1_21[k]
                   + pb_x[k] * fi_36[k];

        t_123[k] = f_3 * fh0_23[k]
                   - f_4 * fh1_23[k]
                   + pb_x[k] * fi_37[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pb_x, fh0_24, fh0_25, fh1_24, \
                         fh1_25, fi_38, fi_39, fi_40, fi_42, fi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_3 * fh0_24[k]
                   - f_4 * fh1_24[k]
                   + pb_x[k] * fi_38[k];

        t_125[k] = f_3 * fh0_25[k]
                   - f_4 * fh1_25[k]
                   + pb_x[k] * fi_39[k];

        t_126[k] = pb_x[k] * fi_40[k];

        t_127[k] = pb_x[k] * fi_42[k];

        t_128[k] = pb_x[k] * fi_43[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_y, pb_z, di_36, fh0_21, \
                         fh1_21, fi_40, fi_41, fi_44, fi_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pb_x[k] * fi_44[k];

        t_130[k] = pb_x[k] * fi_45[k];

        t_131[k] = f_0 * di_36[k]
                   + f_1 * fh0_21[k]
                   - f_2 * fh1_21[k]
                   + pb_y[k] * fi_40[k];

        t_132[k] = pb_z[k] * fi_40[k];

        t_133[k] = f_3 * fh0_21[k]
                   - f_4 * fh1_21[k]
                   + pb_z[k] * fi_41[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pb_z, fh0_22, fh0_23, fh0_24, fh1_22, fh1_23, \
                         fh1_24, fi_42, fi_43, fi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_5 * fh0_22[k]
                   - f_6 * fh1_22[k]
                   + pb_z[k] * fi_42[k];

        t_135[k] = f_7 * fh0_23[k]
                   - f_8 * fh1_23[k]
                   + pb_z[k] * fi_43[k];

        t_136[k] = f_9 * fh0_24[k]
                   - f_10 * fh1_24[k]
                   + pb_z[k] * fi_44[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pa_z, pb_y, pb_z, di_23, di_41, \
                         dk_37, dk_38, fh0_25, fh1_25, fi_45, fi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_0 * di_41[k]
                   + pb_y[k] * fi_45[k];

        t_138[k] = f_1 * fh0_25[k]
                   - f_2 * fh1_25[k]
                   + pb_z[k] * fi_45[k];

        t_139[k] = pa_z[k] * dk_37[k];

        t_140[k] = f_11 * di_23[k]
                   + pb_z[k] * fi_46[k];

        t_141[k] = pa_z[k] * dk_38[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, t_147, pa_z, di_24, di_26, di_28, \
                         dk_39, dk_40, dk_41, dk_42, dk_44, dk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_12 * di_24[k]
                   + pa_z[k] * dk_39[k];

        t_143[k] = pa_z[k] * dk_40[k];

        t_144[k] = f_0 * di_26[k]
                   + pa_z[k] * dk_41[k];

        t_145[k] = pa_z[k] * dk_42[k];

        t_146[k] = f_13 * di_28[k]
                   + pa_z[k] * dk_44[k];

        t_147[k] = pa_z[k] * dk_45[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, pa_z, pb_z, di_31, di_36, di_37, \
                         di_38, dk_48, dk_54, dk_56, dk_57, fi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * di_31[k]
                   + pa_z[k] * dk_48[k];

        t_149[k] = pa_z[k] * dk_54[k];

        t_150[k] = f_11 * di_36[k]
                   + pb_z[k] * fi_47[k];

        t_151[k] = f_12 * di_37[k]
                   + pa_z[k] * dk_56[k];

        t_152[k] = f_0 * di_38[k]
                   + pa_z[k] * dk_57[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_z, pb_y, di_39, di_40, di_41, di_46, \
                         dk_58, dk_59, dk_61, fi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_13 * di_39[k]
                   + pa_z[k] * dk_58[k];

        t_154[k] = f_14 * di_40[k]
                   + pa_z[k] * dk_59[k];

        t_155[k] = f_12 * di_46[k]
                   + pb_y[k] * fi_48[k];

        t_156[k] = f_15 * di_41[k]
                   + pa_z[k] * dk_61[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, pa_y, di_48, di_49, dk_70, \
                         dk_71, dk_72, dk_73, dk_74, dk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pa_y[k] * dk_70[k];

        t_158[k] = pa_y[k] * dk_71[k];

        t_159[k] = f_12 * di_48[k]
                   + pa_y[k] * dk_72[k];

        t_160[k] = pa_y[k] * dk_73[k];

        t_161[k] = f_0 * di_49[k]
                   + pa_y[k] * dk_74[k];

        t_162[k] = pa_y[k] * dk_75[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_y, di_51, di_53, di_60, dk_76, \
                         dk_78, dk_79, dk_82, dk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_13 * di_51[k]
                   + pa_y[k] * dk_76[k];

        t_164[k] = pa_y[k] * dk_78[k];

        t_165[k] = f_14 * di_53[k]
                   + pa_y[k] * dk_79[k];

        t_166[k] = pa_y[k] * dk_82[k];

        t_167[k] = f_15 * di_60[k]
                   + pa_y[k] * dk_88[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pb_z, di_42, di_61, di_62, di_63, \
                         dk_90, dk_91, dk_92, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_12 * di_42[k]
                   + pb_z[k] * fi_49[k];

        t_169[k] = f_14 * di_61[k]
                   + pa_y[k] * dk_90[k];

        t_170[k] = f_13 * di_62[k]
                   + pa_y[k] * dk_91[k];

        t_171[k] = f_0 * di_63[k]
                   + pa_y[k] * dk_92[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_x, pb_y, di_64, di_65, dk_93, \
                         dk_95, fh0_26, fh1_26, fi_50, fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_12 * di_64[k]
                   + pa_y[k] * dk_93[k];

        t_173[k] = f_11 * di_65[k]
                   + pb_y[k] * fi_50[k];

        t_174[k] = pa_y[k] * dk_95[k];

        t_175[k] = f_1 * fh0_26[k]
                   - f_2 * fh1_26[k]
                   + pb_x[k] * fi_51[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pb_x, pb_z, di_47, fh0_27, fh0_28, fh1_27, \
                         fh1_28, fi_51, fi_52, fi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_0 * di_47[k]
                   + pb_z[k] * fi_51[k];

        t_177[k] = f_9 * fh0_27[k]
                   - f_10 * fh1_27[k]
                   + pb_x[k] * fi_52[k];

        t_178[k] = f_9 * fh0_28[k]
                   - f_10 * fh1_28[k]
                   + pb_x[k] * fi_53[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pb_x, fh0_29, fh0_30, fh0_31, fh1_29, fh1_30, \
                         fh1_31, fi_54, fi_55, fi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_7 * fh0_29[k]
                   - f_8 * fh1_29[k]
                   + pb_x[k] * fi_54[k];

        t_180[k] = f_7 * fh0_30[k]
                   - f_8 * fh1_30[k]
                   + pb_x[k] * fi_55[k];

        t_181[k] = f_5 * fh0_31[k]
                   - f_6 * fh1_31[k]
                   + pb_x[k] * fi_56[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pb_x, fh0_32, fh0_33, fh0_34, fh1_32, fh1_33, \
                         fh1_34, fi_57, fi_58, fi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_5 * fh0_32[k]
                   - f_6 * fh1_32[k]
                   + pb_x[k] * fi_57[k];

        t_183[k] = f_5 * fh0_33[k]
                   - f_6 * fh1_33[k]
                   + pb_x[k] * fi_58[k];

        t_184[k] = f_3 * fh0_34[k]
                   - f_4 * fh1_34[k]
                   + pb_x[k] * fi_59[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, fh0_35, fh0_36, fh0_38, fh1_35, \
                         fh1_36, fh1_38, fi_60, fi_61, fi_62, fi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_3 * fh0_35[k]
                   - f_4 * fh1_35[k]
                   + pb_x[k] * fi_60[k];

        t_186[k] = f_3 * fh0_36[k]
                   - f_4 * fh1_36[k]
                   + pb_x[k] * fi_61[k];

        t_187[k] = f_3 * fh0_38[k]
                   - f_4 * fh1_38[k]
                   + pb_x[k] * fi_62[k];

        t_188[k] = pb_x[k] * fi_63[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pb_x, pb_y, fh0_34, fh1_34, fi_63, \
                         fi_64, fi_65, fi_66, fi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_x[k] * fi_64[k];

        t_190[k] = pb_x[k] * fi_65[k];

        t_191[k] = pb_x[k] * fi_66[k];

        t_192[k] = pb_x[k] * fi_68[k];

        t_193[k] = f_1 * fh0_34[k]
                   - f_2 * fh1_34[k]
                   + pb_y[k] * fi_63[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pb_y, pb_z, di_60, fh0_35, fh0_36, fh1_35, \
                         fh1_36, fi_63, fi_64, fi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_0 * di_60[k]
                   + pb_z[k] * fi_63[k];

        t_195[k] = f_9 * fh0_35[k]
                   - f_10 * fh1_35[k]
                   + pb_y[k] * fi_64[k];

        t_196[k] = f_7 * fh0_36[k]
                   - f_8 * fh1_36[k]
                   + pb_y[k] * fi_65[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pb_z, di_65, fh0_37, fh0_38, \
                         fh1_37, fh1_38, fi_66, fi_67, fi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_5 * fh0_37[k]
                   - f_6 * fh1_37[k]
                   + pb_y[k] * fi_66[k];

        t_198[k] = f_3 * fh0_38[k]
                   - f_4 * fh1_38[k]
                   + pb_y[k] * fi_67[k];

        t_199[k] = pb_y[k] * fi_68[k];

        t_200[k] = f_0 * di_65[k]
                   + f_1 * fh0_38[k]
                   - f_2 * fh1_38[k]
                   + pb_z[k] * fi_68[k];
    }
}

auto
compute_prim_fk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 2.0 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_18 = buffer.data(fh0 + 18);
    const auto *fh0_19 = buffer.data(fh0 + 19);
    const auto *fh0_20 = buffer.data(fh0 + 20);
    const auto *fh0_21 = buffer.data(fh0 + 21);
    const auto *fh0_22 = buffer.data(fh0 + 22);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_36 = buffer.data(fh0 + 36);
    const auto *fh0_37 = buffer.data(fh0 + 37);
    const auto *fh0_38 = buffer.data(fh0 + 38);
    const auto *fh0_39 = buffer.data(fh0 + 39);
    const auto *fh0_40 = buffer.data(fh0 + 40);
    const auto *fh0_41 = buffer.data(fh0 + 41);
    const auto *fh0_42 = buffer.data(fh0 + 42);
    const auto *fh0_43 = buffer.data(fh0 + 43);
    const auto *fh0_44 = buffer.data(fh0 + 44);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_13 = buffer.data(fh1 + 13);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_37 = buffer.data(fh1 + 37);
    const auto *fh1_38 = buffer.data(fh1 + 38);
    const auto *fh1_39 = buffer.data(fh1 + 39);
    const auto *fh1_40 = buffer.data(fh1 + 40);
    const auto *fh1_41 = buffer.data(fh1 + 41);
    const auto *fh1_42 = buffer.data(fh1 + 42);
    const auto *fh1_43 = buffer.data(fh1 + 43);
    const auto *fh1_44 = buffer.data(fh1 + 44);
    const auto *fh1_45 = buffer.data(fh1 + 45);
    const auto *fh1_46 = buffer.data(fh1 + 46);
    const auto *fh1_47 = buffer.data(fh1 + 47);
    const auto *fh1_48 = buffer.data(fh1 + 48);
    const auto *fh1_65 = buffer.data(fh1 + 65);
    const auto *fh1_67 = buffer.data(fh1 + 67);
    const auto *fh1_68 = buffer.data(fh1 + 68);
    const auto *fh1_69 = buffer.data(fh1 + 69);
    const auto *fh1_70 = buffer.data(fh1 + 70);
    const auto *fh1_71 = buffer.data(fh1 + 71);
    const auto *fh1_72 = buffer.data(fh1 + 72);
    const auto *fh1_73 = buffer.data(fh1 + 73);
    const auto *fh1_74 = buffer.data(fh1 + 74);
    const auto *fh1_75 = buffer.data(fh1 + 75);
    const auto *fh1_76 = buffer.data(fh1 + 76);
    const auto *fh1_77 = buffer.data(fh1 + 77);
    const auto *fh1_78 = buffer.data(fh1 + 78);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, di_0, fh0_0, fh0_1, fh1_0, \
                         fh1_1, fi_0, fi_1, fi_2, fi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_2[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];

        t_3[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, fh0_2, fh0_3, fh0_4, fh1_2, fh1_3, \
                         fh1_4, fi_4, fi_5, fi_6, fi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_5[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = f_3 * fh0_4[k]
                 - f_4 * fh1_4[k]
                 + pb_y[k] * fi_6[k];

        t_7[k] = f_7 * fh0_4[k]
                 - f_8 * fh1_4[k]
                 + pb_z[k] * fi_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, fh0_5, fh0_6, fh0_7, fh1_5, fh1_6, \
                         fh1_7, fi_8, fi_9, fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * fh0_5[k]
                 - f_10 * fh1_5[k]
                 + pb_y[k] * fi_8[k];

        t_9[k] = f_5 * fh0_6[k]
                 - f_6 * fh1_6[k]
                 + pb_y[k] * fi_9[k];

        t_10[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_11[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, di_12, di_17, fh0_8, fh0_9, \
                         fh1_8, fh1_10, fi_12, fi_13, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * di_12[k]
                  + pb_x[k] * fi_12[k];

        t_13[k] = f_0 * di_17[k]
                  + pb_x[k] * fi_17[k];

        t_14[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_15[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_10[k]
                  + pb_y[k] * fi_13[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, fh0_10, fh0_11, fh0_12, fh1_11, \
                         fh1_12, fh1_13, fi_14, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_11[k]
                  + pb_y[k] * fi_14[k];

        t_17[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_12[k]
                  + pb_y[k] * fi_15[k];

        t_18[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_13[k]
                  + pb_y[k] * fi_16[k];

        t_19[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_13[k]
                  + pb_z[k] * fi_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, di_0, di_1, di_3, di_5, \
                         dk_0, dk_1, dk_3, dk_5, fi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * dk_0[k];

        t_21[k] = f_11 * di_0[k]
                  + pb_y[k] * fi_18[k];

        t_22[k] = f_12 * di_1[k]
                  + pa_y[k] * dk_1[k];

        t_23[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_3[k];

        t_24[k] = f_13 * di_5[k]
                  + pa_y[k] * dk_5[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_x, di_8, di_12, di_22, dk_0, \
                         dk_8, dk_12, fi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_14 * di_8[k]
                  + pa_y[k] * dk_8[k];

        t_26[k] = f_12 * di_22[k]
                  + pb_x[k] * fi_22[k];

        t_27[k] = f_15 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_28[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_z, pb_z, di_0, di_2, di_4, di_6, dk_2, \
                         dk_4, dk_6, fi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * di_0[k]
                  + pb_z[k] * fi_23[k];

        t_30[k] = f_12 * di_2[k]
                  + pa_z[k] * dk_2[k];

        t_31[k] = f_0 * di_4[k]
                  + pa_z[k] * dk_4[k];

        t_32[k] = f_12 * di_6[k]
                  + pa_z[k] * dk_6[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, di_7, di_9, di_10, di_11, dk_7, dk_9, \
                         dk_10, dk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_13 * di_7[k]
                  + pa_z[k] * dk_7[k];

        t_34[k] = f_12 * di_9[k]
                  + pa_z[k] * dk_9[k];

        t_35[k] = f_0 * di_10[k]
                  + pa_z[k] * dk_10[k];

        t_36[k] = f_14 * di_11[k]
                  + pa_z[k] * dk_11[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, di_13, di_14, di_15, di_28, \
                         dk_13, dk_14, dk_15, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_12 * di_28[k]
                  + pb_x[k] * fi_31[k];

        t_38[k] = f_12 * di_13[k]
                  + pa_z[k] * dk_13[k];

        t_39[k] = f_0 * di_14[k]
                  + pa_z[k] * dk_14[k];

        t_40[k] = f_13 * di_15[k]
                  + pa_z[k] * dk_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pa_z, pb_y, di_16, di_17, di_18, di_29, \
                         dk_16, dk_17, dk_18, fi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_14 * di_16[k]
                  + pa_z[k] * dk_16[k];

        t_42[k] = f_15 * di_17[k]
                  + pa_z[k] * dk_17[k];

        t_43[k] = f_15 * di_29[k]
                  + pa_x[k] * dk_18[k];

        t_44[k] = f_12 * di_18[k]
                  + pb_y[k] * fi_32[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, di_31, di_33, di_36, di_40, dk_19, \
                         dk_21, dk_23, dk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_14 * di_31[k]
                  + pa_x[k] * dk_19[k];

        t_46[k] = f_13 * di_33[k]
                  + pa_x[k] * dk_21[k];

        t_47[k] = f_0 * di_36[k]
                  + pa_x[k] * dk_23[k];

        t_48[k] = f_12 * di_40[k]
                  + pa_x[k] * dk_26[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pb_x, pb_z, di_23, di_41, di_55, dk_30, \
                         dk_36, fi_36, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * di_41[k]
                  + pb_x[k] * fi_36[k];

        t_50[k] = pa_x[k] * dk_30[k];

        t_51[k] = f_15 * di_55[k]
                  + pa_x[k] * dk_36[k];

        t_52[k] = f_12 * di_23[k]
                  + pb_z[k] * fi_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, di_59, di_62, di_66, di_67, dk_38, \
                         dk_40, dk_43, dk_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_14 * di_59[k]
                  + pa_x[k] * dk_38[k];

        t_54[k] = f_13 * di_62[k]
                  + pa_x[k] * dk_40[k];

        t_55[k] = f_0 * di_66[k]
                  + pa_x[k] * dk_43[k];

        t_56[k] = f_12 * di_67[k]
                  + pa_x[k] * dk_47[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_x, pb_y, di_29, di_74, dk_53, \
                         fh0_17, fh1_35, fi_42, fi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_11 * di_74[k]
                  + pb_x[k] * fi_42[k];

        t_58[k] = pa_x[k] * dk_53[k];

        t_59[k] = f_1 * fh0_17[k]
                  - f_2 * fh1_35[k]
                  + pb_x[k] * fi_43[k];

        t_60[k] = f_0 * di_29[k]
                  + pb_y[k] * fi_43[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, fh0_18, fh0_19, fh0_20, fh1_37, fh1_38, \
                         fh1_39, fi_44, fi_45, fi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * fh0_18[k]
                  - f_10 * fh1_37[k]
                  + pb_x[k] * fi_44[k];

        t_62[k] = f_9 * fh0_19[k]
                  - f_10 * fh1_38[k]
                  + pb_x[k] * fi_45[k];

        t_63[k] = f_7 * fh0_20[k]
                  - f_8 * fh1_39[k]
                  + pb_x[k] * fi_46[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, fh0_21, fh0_22, fh0_23, fh1_40, fh1_41, \
                         fh1_42, fi_47, fi_48, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_7 * fh0_21[k]
                  - f_8 * fh1_40[k]
                  + pb_x[k] * fi_47[k];

        t_65[k] = f_5 * fh0_22[k]
                  - f_6 * fh1_41[k]
                  + pb_x[k] * fi_48[k];

        t_66[k] = f_5 * fh0_23[k]
                  - f_6 * fh1_42[k]
                  + pb_x[k] * fi_49[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, fh0_24, fh0_25, fh0_27, fh1_43, fh1_44, \
                         fh1_46, fi_50, fi_51, fi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * fh0_24[k]
                  - f_6 * fh1_43[k]
                  + pb_x[k] * fi_50[k];

        t_68[k] = f_3 * fh0_25[k]
                  - f_4 * fh1_44[k]
                  + pb_x[k] * fi_51[k];

        t_69[k] = f_3 * fh0_27[k]
                  - f_4 * fh1_46[k]
                  + pb_x[k] * fi_52[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_x, pb_y, di_41, fh0_25, fh0_28, fh0_29, fh1_44, \
                         fh1_47, fh1_48, fi_53, fi_54, fi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * fh0_28[k]
                  - f_4 * fh1_47[k]
                  + pb_x[k] * fi_53[k];

        t_71[k] = f_3 * fh0_29[k]
                  - f_4 * fh1_48[k]
                  + pb_x[k] * fi_54[k];

        t_72[k] = f_0 * di_41[k]
                  + f_1 * fh0_25[k]
                  - f_2 * fh1_44[k]
                  + pb_y[k] * fi_55[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_z, fh0_25, fh0_26, fh0_27, fh1_44, fh1_45, \
                         fh1_46, fi_56, fi_57, fi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * fh0_25[k]
                  - f_4 * fh1_44[k]
                  + pb_z[k] * fi_56[k];

        t_74[k] = f_5 * fh0_26[k]
                  - f_6 * fh1_45[k]
                  + pb_z[k] * fi_57[k];

        t_75[k] = f_7 * fh0_27[k]
                  - f_8 * fh1_46[k]
                  + pb_z[k] * fi_58[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, di_30, di_47, dk_20, \
                         fh0_28, fh0_29, fh1_47, fh1_48, fi_59, fi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_9 * fh0_28[k]
                  - f_10 * fh1_47[k]
                  + pb_z[k] * fi_59[k];

        t_77[k] = f_0 * di_47[k]
                  + pb_y[k] * fi_61[k];

        t_78[k] = f_1 * fh0_29[k]
                  - f_2 * fh1_48[k]
                  + pb_z[k] * fi_61[k];

        t_79[k] = f_12 * di_30[k]
                  + pa_z[k] * dk_20[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_z, di_32, di_34, di_35, di_37, \
                         di_38, dk_22, dk_24, dk_25, dk_27, dk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * di_32[k]
                  + pa_z[k] * dk_22[k];

        t_81[k] = f_12 * di_34[k]
                  + pa_z[k] * dk_24[k];

        t_82[k] = f_13 * di_35[k]
                  + pa_z[k] * dk_25[k];

        t_83[k] = f_12 * di_37[k]
                  + pa_z[k] * dk_27[k];

        t_84[k] = f_0 * di_38[k]
                  + pa_z[k] * dk_28[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_z, di_39, di_41, di_42, di_43, \
                         dk_29, dk_30, dk_31, dk_32, fi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_14 * di_39[k]
                  + pa_z[k] * dk_29[k];

        t_86[k] = pa_z[k] * dk_30[k];

        t_87[k] = f_11 * di_41[k]
                  + pb_z[k] * fi_65[k];

        t_88[k] = f_12 * di_42[k]
                  + pa_z[k] * dk_31[k];

        t_89[k] = f_0 * di_43[k]
                  + pa_z[k] * dk_32[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_z, pb_y, di_44, di_45, di_47, di_54, \
                         dk_33, dk_34, dk_35, fi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_13 * di_44[k]
                  + pa_z[k] * dk_33[k];

        t_91[k] = f_14 * di_45[k]
                  + pa_z[k] * dk_34[k];

        t_92[k] = f_12 * di_54[k]
                  + pb_y[k] * fi_71[k];

        t_93[k] = f_15 * di_47[k]
                  + pa_z[k] * dk_35[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pa_y, di_56, di_58, di_60, di_61, \
                         di_63, dk_37, dk_39, dk_41, dk_42, dk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_12 * di_56[k]
                  + pa_y[k] * dk_37[k];

        t_95[k] = f_0 * di_58[k]
                  + pa_y[k] * dk_39[k];

        t_96[k] = f_13 * di_60[k]
                  + pa_y[k] * dk_41[k];

        t_97[k] = f_12 * di_61[k]
                  + pa_y[k] * dk_42[k];

        t_98[k] = f_14 * di_63[k]
                  + pa_y[k] * dk_44[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_y, pb_z, di_48, di_64, di_65, di_68, \
                         dk_45, dk_46, dk_48, fi_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_0 * di_64[k]
                  + pa_y[k] * dk_45[k];

        t_100[k] = f_12 * di_65[k]
                   + pa_y[k] * dk_46[k];

        t_101[k] = f_15 * di_68[k]
                   + pa_y[k] * dk_48[k];

        t_102[k] = f_12 * di_48[k]
                   + pb_z[k] * fi_75[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_y, di_70, di_71, di_72, di_73, dk_49, \
                         dk_50, dk_51, dk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_14 * di_70[k]
                   + pa_y[k] * dk_49[k];

        t_104[k] = f_13 * di_71[k]
                   + pa_y[k] * dk_50[k];

        t_105[k] = f_0 * di_72[k]
                   + pa_y[k] * dk_51[k];

        t_106[k] = f_12 * di_73[k]
                   + pa_y[k] * dk_52[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pb_x, pb_y, pb_z, di_55, di_74, \
                         dk_53, fh0_32, fh1_65, fi_81, fi_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_11 * di_74[k]
                   + pb_y[k] * fi_81[k];

        t_108[k] = pa_y[k] * dk_53[k];

        t_109[k] = f_1 * fh0_32[k]
                   - f_2 * fh1_65[k]
                   + pb_x[k] * fi_82[k];

        t_110[k] = f_0 * di_55[k]
                   + pb_z[k] * fi_82[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, fh0_33, fh0_34, fh0_35, fh1_67, fh1_68, \
                         fh1_69, fi_84, fi_85, fi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_9 * fh0_33[k]
                   - f_10 * fh1_67[k]
                   + pb_x[k] * fi_84[k];

        t_112[k] = f_9 * fh0_34[k]
                   - f_10 * fh1_68[k]
                   + pb_x[k] * fi_85[k];

        t_113[k] = f_7 * fh0_35[k]
                   - f_8 * fh1_69[k]
                   + pb_x[k] * fi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, fh0_36, fh0_37, fh0_38, fh1_70, fh1_71, \
                         fh1_72, fi_87, fi_88, fi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * fh0_36[k]
                   - f_8 * fh1_70[k]
                   + pb_x[k] * fi_87[k];

        t_115[k] = f_5 * fh0_37[k]
                   - f_6 * fh1_71[k]
                   + pb_x[k] * fi_88[k];

        t_116[k] = f_5 * fh0_38[k]
                   - f_6 * fh1_72[k]
                   + pb_x[k] * fi_89[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_x, fh0_39, fh0_40, fh0_41, fh1_73, fh1_74, \
                         fh1_75, fi_90, fi_91, fi_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * fh0_39[k]
                   - f_6 * fh1_73[k]
                   + pb_x[k] * fi_90[k];

        t_118[k] = f_3 * fh0_40[k]
                   - f_4 * fh1_74[k]
                   + pb_x[k] * fi_91[k];

        t_119[k] = f_3 * fh0_41[k]
                   - f_4 * fh1_75[k]
                   + pb_x[k] * fi_92[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pb_x, pb_y, fh0_40, fh0_42, fh0_44, fh1_74, \
                         fh1_76, fh1_78, fi_93, fi_94, fi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_3 * fh0_42[k]
                   - f_4 * fh1_76[k]
                   + pb_x[k] * fi_93[k];

        t_121[k] = f_3 * fh0_44[k]
                   - f_4 * fh1_78[k]
                   + pb_x[k] * fi_94[k];

        t_122[k] = f_1 * fh0_40[k]
                   - f_2 * fh1_74[k]
                   + pb_y[k] * fi_95[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_y, pb_z, di_68, fh0_41, fh0_42, fh1_75, \
                         fh1_76, fi_95, fi_97, fi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * di_68[k]
                   + pb_z[k] * fi_95[k];

        t_124[k] = f_9 * fh0_41[k]
                   - f_10 * fh1_75[k]
                   + pb_y[k] * fi_97[k];

        t_125[k] = f_7 * fh0_42[k]
                   - f_8 * fh1_76[k]
                   + pb_y[k] * fi_98[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pb_y, pb_z, di_74, fh0_43, fh0_44, fh1_77, \
                         fh1_78, fi_99, fi_100, fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_5 * fh0_43[k]
                   - f_6 * fh1_77[k]
                   + pb_y[k] * fi_99[k];

        t_127[k] = f_3 * fh0_44[k]
                   - f_4 * fh1_78[k]
                   + pb_y[k] * fi_100[k];

        t_128[k] = f_0 * di_74[k]
                   + f_1 * fh0_44[k]
                   - f_2 * fh1_78[k]
                   + pb_z[k] * fi_101[k];
    }
}

auto
compute_prim_fk_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 1.0 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_95 = buffer.data(dk + 95);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_21 = buffer.data(fh0 + 21);
    const auto *fh0_22 = buffer.data(fh0 + 22);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_38 = buffer.data(fh0 + 38);
    const auto *fh0_39 = buffer.data(fh0 + 39);
    const auto *fh0_40 = buffer.data(fh0 + 40);
    const auto *fh0_41 = buffer.data(fh0 + 41);
    const auto *fh0_42 = buffer.data(fh0 + 42);
    const auto *fh0_43 = buffer.data(fh0 + 43);
    const auto *fh0_44 = buffer.data(fh0 + 44);
    const auto *fh0_45 = buffer.data(fh0 + 45);
    const auto *fh0_46 = buffer.data(fh0 + 46);
    const auto *fh0_47 = buffer.data(fh0 + 47);
    const auto *fh0_48 = buffer.data(fh0 + 48);
    const auto *fh0_49 = buffer.data(fh0 + 49);
    const auto *fh0_50 = buffer.data(fh0 + 50);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_23 = buffer.data(fh1 + 23);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_41 = buffer.data(fh1 + 41);
    const auto *fh1_42 = buffer.data(fh1 + 42);
    const auto *fh1_43 = buffer.data(fh1 + 43);
    const auto *fh1_44 = buffer.data(fh1 + 44);
    const auto *fh1_45 = buffer.data(fh1 + 45);
    const auto *fh1_46 = buffer.data(fh1 + 46);
    const auto *fh1_47 = buffer.data(fh1 + 47);
    const auto *fh1_48 = buffer.data(fh1 + 48);
    const auto *fh1_49 = buffer.data(fh1 + 49);
    const auto *fh1_50 = buffer.data(fh1 + 50);
    const auto *fh1_51 = buffer.data(fh1 + 51);
    const auto *fh1_52 = buffer.data(fh1 + 52);
    const auto *fh1_53 = buffer.data(fh1 + 53);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, di_0, fh0_0, fh1_0, fi_0, \
                         fi_1, fi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pb_y[k] * fi_0[k];

        t_2[k] = pb_z[k] * fi_0[k];

        t_3[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_4[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fh0_1, fh0_2, fh0_3, fh1_1, \
                         fh1_2, fh1_3, fi_3, fi_4, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_6[k] = pb_z[k] * fi_3[k];

        t_7[k] = pb_y[k] * fi_4[k];

        t_8[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_9[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pb_y, pb_z, fh0_4, fh0_5, fh1_4, \
                         fh1_5, fi_5, fi_6, fi_7, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * fi_5[k];

        t_11[k] = f_3 * fh0_4[k]
                  - f_4 * fh1_4[k]
                  + pb_y[k] * fi_6[k];

        t_12[k] = pb_y[k] * fi_7[k];

        t_13[k] = f_7 * fh0_4[k]
                  - f_8 * fh1_4[k]
                  + pb_z[k] * fi_7[k];

        t_14[k] = f_9 * fh0_5[k]
                  - f_10 * fh1_5[k]
                  + pb_y[k] * fi_8[k];

        t_15[k] = pb_z[k] * fi_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, fh0_6, fh0_7, fh1_6, fh1_7, fi_9, \
                         fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * fh0_6[k]
                  - f_6 * fh1_6[k]
                  + pb_y[k] * fi_9[k];

        t_17[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_18[k] = pb_y[k] * fi_11[k];

        t_19[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, fh0_8, fh0_9, fh0_10, fh1_8, \
                         fh1_9, fh1_10, fi_12, fi_14, fi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_21[k] = pb_z[k] * fi_12[k];

        t_22[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_9[k]
                  + pb_y[k] * fi_14[k];

        t_23[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_10[k]
                  + pb_y[k] * fi_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_y, pb_z, dk_0, fh0_11, fh0_12, \
                         fh1_11, fh1_12, fi_16, fi_17, fi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_11[k]
                  + pb_y[k] * fi_16[k];

        t_25[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_12[k]
                  + pb_y[k] * fi_17[k];

        t_26[k] = pb_y[k] * fi_18[k];

        t_27[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_12[k]
                  + pb_z[k] * fi_18[k];

        t_28[k] = pa_y[k] * dk_0[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_y, di_1, di_3, di_5, dk_3, \
                         dk_4, dk_5, dk_7, dk_8, dk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * di_1[k]
                  + pa_y[k] * dk_3[k];

        t_30[k] = pa_y[k] * dk_4[k];

        t_31[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_5[k];

        t_32[k] = pa_y[k] * dk_7[k];

        t_33[k] = f_12 * di_5[k]
                  + pa_y[k] * dk_8[k];

        t_34[k] = pa_y[k] * dk_11[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, di_9, di_14, di_16, di_17, dk_12, \
                         dk_16, dk_17, dk_19, dk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * di_9[k]
                  + pa_y[k] * dk_12[k];

        t_36[k] = pa_y[k] * dk_16[k];

        t_37[k] = f_14 * di_14[k]
                  + pa_y[k] * dk_17[k];

        t_38[k] = f_13 * di_16[k]
                  + pa_y[k] * dk_19[k];

        t_39[k] = f_12 * di_17[k]
                  + pa_y[k] * dk_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, di_18, di_19, di_20, \
                         dk_0, dk_21, dk_22, dk_23, fi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * di_18[k]
                  + pa_y[k] * dk_21[k];

        t_41[k] = f_11 * di_19[k]
                  + pa_y[k] * dk_22[k];

        t_42[k] = f_15 * di_20[k]
                  + pb_y[k] * fi_21[k];

        t_43[k] = pa_y[k] * dk_23[k];

        t_44[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_z, pb_y, pb_z, di_0, di_2, dk_3, \
                         dk_4, dk_5, fi_22, fi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_15 * di_0[k]
                  + pb_z[k] * fi_22[k];

        t_46[k] = pa_z[k] * dk_3[k];

        t_47[k] = f_11 * di_2[k]
                  + pa_z[k] * dk_4[k];

        t_48[k] = pa_z[k] * dk_5[k];

        t_49[k] = pb_y[k] * fi_23[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_y, di_4, di_6, di_8, dk_7, \
                         dk_8, dk_10, dk_11, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * di_4[k]
                  + pa_z[k] * dk_7[k];

        t_51[k] = pa_z[k] * dk_8[k];

        t_52[k] = f_11 * di_6[k]
                  + pa_z[k] * dk_10[k];

        t_53[k] = pb_y[k] * fi_24[k];

        t_54[k] = f_12 * di_8[k]
                  + pa_z[k] * dk_11[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_z, pb_y, di_10, di_11, di_13, dk_12, \
                         dk_14, dk_15, dk_16, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * dk_12[k];

        t_56[k] = f_11 * di_10[k]
                  + pa_z[k] * dk_14[k];

        t_57[k] = f_0 * di_11[k]
                  + pa_z[k] * dk_15[k];

        t_58[k] = pb_y[k] * fi_25[k];

        t_59[k] = f_13 * di_13[k]
                  + pa_z[k] * dk_16[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_z, pb_z, di_14, di_15, di_16, di_17, \
                         dk_17, dk_19, dk_20, dk_21, fi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_z[k] * dk_17[k];

        t_61[k] = f_15 * di_14[k]
                  + pb_z[k] * fi_26[k];

        t_62[k] = f_11 * di_15[k]
                  + pa_z[k] * dk_19[k];

        t_63[k] = f_0 * di_16[k]
                  + pa_z[k] * dk_20[k];

        t_64[k] = f_12 * di_17[k]
                  + pa_z[k] * dk_21[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pa_z, pb_y, di_18, di_20, di_22, dk_22, \
                         dk_23, dk_37, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_13 * di_18[k]
                  + pa_z[k] * dk_22[k];

        t_66[k] = pb_y[k] * fi_31[k];

        t_67[k] = f_14 * di_20[k]
                  + pa_z[k] * dk_23[k];

        t_68[k] = f_14 * di_22[k]
                  + pa_x[k] * dk_37[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_x, di_24, di_25, di_26, di_28, \
                         di_29, dk_39, dk_40, dk_41, dk_43, dk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_13 * di_24[k]
                  + pa_x[k] * dk_39[k];

        t_70[k] = f_13 * di_25[k]
                  + pa_x[k] * dk_40[k];

        t_71[k] = f_12 * di_26[k]
                  + pa_x[k] * dk_41[k];

        t_72[k] = f_12 * di_28[k]
                  + pa_x[k] * dk_43[k];

        t_73[k] = f_0 * di_29[k]
                  + pa_x[k] * dk_44[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pa_x, pb_x, di_32, di_33, di_36, di_37, \
                         dk_47, dk_48, dk_52, dk_54, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * di_32[k]
                  + pa_x[k] * dk_47[k];

        t_75[k] = f_11 * di_33[k]
                  + pa_x[k] * dk_48[k];

        t_76[k] = f_11 * di_36[k]
                  + pa_x[k] * dk_52[k];

        t_77[k] = f_15 * di_37[k]
                  + pb_x[k] * fi_37[k];

        t_78[k] = pa_x[k] * dk_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, t_84, t_85, pa_x, pa_y, dk_31, dk_56, \
                         dk_57, dk_58, dk_59, dk_60, dk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pa_x[k] * dk_56[k];

        t_80[k] = pa_x[k] * dk_57[k];

        t_81[k] = pa_x[k] * dk_58[k];

        t_82[k] = pa_x[k] * dk_59[k];

        t_83[k] = pa_x[k] * dk_60[k];

        t_84[k] = pa_x[k] * dk_61[k];

        t_85[k] = pa_y[k] * dk_31[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, dk_25, dk_26, \
                         dk_27, dk_28, dk_32, dk_33, dk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_z[k] * dk_25[k];

        t_87[k] = pa_y[k] * dk_32[k];

        t_88[k] = pa_z[k] * dk_26[k];

        t_89[k] = pa_y[k] * dk_33[k];

        t_90[k] = pa_z[k] * dk_27[k];

        t_91[k] = pa_y[k] * dk_34[k];

        t_92[k] = pa_z[k] * dk_28[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, pa_x, pa_y, dk_35, dk_63, \
                         dk_64, dk_65, dk_66, dk_67, dk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * dk_35[k];

        t_94[k] = pa_x[k] * dk_63[k];

        t_95[k] = pa_x[k] * dk_64[k];

        t_96[k] = pa_x[k] * dk_65[k];

        t_97[k] = pa_x[k] * dk_66[k];

        t_98[k] = pa_x[k] * dk_67[k];

        t_99[k] = pa_x[k] * dk_68[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_z, di_21, di_45, di_47, di_48, \
                         dk_70, dk_73, dk_74, fi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_14 * di_45[k]
                   + pa_x[k] * dk_70[k];

        t_101[k] = f_11 * di_21[k]
                   + pb_z[k] * fi_38[k];

        t_102[k] = f_13 * di_47[k]
                   + pa_x[k] * dk_73[k];

        t_103[k] = f_13 * di_48[k]
                   + pa_x[k] * dk_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_x, di_49, di_51, di_52, di_55, \
                         di_56, dk_75, dk_77, dk_78, dk_81, dk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_12 * di_49[k]
                   + pa_x[k] * dk_75[k];

        t_105[k] = f_12 * di_51[k]
                   + pa_x[k] * dk_77[k];

        t_106[k] = f_0 * di_52[k]
                   + pa_x[k] * dk_78[k];

        t_107[k] = f_0 * di_55[k]
                   + pa_x[k] * dk_81[k];

        t_108[k] = f_11 * di_56[k]
                   + pa_x[k] * dk_82[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, pa_x, pb_x, di_59, di_65, \
                         dk_86, dk_88, dk_89, dk_90, dk_91, fi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_11 * di_59[k]
                   + pa_x[k] * dk_86[k];

        t_110[k] = f_15 * di_65[k]
                   + pb_x[k] * fi_43[k];

        t_111[k] = pa_x[k] * dk_88[k];

        t_112[k] = pa_x[k] * dk_89[k];

        t_113[k] = pa_x[k] * dk_90[k];

        t_114[k] = pa_x[k] * dk_91[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_x, pb_x, pb_z, dk_92, dk_93, \
                         dk_95, fh0_21, fh1_23, fi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_x[k] * dk_92[k];

        t_116[k] = pa_x[k] * dk_93[k];

        t_117[k] = pa_x[k] * dk_95[k];

        t_118[k] = f_1 * fh0_21[k]
                   - f_2 * fh1_23[k]
                   + pb_x[k] * fi_44[k];

        t_119[k] = pb_z[k] * fi_44[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_x, pb_z, fh0_22, fh0_23, fh0_24, \
                         fh1_24, fh1_25, fh1_26, fi_46, fi_47, fi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_9 * fh0_22[k]
                   - f_10 * fh1_24[k]
                   + pb_x[k] * fi_46[k];

        t_121[k] = f_9 * fh0_23[k]
                   - f_10 * fh1_25[k]
                   + pb_x[k] * fi_47[k];

        t_122[k] = f_7 * fh0_24[k]
                   - f_8 * fh1_26[k]
                   + pb_x[k] * fi_48[k];

        t_123[k] = pb_z[k] * fi_46[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_x, pb_z, fh0_25, fh0_26, fh0_27, \
                         fh1_27, fh1_28, fh1_29, fi_48, fi_49, fi_50, \
                         fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_7 * fh0_25[k]
                   - f_8 * fh1_27[k]
                   + pb_x[k] * fi_49[k];

        t_125[k] = f_5 * fh0_26[k]
                   - f_6 * fh1_28[k]
                   + pb_x[k] * fi_50[k];

        t_126[k] = pb_z[k] * fi_48[k];

        t_127[k] = f_5 * fh0_27[k]
                   - f_6 * fh1_29[k]
                   + pb_x[k] * fi_51[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_x, pb_z, fh0_28, fh0_29, fh0_31, \
                         fh1_30, fh1_31, fh1_33, fi_50, fi_52, fi_53, \
                         fi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_5 * fh0_28[k]
                   - f_6 * fh1_30[k]
                   + pb_x[k] * fi_52[k];

        t_129[k] = f_3 * fh0_29[k]
                   - f_4 * fh1_31[k]
                   + pb_x[k] * fi_53[k];

        t_130[k] = pb_z[k] * fi_50[k];

        t_131[k] = f_3 * fh0_31[k]
                   - f_4 * fh1_33[k]
                   + pb_x[k] * fi_54[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pb_x, fh0_32, fh0_33, fh1_34, \
                         fh1_35, fi_55, fi_56, fi_57, fi_59, fi_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_3 * fh0_32[k]
                   - f_4 * fh1_34[k]
                   + pb_x[k] * fi_55[k];

        t_133[k] = f_3 * fh0_33[k]
                   - f_4 * fh1_35[k]
                   + pb_x[k] * fi_56[k];

        t_134[k] = pb_x[k] * fi_57[k];

        t_135[k] = pb_x[k] * fi_59[k];

        t_136[k] = pb_x[k] * fi_60[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pb_x, pb_y, pb_z, di_37, fh0_29, \
                         fh1_31, fi_57, fi_58, fi_61, fi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = pb_x[k] * fi_61[k];

        t_138[k] = pb_x[k] * fi_62[k];

        t_139[k] = f_0 * di_37[k]
                   + f_1 * fh0_29[k]
                   - f_2 * fh1_31[k]
                   + pb_y[k] * fi_57[k];

        t_140[k] = pb_z[k] * fi_57[k];

        t_141[k] = f_3 * fh0_29[k]
                   - f_4 * fh1_31[k]
                   + pb_z[k] * fi_58[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_z, fh0_30, fh0_31, fh0_32, fh1_32, fh1_33, \
                         fh1_34, fi_59, fi_60, fi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_5 * fh0_30[k]
                   - f_6 * fh1_32[k]
                   + pb_z[k] * fi_59[k];

        t_143[k] = f_7 * fh0_31[k]
                   - f_8 * fh1_33[k]
                   + pb_z[k] * fi_60[k];

        t_144[k] = f_9 * fh0_32[k]
                   - f_10 * fh1_34[k]
                   + pb_z[k] * fi_61[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, di_22, di_42, \
                         dk_37, dk_39, fh0_33, fh1_35, fi_62, fi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_0 * di_42[k]
                   + pb_y[k] * fi_62[k];

        t_146[k] = f_1 * fh0_33[k]
                   - f_2 * fh1_35[k]
                   + pb_z[k] * fi_62[k];

        t_147[k] = pa_z[k] * dk_37[k];

        t_148[k] = f_15 * di_22[k]
                   + pb_z[k] * fi_63[k];

        t_149[k] = pa_z[k] * dk_39[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_z, di_23, di_25, di_27, dk_40, \
                         dk_41, dk_43, dk_44, dk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_11 * di_23[k]
                   + pa_z[k] * dk_40[k];

        t_151[k] = pa_z[k] * dk_41[k];

        t_152[k] = f_0 * di_25[k]
                   + pa_z[k] * dk_43[k];

        t_153[k] = pa_z[k] * dk_44[k];

        t_154[k] = f_11 * di_27[k]
                   + pa_z[k] * dk_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_z, di_28, di_30, di_31, di_32, \
                         dk_47, dk_48, dk_50, dk_51, dk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * di_28[k]
                   + pa_z[k] * dk_47[k];

        t_156[k] = pa_z[k] * dk_48[k];

        t_157[k] = f_11 * di_30[k]
                   + pa_z[k] * dk_50[k];

        t_158[k] = f_0 * di_31[k]
                   + pa_z[k] * dk_51[k];

        t_159[k] = f_13 * di_32[k]
                   + pa_z[k] * dk_52[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, t_165, pa_z, pb_x, pb_z, di_37, \
                         dk_54, fi_68, fi_69, fi_70, fi_71, fi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_x[k] * fi_69[k];

        t_161[k] = pb_x[k] * fi_70[k];

        t_162[k] = pb_x[k] * fi_71[k];

        t_163[k] = pb_x[k] * fi_72[k];

        t_164[k] = pa_z[k] * dk_54[k];

        t_165[k] = f_15 * di_37[k]
                   + pb_z[k] * fi_68[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, di_38, di_39, di_40, di_41, dk_56, \
                         dk_57, dk_58, dk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_11 * di_38[k]
                   + pa_z[k] * dk_56[k];

        t_167[k] = f_0 * di_39[k]
                   + pa_z[k] * dk_57[k];

        t_168[k] = f_12 * di_40[k]
                   + pa_z[k] * dk_58[k];

        t_169[k] = f_13 * di_41[k]
                   + pa_z[k] * dk_59[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pa_y, pa_z, pb_y, di_42, di_44, \
                         di_46, dk_61, dk_70, dk_72, dk_73, fi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_11 * di_44[k]
                   + pb_y[k] * fi_72[k];

        t_171[k] = f_14 * di_42[k]
                   + pa_z[k] * dk_61[k];

        t_172[k] = pa_y[k] * dk_70[k];

        t_173[k] = pa_y[k] * dk_72[k];

        t_174[k] = f_11 * di_46[k]
                   + pa_y[k] * dk_73[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, pa_y, di_47, di_49, di_50, \
                         dk_74, dk_75, dk_77, dk_78, dk_79, dk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_y[k] * dk_74[k];

        t_176[k] = f_0 * di_47[k]
                   + pa_y[k] * dk_75[k];

        t_177[k] = pa_y[k] * dk_77[k];

        t_178[k] = f_12 * di_49[k]
                   + pa_y[k] * dk_78[k];

        t_179[k] = f_11 * di_50[k]
                   + pa_y[k] * dk_79[k];

        t_180[k] = pa_y[k] * dk_81[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, pa_y, pb_x, di_52, di_53, di_54, \
                         dk_82, dk_83, dk_84, dk_86, fi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * di_52[k]
                   + pa_y[k] * dk_82[k];

        t_182[k] = f_0 * di_53[k]
                   + pa_y[k] * dk_83[k];

        t_183[k] = f_11 * di_54[k]
                   + pa_y[k] * dk_84[k];

        t_184[k] = pa_y[k] * dk_86[k];

        t_185[k] = pb_x[k] * fi_77[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pb_x, pb_z, di_43, di_60, \
                         dk_88, fi_77, fi_78, fi_79, fi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pb_x[k] * fi_78[k];

        t_187[k] = pb_x[k] * fi_79[k];

        t_188[k] = pb_x[k] * fi_80[k];

        t_189[k] = f_14 * di_60[k]
                   + pa_y[k] * dk_88[k];

        t_190[k] = f_11 * di_43[k]
                   + pb_z[k] * fi_77[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, di_61, di_62, di_63, di_64, dk_90, \
                         dk_91, dk_92, dk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * di_61[k]
                   + pa_y[k] * dk_90[k];

        t_192[k] = f_12 * di_62[k]
                   + pa_y[k] * dk_91[k];

        t_193[k] = f_0 * di_63[k]
                   + pa_y[k] * dk_92[k];

        t_194[k] = f_11 * di_64[k]
                   + pa_y[k] * dk_93[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pa_y, pb_x, pb_y, pb_z, di_45, \
                         di_65, dk_95, fh0_38, fh1_41, fi_82, fi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_15 * di_65[k]
                   + pb_y[k] * fi_82[k];

        t_196[k] = pa_y[k] * dk_95[k];

        t_197[k] = f_1 * fh0_38[k]
                   - f_2 * fh1_41[k]
                   + pb_x[k] * fi_83[k];

        t_198[k] = pb_y[k] * fi_83[k];

        t_199[k] = f_0 * di_45[k]
                   + pb_z[k] * fi_83[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pb_x, pb_y, fh0_39, fh0_40, fh0_41, \
                         fh1_42, fh1_43, fh1_44, fi_85, fi_86, fi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_9 * fh0_39[k]
                   - f_10 * fh1_42[k]
                   + pb_x[k] * fi_85[k];

        t_201[k] = f_9 * fh0_40[k]
                   - f_10 * fh1_43[k]
                   + pb_x[k] * fi_86[k];

        t_202[k] = f_7 * fh0_41[k]
                   - f_8 * fh1_44[k]
                   + pb_x[k] * fi_87[k];

        t_203[k] = pb_y[k] * fi_86[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_x, pb_y, fh0_42, fh0_43, fh0_44, \
                         fh1_45, fh1_46, fh1_47, fi_88, fi_89, fi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_7 * fh0_42[k]
                   - f_8 * fh1_45[k]
                   + pb_x[k] * fi_88[k];

        t_205[k] = f_5 * fh0_43[k]
                   - f_6 * fh1_46[k]
                   + pb_x[k] * fi_89[k];

        t_206[k] = f_5 * fh0_44[k]
                   - f_6 * fh1_47[k]
                   + pb_x[k] * fi_90[k];

        t_207[k] = pb_y[k] * fi_88[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pb_x, fh0_45, fh0_46, fh0_47, fh1_48, fh1_49, \
                         fh1_50, fi_91, fi_92, fi_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_5 * fh0_45[k]
                   - f_6 * fh1_48[k]
                   + pb_x[k] * fi_91[k];

        t_209[k] = f_3 * fh0_46[k]
                   - f_4 * fh1_49[k]
                   + pb_x[k] * fi_92[k];

        t_210[k] = f_3 * fh0_47[k]
                   - f_4 * fh1_50[k]
                   + pb_x[k] * fi_93[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pb_x, pb_y, fh0_48, fh0_50, \
                         fh1_51, fh1_53, fi_91, fi_94, fi_95, fi_96, \
                         fi_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_3 * fh0_48[k]
                   - f_4 * fh1_51[k]
                   + pb_x[k] * fi_94[k];

        t_212[k] = pb_y[k] * fi_91[k];

        t_213[k] = f_3 * fh0_50[k]
                   - f_4 * fh1_53[k]
                   + pb_x[k] * fi_95[k];

        t_214[k] = pb_x[k] * fi_96[k];

        t_215[k] = pb_x[k] * fi_97[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pb_x, pb_y, pb_z, di_60, fh0_46, \
                         fh1_49, fi_96, fi_98, fi_99, fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * fi_98[k];

        t_217[k] = pb_x[k] * fi_99[k];

        t_218[k] = pb_x[k] * fi_101[k];

        t_219[k] = f_1 * fh0_46[k]
                   - f_2 * fh1_49[k]
                   + pb_y[k] * fi_96[k];

        t_220[k] = f_0 * di_60[k]
                   + pb_z[k] * fi_96[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pb_y, fh0_47, fh0_48, fh0_49, fh1_50, fh1_51, \
                         fh1_52, fi_97, fi_98, fi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * fh0_47[k]
                   - f_10 * fh1_50[k]
                   + pb_y[k] * fi_97[k];

        t_222[k] = f_7 * fh0_48[k]
                   - f_8 * fh1_51[k]
                   + pb_y[k] * fi_98[k];

        t_223[k] = f_5 * fh0_49[k]
                   - f_6 * fh1_52[k]
                   + pb_y[k] * fi_99[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, pb_y, pb_z, di_65, fh0_50, fh1_53, fi_100, \
                         fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_3 * fh0_50[k]
                   - f_4 * fh1_53[k]
                   + pb_y[k] * fi_100[k];

        t_225[k] = pb_y[k] * fi_101[k];

        t_226[k] = f_0 * di_65[k]
                   + f_1 * fh0_50[k]
                   - f_2 * fh1_53[k]
                   + pb_z[k] * fi_101[k];
    }
}

auto
compute_prim_fk_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 1.0 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_41 = buffer.data(dk + 41);
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
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_71 = buffer.data(dk + 71);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_41 = buffer.data(fh0 + 41);
    const auto *fh0_42 = buffer.data(fh0 + 42);
    const auto *fh0_43 = buffer.data(fh0 + 43);
    const auto *fh0_44 = buffer.data(fh0 + 44);
    const auto *fh0_45 = buffer.data(fh0 + 45);
    const auto *fh0_46 = buffer.data(fh0 + 46);
    const auto *fh0_47 = buffer.data(fh0 + 47);
    const auto *fh0_48 = buffer.data(fh0 + 48);
    const auto *fh0_49 = buffer.data(fh0 + 49);
    const auto *fh0_50 = buffer.data(fh0 + 50);
    const auto *fh0_51 = buffer.data(fh0 + 51);
    const auto *fh0_52 = buffer.data(fh0 + 52);
    const auto *fh0_53 = buffer.data(fh0 + 53);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_36 = buffer.data(fh1 + 36);
    const auto *fh1_44 = buffer.data(fh1 + 44);
    const auto *fh1_45 = buffer.data(fh1 + 45);
    const auto *fh1_46 = buffer.data(fh1 + 46);
    const auto *fh1_47 = buffer.data(fh1 + 47);
    const auto *fh1_48 = buffer.data(fh1 + 48);
    const auto *fh1_49 = buffer.data(fh1 + 49);
    const auto *fh1_50 = buffer.data(fh1 + 50);
    const auto *fh1_51 = buffer.data(fh1 + 51);
    const auto *fh1_52 = buffer.data(fh1 + 52);
    const auto *fh1_53 = buffer.data(fh1 + 53);
    const auto *fh1_54 = buffer.data(fh1 + 54);
    const auto *fh1_55 = buffer.data(fh1 + 55);
    const auto *fh1_56 = buffer.data(fh1 + 56);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, di_0, fh0_0, fh1_0, fi_0, \
                         fi_1, fi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pb_y[k] * fi_0[k];

        t_2[k] = pb_z[k] * fi_0[k];

        t_3[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_4[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, fh0_1, fh0_2, fh0_3, fh1_1, fh1_2, \
                         fh1_3, fi_3, fi_4, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_6[k] = pb_y[k] * fi_4[k];

        t_7[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_8[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, fh0_4, fh0_5, fh1_4, fh1_5, fi_6, \
                         fi_7, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * fh0_4[k]
                 - f_4 * fh1_4[k]
                 + pb_y[k] * fi_6[k];

        t_10[k] = pb_y[k] * fi_7[k];

        t_11[k] = f_7 * fh0_4[k]
                  - f_8 * fh1_4[k]
                  + pb_z[k] * fi_7[k];

        t_12[k] = f_9 * fh0_5[k]
                  - f_10 * fh1_5[k]
                  + pb_y[k] * fi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, fh0_6, fh0_7, fh1_6, fh1_7, fi_9, \
                         fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fh0_6[k]
                  - f_6 * fh1_6[k]
                  + pb_y[k] * fi_9[k];

        t_14[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_15[k] = pb_y[k] * fi_11[k];

        t_16[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, fh0_8, fh0_9, fh0_10, fh1_8, fh1_9, fh1_10, \
                         fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_18[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_9[k]
                  + pb_y[k] * fi_13[k];

        t_19[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_10[k]
                  + pb_y[k] * fi_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, pb_z, dk_0, fh0_11, fh0_12, \
                         fh1_11, fh1_12, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_11[k]
                  + pb_y[k] * fi_15[k];

        t_21[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_12[k]
                  + pb_y[k] * fi_16[k];

        t_22[k] = pb_y[k] * fi_17[k];

        t_23[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_12[k]
                  + pb_z[k] * fi_17[k];

        t_24[k] = pa_y[k] * dk_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, di_1, di_3, di_5, di_8, di_12, \
                         dk_3, dk_5, dk_8, dk_12, dk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * di_1[k]
                  + pa_y[k] * dk_3[k];

        t_26[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_5[k];

        t_27[k] = f_12 * di_5[k]
                  + pa_y[k] * dk_8[k];

        t_28[k] = f_13 * di_8[k]
                  + pa_y[k] * dk_12[k];

        t_29[k] = f_14 * di_12[k]
                  + pa_y[k] * dk_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_z, pb_z, di_0, di_2, di_4, di_7, \
                         dk_0, dk_4, dk_7, dk_11, fi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * dk_0[k];

        t_31[k] = f_15 * di_0[k]
                  + pb_z[k] * fi_20[k];

        t_32[k] = f_11 * di_2[k]
                  + pa_z[k] * dk_4[k];

        t_33[k] = f_0 * di_4[k]
                  + pa_z[k] * dk_7[k];

        t_34[k] = f_12 * di_7[k]
                  + pa_z[k] * dk_11[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, di_11, di_13, di_14, di_15, \
                         di_16, dk_16, dk_18, dk_19, dk_20, dk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * di_11[k]
                  + pa_z[k] * dk_16[k];

        t_36[k] = f_11 * di_13[k]
                  + pa_z[k] * dk_18[k];

        t_37[k] = f_0 * di_14[k]
                  + pa_z[k] * dk_19[k];

        t_38[k] = f_12 * di_15[k]
                  + pa_z[k] * dk_20[k];

        t_39[k] = f_13 * di_16[k]
                  + pa_z[k] * dk_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_z, di_18, di_20, di_22, di_24, \
                         dk_23, dk_24, dk_25, dk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_14 * di_18[k]
                  + pa_z[k] * dk_23[k];

        t_41[k] = f_14 * di_20[k]
                  + pa_x[k] * dk_24[k];

        t_42[k] = f_13 * di_22[k]
                  + pa_x[k] * dk_25[k];

        t_43[k] = f_12 * di_24[k]
                  + pa_x[k] * dk_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_x, pb_z, di_19, di_26, di_29, di_41, \
                         dk_29, dk_32, dk_41, dk_48, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * di_26[k]
                  + pa_x[k] * dk_29[k];

        t_45[k] = f_11 * di_29[k]
                  + pa_x[k] * dk_32[k];

        t_46[k] = pa_x[k] * dk_41[k];

        t_47[k] = f_14 * di_41[k]
                  + pa_x[k] * dk_48[k];

        t_48[k] = f_11 * di_19[k]
                  + pb_z[k] * fi_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_x, di_44, di_46, di_49, di_53, \
                         dk_50, dk_52, dk_55, dk_59, dk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_13 * di_44[k]
                  + pa_x[k] * dk_50[k];

        t_50[k] = f_12 * di_46[k]
                  + pa_x[k] * dk_52[k];

        t_51[k] = f_0 * di_49[k]
                  + pa_x[k] * dk_55[k];

        t_52[k] = f_11 * di_53[k]
                  + pa_x[k] * dk_59[k];

        t_53[k] = pa_x[k] * dk_71[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, fh0_23, fh0_24, fh0_25, fh1_24, fh1_25, \
                         fh1_26, fi_26, fi_27, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * fh0_23[k]
                  - f_2 * fh1_24[k]
                  + pb_x[k] * fi_26[k];

        t_55[k] = f_9 * fh0_24[k]
                  - f_10 * fh1_25[k]
                  + pb_x[k] * fi_27[k];

        t_56[k] = f_9 * fh0_25[k]
                  - f_10 * fh1_26[k]
                  + pb_x[k] * fi_28[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, fh0_26, fh0_27, fh0_28, fh1_27, fh1_28, \
                         fh1_29, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_7 * fh0_26[k]
                  - f_8 * fh1_27[k]
                  + pb_x[k] * fi_29[k];

        t_58[k] = f_7 * fh0_27[k]
                  - f_8 * fh1_28[k]
                  + pb_x[k] * fi_30[k];

        t_59[k] = f_5 * fh0_28[k]
                  - f_6 * fh1_29[k]
                  + pb_x[k] * fi_31[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_x, fh0_29, fh0_30, fh0_31, fh1_30, fh1_31, \
                         fh1_32, fi_32, fi_33, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * fh0_29[k]
                  - f_6 * fh1_30[k]
                  + pb_x[k] * fi_32[k];

        t_61[k] = f_5 * fh0_30[k]
                  - f_6 * fh1_31[k]
                  + pb_x[k] * fi_33[k];

        t_62[k] = f_3 * fh0_31[k]
                  - f_4 * fh1_32[k]
                  + pb_x[k] * fi_34[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_x, fh0_33, fh0_34, fh0_35, fh1_34, fh1_35, \
                         fh1_36, fi_35, fi_36, fi_37, fi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * fh0_33[k]
                  - f_4 * fh1_34[k]
                  + pb_x[k] * fi_35[k];

        t_64[k] = f_3 * fh0_34[k]
                  - f_4 * fh1_35[k]
                  + pb_x[k] * fi_36[k];

        t_65[k] = f_3 * fh0_35[k]
                  - f_4 * fh1_36[k]
                  + pb_x[k] * fi_37[k];

        t_66[k] = pb_x[k] * fi_38[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pb_x, pb_y, di_33, fh0_31, fh1_32, \
                         fi_38, fi_40, fi_41, fi_42, fi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * fi_40[k];

        t_68[k] = pb_x[k] * fi_41[k];

        t_69[k] = pb_x[k] * fi_42[k];

        t_70[k] = pb_x[k] * fi_43[k];

        t_71[k] = f_0 * di_33[k]
                  + f_1 * fh0_31[k]
                  - f_2 * fh1_32[k]
                  + pb_y[k] * fi_38[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_z, fh0_31, fh0_32, fh0_33, fh1_32, fh1_33, \
                         fh1_34, fi_38, fi_39, fi_40, fi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pb_z[k] * fi_38[k];

        t_73[k] = f_3 * fh0_31[k]
                  - f_4 * fh1_32[k]
                  + pb_z[k] * fi_39[k];

        t_74[k] = f_5 * fh0_32[k]
                  - f_6 * fh1_33[k]
                  + pb_z[k] * fi_40[k];

        t_75[k] = f_7 * fh0_33[k]
                  - f_8 * fh1_34[k]
                  + pb_z[k] * fi_41[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, di_21, di_38, dk_26, \
                         fh0_34, fh0_35, fh1_35, fh1_36, fi_42, fi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_9 * fh0_34[k]
                  - f_10 * fh1_35[k]
                  + pb_z[k] * fi_42[k];

        t_77[k] = f_0 * di_38[k]
                  + pb_y[k] * fi_43[k];

        t_78[k] = f_1 * fh0_35[k]
                  - f_2 * fh1_36[k]
                  + pb_z[k] * fi_43[k];

        t_79[k] = f_11 * di_21[k]
                  + pa_z[k] * dk_26[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_z, pb_z, di_23, di_25, di_28, di_33, \
                         dk_28, dk_31, dk_35, dk_41, fi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * di_23[k]
                  + pa_z[k] * dk_28[k];

        t_81[k] = f_12 * di_25[k]
                  + pa_z[k] * dk_31[k];

        t_82[k] = f_13 * di_28[k]
                  + pa_z[k] * dk_35[k];

        t_83[k] = pa_z[k] * dk_41[k];

        t_84[k] = f_15 * di_33[k]
                  + pb_z[k] * fi_44[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_z, di_34, di_35, di_36, di_37, dk_43, \
                         dk_44, dk_45, dk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * di_34[k]
                  + pa_z[k] * dk_43[k];

        t_86[k] = f_0 * di_35[k]
                  + pa_z[k] * dk_44[k];

        t_87[k] = f_12 * di_36[k]
                  + pa_z[k] * dk_45[k];

        t_88[k] = f_13 * di_37[k]
                  + pa_z[k] * dk_46[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, di_38, di_40, di_42, di_43, \
                         dk_47, dk_49, dk_51, fi_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_11 * di_40[k]
                  + pb_y[k] * fi_45[k];

        t_90[k] = f_14 * di_38[k]
                  + pa_z[k] * dk_47[k];

        t_91[k] = f_11 * di_42[k]
                  + pa_y[k] * dk_49[k];

        t_92[k] = f_0 * di_43[k]
                  + pa_y[k] * dk_51[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_z, di_39, di_45, di_47, di_54, \
                         dk_53, dk_56, dk_65, fi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_12 * di_45[k]
                  + pa_y[k] * dk_53[k];

        t_94[k] = f_13 * di_47[k]
                  + pa_y[k] * dk_56[k];

        t_95[k] = f_14 * di_54[k]
                  + pa_y[k] * dk_65[k];

        t_96[k] = f_11 * di_39[k]
                  + pb_z[k] * fi_46[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, di_55, di_56, di_57, di_58, dk_66, \
                         dk_67, dk_68, dk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * di_55[k]
                  + pa_y[k] * dk_66[k];

        t_98[k] = f_12 * di_56[k]
                  + pa_y[k] * dk_67[k];

        t_99[k] = f_0 * di_57[k]
                  + pa_y[k] * dk_68[k];

        t_100[k] = f_11 * di_58[k]
                   + pa_y[k] * dk_69[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pb_x, pb_y, pb_z, di_41, di_59, \
                         dk_71, fh0_41, fh1_44, fi_47, fi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_15 * di_59[k]
                   + pb_y[k] * fi_47[k];

        t_102[k] = pa_y[k] * dk_71[k];

        t_103[k] = f_1 * fh0_41[k]
                   - f_2 * fh1_44[k]
                   + pb_x[k] * fi_48[k];

        t_104[k] = f_0 * di_41[k]
                   + pb_z[k] * fi_48[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pb_x, fh0_42, fh0_43, fh0_44, fh1_45, fh1_46, \
                         fh1_47, fi_49, fi_50, fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_9 * fh0_42[k]
                   - f_10 * fh1_45[k]
                   + pb_x[k] * fi_49[k];

        t_106[k] = f_9 * fh0_43[k]
                   - f_10 * fh1_46[k]
                   + pb_x[k] * fi_50[k];

        t_107[k] = f_7 * fh0_44[k]
                   - f_8 * fh1_47[k]
                   + pb_x[k] * fi_51[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_x, fh0_45, fh0_46, fh0_47, fh1_48, fh1_49, \
                         fh1_50, fi_52, fi_53, fi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_7 * fh0_45[k]
                   - f_8 * fh1_48[k]
                   + pb_x[k] * fi_52[k];

        t_109[k] = f_5 * fh0_46[k]
                   - f_6 * fh1_49[k]
                   + pb_x[k] * fi_53[k];

        t_110[k] = f_5 * fh0_47[k]
                   - f_6 * fh1_50[k]
                   + pb_x[k] * fi_54[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, fh0_48, fh0_49, fh0_50, fh1_51, fh1_52, \
                         fh1_53, fi_55, fi_56, fi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * fh0_48[k]
                   - f_6 * fh1_51[k]
                   + pb_x[k] * fi_55[k];

        t_112[k] = f_3 * fh0_49[k]
                   - f_4 * fh1_52[k]
                   + pb_x[k] * fi_56[k];

        t_113[k] = f_3 * fh0_50[k]
                   - f_4 * fh1_53[k]
                   + pb_x[k] * fi_57[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pb_x, fh0_51, fh0_53, fh1_54, \
                         fh1_56, fi_58, fi_59, fi_60, fi_61, fi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * fh0_51[k]
                   - f_4 * fh1_54[k]
                   + pb_x[k] * fi_58[k];

        t_115[k] = f_3 * fh0_53[k]
                   - f_4 * fh1_56[k]
                   + pb_x[k] * fi_59[k];

        t_116[k] = pb_x[k] * fi_60[k];

        t_117[k] = pb_x[k] * fi_61[k];

        t_118[k] = pb_x[k] * fi_62[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, pb_z, di_54, fh0_49, fh1_52, \
                         fi_60, fi_63, fi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_x[k] * fi_63[k];

        t_120[k] = pb_x[k] * fi_65[k];

        t_121[k] = f_1 * fh0_49[k]
                   - f_2 * fh1_52[k]
                   + pb_y[k] * fi_60[k];

        t_122[k] = f_0 * di_54[k]
                   + pb_z[k] * fi_60[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_y, fh0_50, fh0_51, fh0_52, fh1_53, fh1_54, \
                         fh1_55, fi_61, fi_62, fi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_9 * fh0_50[k]
                   - f_10 * fh1_53[k]
                   + pb_y[k] * fi_61[k];

        t_124[k] = f_7 * fh0_51[k]
                   - f_8 * fh1_54[k]
                   + pb_y[k] * fi_62[k];

        t_125[k] = f_5 * fh0_52[k]
                   - f_6 * fh1_55[k]
                   + pb_y[k] * fi_63[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pb_y, pb_z, di_59, fh0_53, fh1_56, fi_64, \
                         fi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_3 * fh0_53[k]
                   - f_4 * fh1_56[k]
                   + pb_y[k] * fi_64[k];

        t_127[k] = pb_y[k] * fi_65[k];

        t_128[k] = f_0 * di_59[k]
                   + f_1 * fh0_53[k]
                   - f_2 * fh1_56[k]
                   + pb_z[k] * fi_65[k];
    }
}

auto
compute_prim_fk_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_8 = buffer.data(fh0 + 8);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_8 = buffer.data(fh1 + 8);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_1, dk_1, dk_2, fh0_5, fh1_5, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_1[k]
                 + f_1 * fh0_5[k]
                 - f_2 * fh1_5[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_2, fh0_8, fh1_8, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_2[k]
                 + f_1 * fh0_8[k]
                 - f_2 * fh1_8[k]
                 + pb_z[k] * fi_8[k];
    }
}

auto
compute_prim_fk_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_8 = buffer.data(fh0 + 8);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_44 = buffer.data(fh1 + 44);
    const auto *fh1_77 = buffer.data(fh1 + 77);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_103 = buffer.data(fi + 103);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_1, dk_1, dk_2, fh0_5, fh1_44, \
                         fi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_1[k]
                 + f_1 * fh0_5[k]
                 - f_2 * fh1_44[k]
                 + pb_y[k] * fi_57[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_2, fh0_8, fh1_77, fi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_2[k]
                 + f_1 * fh0_8[k]
                 - f_2 * fh1_77[k]
                 + pb_z[k] * fi_103[k];
    }
}

auto
compute_prim_fk_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_35 = buffer.data(di + 35);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_13 = buffer.data(fh0 + 13);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_36 = buffer.data(fh0 + 36);
    const auto *fh0_37 = buffer.data(fh0 + 37);
    const auto *fh0_38 = buffer.data(fh0 + 38);
    const auto *fh0_40 = buffer.data(fh0 + 40);
    const auto *fh0_41 = buffer.data(fh0 + 41);
    const auto *fh0_42 = buffer.data(fh0 + 42);
    const auto *fh0_43 = buffer.data(fh0 + 43);
    const auto *fh0_44 = buffer.data(fh0 + 44);
    const auto *fh0_45 = buffer.data(fh0 + 45);
    const auto *fh0_46 = buffer.data(fh0 + 46);
    const auto *fh0_47 = buffer.data(fh0 + 47);
    const auto *fh0_48 = buffer.data(fh0 + 48);
    const auto *fh0_63 = buffer.data(fh0 + 63);
    const auto *fh0_65 = buffer.data(fh0 + 65);
    const auto *fh0_66 = buffer.data(fh0 + 66);
    const auto *fh0_67 = buffer.data(fh0 + 67);
    const auto *fh0_69 = buffer.data(fh0 + 69);
    const auto *fh0_70 = buffer.data(fh0 + 70);
    const auto *fh0_71 = buffer.data(fh0 + 71);
    const auto *fh0_72 = buffer.data(fh0 + 72);
    const auto *fh0_73 = buffer.data(fh0 + 73);
    const auto *fh0_74 = buffer.data(fh0 + 74);
    const auto *fh0_75 = buffer.data(fh0 + 75);
    const auto *fh0_76 = buffer.data(fh0 + 76);
    const auto *fh0_77 = buffer.data(fh0 + 77);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_13 = buffer.data(fh1 + 13);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_36 = buffer.data(fh1 + 36);
    const auto *fh1_37 = buffer.data(fh1 + 37);
    const auto *fh1_38 = buffer.data(fh1 + 38);
    const auto *fh1_39 = buffer.data(fh1 + 39);
    const auto *fh1_40 = buffer.data(fh1 + 40);
    const auto *fh1_41 = buffer.data(fh1 + 41);
    const auto *fh1_49 = buffer.data(fh1 + 49);
    const auto *fh1_51 = buffer.data(fh1 + 51);
    const auto *fh1_52 = buffer.data(fh1 + 52);
    const auto *fh1_53 = buffer.data(fh1 + 53);
    const auto *fh1_54 = buffer.data(fh1 + 54);
    const auto *fh1_55 = buffer.data(fh1 + 55);
    const auto *fh1_56 = buffer.data(fh1 + 56);
    const auto *fh1_57 = buffer.data(fh1 + 57);
    const auto *fh1_58 = buffer.data(fh1 + 58);
    const auto *fh1_59 = buffer.data(fh1 + 59);
    const auto *fh1_60 = buffer.data(fh1 + 60);
    const auto *fh1_61 = buffer.data(fh1 + 61);
    const auto *fh1_62 = buffer.data(fh1 + 62);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, di_0, fh0_0, fh0_1, fh1_0, \
                         fh1_1, fi_0, fi_1, fi_2, fi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_2[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];

        t_3[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, fh0_2, fh0_3, fh0_4, fh1_2, fh1_3, \
                         fh1_4, fi_4, fi_5, fi_6, fi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_5[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = f_3 * fh0_4[k]
                 - f_4 * fh1_4[k]
                 + pb_y[k] * fi_6[k];

        t_7[k] = f_7 * fh0_4[k]
                 - f_8 * fh1_4[k]
                 + pb_z[k] * fi_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, fh0_5, fh0_7, fh0_8, fh1_5, fh1_6, \
                         fh1_7, fi_8, fi_9, fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * fh0_5[k]
                 - f_10 * fh1_5[k]
                 + pb_y[k] * fi_8[k];

        t_9[k] = f_5 * fh0_7[k]
                 - f_6 * fh1_6[k]
                 + pb_y[k] * fi_9[k];

        t_10[k] = f_3 * fh0_8[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_11[k] = f_9 * fh0_8[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, fh0_9, fh0_11, fh0_12, fh1_8, fh1_10, fh1_11, \
                         fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * fh0_9[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_13[k] = f_9 * fh0_11[k]
                  - f_10 * fh1_10[k]
                  + pb_y[k] * fi_13[k];

        t_14[k] = f_7 * fh0_12[k]
                  - f_8 * fh1_11[k]
                  + pb_y[k] * fi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pb_y, pb_z, dk_0, fh0_13, fh0_14, \
                         fh1_12, fh1_13, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * fh0_13[k]
                  - f_6 * fh1_12[k]
                  + pb_y[k] * fi_15[k];

        t_16[k] = f_3 * fh0_14[k]
                  - f_4 * fh1_13[k]
                  + pb_y[k] * fi_16[k];

        t_17[k] = f_1 * fh0_14[k]
                  - f_2 * fh1_13[k]
                  + pb_z[k] * fi_17[k];

        t_18[k] = pa_y[k] * dk_0[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pb_x, dk_0, dk_1, dk_2, fh0_34, \
                         fh1_28, fi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * dk_0[k];

        t_20[k] = pa_x[k] * dk_1[k];

        t_21[k] = pa_x[k] * dk_2[k];

        t_22[k] = f_1 * fh0_34[k]
                  - f_2 * fh1_28[k]
                  + pb_x[k] * fi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, fh0_36, fh0_37, fh0_38, fh1_30, fh1_31, \
                         fh1_32, fi_23, fi_24, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_9 * fh0_36[k]
                  - f_10 * fh1_30[k]
                  + pb_x[k] * fi_23[k];

        t_24[k] = f_9 * fh0_37[k]
                  - f_10 * fh1_31[k]
                  + pb_x[k] * fi_24[k];

        t_25[k] = f_7 * fh0_38[k]
                  - f_8 * fh1_32[k]
                  + pb_x[k] * fi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_x, fh0_40, fh0_41, fh0_42, fh1_33, fh1_34, \
                         fh1_35, fi_26, fi_27, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * fh0_40[k]
                  - f_8 * fh1_33[k]
                  + pb_x[k] * fi_26[k];

        t_27[k] = f_5 * fh0_41[k]
                  - f_6 * fh1_34[k]
                  + pb_x[k] * fi_27[k];

        t_28[k] = f_5 * fh0_42[k]
                  - f_6 * fh1_35[k]
                  + pb_x[k] * fi_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, fh0_43, fh0_44, fh0_46, fh1_36, fh1_37, \
                         fh1_39, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_5 * fh0_43[k]
                  - f_6 * fh1_36[k]
                  + pb_x[k] * fi_29[k];

        t_30[k] = f_3 * fh0_44[k]
                  - f_4 * fh1_37[k]
                  + pb_x[k] * fi_30[k];

        t_31[k] = f_3 * fh0_46[k]
                  - f_4 * fh1_39[k]
                  + pb_x[k] * fi_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, di_19, fh0_44, fh0_47, fh0_48, fh1_37, \
                         fh1_40, fh1_41, fi_32, fi_33, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * fh0_47[k]
                  - f_4 * fh1_40[k]
                  + pb_x[k] * fi_32[k];

        t_33[k] = f_3 * fh0_48[k]
                  - f_4 * fh1_41[k]
                  + pb_x[k] * fi_33[k];

        t_34[k] = f_0 * di_19[k]
                  + f_1 * fh0_44[k]
                  - f_2 * fh1_37[k]
                  + pb_y[k] * fi_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_z, fh0_44, fh0_45, fh0_46, fh1_37, fh1_38, \
                         fh1_39, fi_35, fi_36, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * fh0_44[k]
                  - f_4 * fh1_37[k]
                  + pb_z[k] * fi_35[k];

        t_36[k] = f_5 * fh0_45[k]
                  - f_6 * fh1_38[k]
                  + pb_z[k] * fi_36[k];

        t_37[k] = f_7 * fh0_46[k]
                  - f_8 * fh1_39[k]
                  + pb_z[k] * fi_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pa_z, pb_z, dk_1, dk_2, fh0_47, fh0_48, \
                         fh1_40, fh1_41, fi_38, fi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * fh0_47[k]
                  - f_10 * fh1_40[k]
                  + pb_z[k] * fi_38[k];

        t_39[k] = f_1 * fh0_48[k]
                  - f_2 * fh1_41[k]
                  + pb_z[k] * fi_39[k];

        t_40[k] = pa_z[k] * dk_1[k];

        t_41[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, fh0_63, fh0_65, fh0_66, fh1_49, fh1_51, \
                         fh1_52, fi_42, fi_43, fi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * fh0_63[k]
                  - f_2 * fh1_49[k]
                  + pb_x[k] * fi_42[k];

        t_43[k] = f_9 * fh0_65[k]
                  - f_10 * fh1_51[k]
                  + pb_x[k] * fi_43[k];

        t_44[k] = f_9 * fh0_66[k]
                  - f_10 * fh1_52[k]
                  + pb_x[k] * fi_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, fh0_67, fh0_69, fh0_70, fh1_53, fh1_54, \
                         fh1_55, fi_45, fi_46, fi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_7 * fh0_67[k]
                  - f_8 * fh1_53[k]
                  + pb_x[k] * fi_45[k];

        t_46[k] = f_7 * fh0_69[k]
                  - f_8 * fh1_54[k]
                  + pb_x[k] * fi_46[k];

        t_47[k] = f_5 * fh0_70[k]
                  - f_6 * fh1_55[k]
                  + pb_x[k] * fi_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, fh0_71, fh0_72, fh0_73, fh1_56, fh1_57, \
                         fh1_58, fi_48, fi_49, fi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * fh0_71[k]
                  - f_6 * fh1_56[k]
                  + pb_x[k] * fi_48[k];

        t_49[k] = f_5 * fh0_72[k]
                  - f_6 * fh1_57[k]
                  + pb_x[k] * fi_49[k];

        t_50[k] = f_3 * fh0_73[k]
                  - f_4 * fh1_58[k]
                  + pb_x[k] * fi_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, fh0_74, fh0_75, fh0_77, fh1_59, fh1_60, \
                         fh1_62, fi_51, fi_52, fi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fh0_74[k]
                  - f_4 * fh1_59[k]
                  + pb_x[k] * fi_51[k];

        t_52[k] = f_3 * fh0_75[k]
                  - f_4 * fh1_60[k]
                  + pb_x[k] * fi_52[k];

        t_53[k] = f_3 * fh0_77[k]
                  - f_4 * fh1_62[k]
                  + pb_x[k] * fi_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, fh0_73, fh0_74, fh0_75, fh1_58, fh1_59, \
                         fh1_60, fi_54, fi_55, fi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * fh0_73[k]
                  - f_2 * fh1_58[k]
                  + pb_y[k] * fi_54[k];

        t_55[k] = f_9 * fh0_74[k]
                  - f_10 * fh1_59[k]
                  + pb_y[k] * fi_55[k];

        t_56[k] = f_7 * fh0_75[k]
                  - f_8 * fh1_60[k]
                  + pb_y[k] * fi_56[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_y, pb_z, di_35, fh0_76, fh0_77, fh1_61, fh1_62, \
                         fi_57, fi_58, fi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * fh0_76[k]
                  - f_6 * fh1_61[k]
                  + pb_y[k] * fi_57[k];

        t_58[k] = f_3 * fh0_77[k]
                  - f_4 * fh1_62[k]
                  + pb_y[k] * fi_58[k];

        t_59[k] = f_0 * di_35[k]
                  + f_1 * fh0_77[k]
                  - f_2 * fh1_62[k]
                  + pb_z[k] * fi_59[k];
    }
}

auto
compute_prim_fk_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_8 = buffer.data(fh0 + 8);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_53 = buffer.data(fh1 + 53);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_101 = buffer.data(fi + 101);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_1, dk_1, dk_2, fh0_5, fh1_31, \
                         fi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_1[k]
                 + f_1 * fh0_5[k]
                 - f_2 * fh1_31[k]
                 + pb_y[k] * fi_58[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_2, fh0_8, fh1_53, fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_2[k]
                 + f_1 * fh0_8[k]
                 - f_2 * fh1_53[k]
                 + pb_z[k] * fi_101[k];
    }
}

auto
compute_prim_fk_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di, const size_t dk,
                                     const size_t fh0, const size_t fh1, const size_t fi,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 1.0 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_41 = buffer.data(fh0 + 41);
    const auto *fh0_42 = buffer.data(fh0 + 42);
    const auto *fh0_43 = buffer.data(fh0 + 43);
    const auto *fh0_44 = buffer.data(fh0 + 44);
    const auto *fh0_45 = buffer.data(fh0 + 45);
    const auto *fh0_46 = buffer.data(fh0 + 46);
    const auto *fh0_47 = buffer.data(fh0 + 47);
    const auto *fh0_48 = buffer.data(fh0 + 48);
    const auto *fh0_49 = buffer.data(fh0 + 49);
    const auto *fh0_50 = buffer.data(fh0 + 50);
    const auto *fh0_51 = buffer.data(fh0 + 51);
    const auto *fh0_52 = buffer.data(fh0 + 52);
    const auto *fh0_53 = buffer.data(fh0 + 53);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_13 = buffer.data(fh1 + 13);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_36 = buffer.data(fh1 + 36);
    const auto *fh1_37 = buffer.data(fh1 + 37);
    const auto *fh1_38 = buffer.data(fh1 + 38);
    const auto *fh1_46 = buffer.data(fh1 + 46);
    const auto *fh1_48 = buffer.data(fh1 + 48);
    const auto *fh1_49 = buffer.data(fh1 + 49);
    const auto *fh1_50 = buffer.data(fh1 + 50);
    const auto *fh1_51 = buffer.data(fh1 + 51);
    const auto *fh1_52 = buffer.data(fh1 + 52);
    const auto *fh1_53 = buffer.data(fh1 + 53);
    const auto *fh1_54 = buffer.data(fh1 + 54);
    const auto *fh1_55 = buffer.data(fh1 + 55);
    const auto *fh1_56 = buffer.data(fh1 + 56);
    const auto *fh1_57 = buffer.data(fh1 + 57);
    const auto *fh1_58 = buffer.data(fh1 + 58);
    const auto *fh1_59 = buffer.data(fh1 + 59);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, di_0, fh0_0, fh1_0, fi_0, \
                         fi_1, fi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pb_y[k] * fi_0[k];

        t_2[k] = pb_z[k] * fi_0[k];

        t_3[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_4[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fh0_1, fh0_2, fh0_3, fh1_1, \
                         fh1_2, fh1_3, fi_3, fi_4, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_6[k] = pb_z[k] * fi_3[k];

        t_7[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_8[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];

        t_9[k] = pb_z[k] * fi_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, fh0_4, fh0_5, fh1_4, fh1_5, fi_6, \
                         fi_7, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fh0_4[k]
                  - f_4 * fh1_4[k]
                  + pb_y[k] * fi_6[k];

        t_11[k] = f_7 * fh0_4[k]
                  - f_8 * fh1_4[k]
                  + pb_z[k] * fi_7[k];

        t_12[k] = f_9 * fh0_5[k]
                  - f_10 * fh1_5[k]
                  + pb_y[k] * fi_8[k];

        t_13[k] = pb_z[k] * fi_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, pb_z, fh0_6, fh0_7, fh0_8, fh1_6, \
                         fh1_7, fh1_8, fi_9, fi_10, fi_11, fi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * fh0_6[k]
                  - f_6 * fh1_6[k]
                  + pb_y[k] * fi_9[k];

        t_15[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_16[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];

        t_17[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, fh0_9, fh0_10, fh0_11, fh1_10, \
                         fh1_11, fh1_12, fi_12, fi_14, fi_15, fi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_z[k] * fi_12[k];

        t_19[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_10[k]
                  + pb_y[k] * fi_14[k];

        t_20[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_11[k]
                  + pb_y[k] * fi_15[k];

        t_21[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_12[k]
                  + pb_y[k] * fi_16[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pb_y, pb_z, di_1, dk_0, dk_1, fh0_12, \
                         fh1_13, fi_17, fi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_13[k]
                  + pb_y[k] * fi_17[k];

        t_23[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_13[k]
                  + pb_z[k] * fi_18[k];

        t_24[k] = pa_y[k] * dk_0[k];

        t_25[k] = f_11 * di_1[k]
                  + pa_y[k] * dk_1[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, di_3, di_5, di_7, di_9, \
                         dk_0, dk_3, dk_5, dk_7, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_3[k];

        t_27[k] = f_12 * di_5[k]
                  + pa_y[k] * dk_5[k];

        t_28[k] = f_13 * di_7[k]
                  + pa_y[k] * dk_7[k];

        t_29[k] = f_14 * di_9[k]
                  + pa_y[k] * dk_9[k];

        t_30[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_z, pb_z, di_0, di_2, di_4, di_6, dk_2, \
                         dk_4, dk_6, fi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_15 * di_0[k]
                  + pb_z[k] * fi_21[k];

        t_32[k] = f_11 * di_2[k]
                  + pa_z[k] * dk_2[k];

        t_33[k] = f_0 * di_4[k]
                  + pa_z[k] * dk_4[k];

        t_34[k] = f_12 * di_6[k]
                  + pa_z[k] * dk_6[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_z, di_8, di_10, di_11, di_12, di_13, \
                         dk_8, dk_10, dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * di_8[k]
                  + pa_z[k] * dk_8[k];

        t_36[k] = f_11 * di_10[k]
                  + pa_z[k] * dk_10[k];

        t_37[k] = f_0 * di_11[k]
                  + pa_z[k] * dk_11[k];

        t_38[k] = f_12 * di_12[k]
                  + pa_z[k] * dk_12[k];

        t_39[k] = f_13 * di_13[k]
                  + pa_z[k] * dk_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_z, di_14, di_16, di_18, di_20, \
                         dk_14, dk_15, dk_16, dk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_14 * di_14[k]
                  + pa_z[k] * dk_14[k];

        t_41[k] = f_14 * di_16[k]
                  + pa_x[k] * dk_15[k];

        t_42[k] = f_13 * di_18[k]
                  + pa_x[k] * dk_16[k];

        t_43[k] = f_12 * di_20[k]
                  + pa_x[k] * dk_18[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_x, pb_x, di_22, di_24, di_25, di_33, \
                         dk_20, dk_22, dk_24, dk_30, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * di_22[k]
                  + pa_x[k] * dk_20[k];

        t_45[k] = f_11 * di_24[k]
                  + pa_x[k] * dk_22[k];

        t_46[k] = f_15 * di_25[k]
                  + pb_x[k] * fi_28[k];

        t_47[k] = pa_x[k] * dk_24[k];

        t_48[k] = f_14 * di_33[k]
                  + pa_x[k] * dk_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pb_z, di_15, di_36, di_38, di_40, \
                         dk_32, dk_34, dk_36, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * di_15[k]
                  + pb_z[k] * fi_29[k];

        t_50[k] = f_13 * di_36[k]
                  + pa_x[k] * dk_32[k];

        t_51[k] = f_12 * di_38[k]
                  + pa_x[k] * dk_34[k];

        t_52[k] = f_0 * di_40[k]
                  + pa_x[k] * dk_36[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_x, pb_x, pb_z, di_41, di_47, dk_38, \
                         dk_44, fh0_23, fh1_25, fi_34, fi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * di_41[k]
                  + pa_x[k] * dk_38[k];

        t_54[k] = f_15 * di_47[k]
                  + pb_x[k] * fi_34[k];

        t_55[k] = pa_x[k] * dk_44[k];

        t_56[k] = f_1 * fh0_23[k]
                  - f_2 * fh1_25[k]
                  + pb_x[k] * fi_35[k];

        t_57[k] = pb_z[k] * fi_35[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_z, fh0_24, fh0_25, fh0_26, fh1_27, \
                         fh1_28, fh1_29, fi_37, fi_38, fi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_9 * fh0_24[k]
                  - f_10 * fh1_27[k]
                  + pb_x[k] * fi_37[k];

        t_59[k] = f_9 * fh0_25[k]
                  - f_10 * fh1_28[k]
                  + pb_x[k] * fi_38[k];

        t_60[k] = f_7 * fh0_26[k]
                  - f_8 * fh1_29[k]
                  + pb_x[k] * fi_39[k];

        t_61[k] = pb_z[k] * fi_37[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_x, pb_z, fh0_27, fh0_28, fh0_29, fh1_30, \
                         fh1_31, fh1_32, fi_39, fi_40, fi_41, fi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_7 * fh0_27[k]
                  - f_8 * fh1_30[k]
                  + pb_x[k] * fi_40[k];

        t_63[k] = f_5 * fh0_28[k]
                  - f_6 * fh1_31[k]
                  + pb_x[k] * fi_41[k];

        t_64[k] = pb_z[k] * fi_39[k];

        t_65[k] = f_5 * fh0_29[k]
                  - f_6 * fh1_32[k]
                  + pb_x[k] * fi_42[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pb_z, fh0_30, fh0_31, fh0_33, fh1_33, \
                         fh1_34, fh1_36, fi_41, fi_43, fi_44, fi_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * fh0_30[k]
                  - f_6 * fh1_33[k]
                  + pb_x[k] * fi_43[k];

        t_67[k] = f_3 * fh0_31[k]
                  - f_4 * fh1_34[k]
                  + pb_x[k] * fi_44[k];

        t_68[k] = pb_z[k] * fi_41[k];

        t_69[k] = f_3 * fh0_33[k]
                  - f_4 * fh1_36[k]
                  + pb_x[k] * fi_45[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pb_y, di_25, fh0_31, fh0_34, fh0_35, \
                         fh1_34, fh1_37, fh1_38, fi_46, fi_47, fi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * fh0_34[k]
                  - f_4 * fh1_37[k]
                  + pb_x[k] * fi_46[k];

        t_71[k] = f_3 * fh0_35[k]
                  - f_4 * fh1_38[k]
                  + pb_x[k] * fi_47[k];

        t_72[k] = pb_x[k] * fi_48[k];

        t_73[k] = f_0 * di_25[k]
                  + f_1 * fh0_31[k]
                  - f_2 * fh1_34[k]
                  + pb_y[k] * fi_48[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, fh0_31, fh0_32, fh0_33, fh1_34, fh1_35, \
                         fh1_36, fi_48, fi_49, fi_50, fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_z[k] * fi_48[k];

        t_75[k] = f_3 * fh0_31[k]
                  - f_4 * fh1_34[k]
                  + pb_z[k] * fi_49[k];

        t_76[k] = f_5 * fh0_32[k]
                  - f_6 * fh1_35[k]
                  + pb_z[k] * fi_50[k];

        t_77[k] = f_7 * fh0_33[k]
                  - f_8 * fh1_36[k]
                  + pb_z[k] * fi_51[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_z, pb_y, pb_z, di_17, di_30, dk_17, \
                         fh0_34, fh0_35, fh1_37, fh1_38, fi_52, fi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * fh0_34[k]
                  - f_10 * fh1_37[k]
                  + pb_z[k] * fi_52[k];

        t_79[k] = f_0 * di_30[k]
                  + pb_y[k] * fi_53[k];

        t_80[k] = f_1 * fh0_35[k]
                  - f_2 * fh1_38[k]
                  + pb_z[k] * fi_53[k];

        t_81[k] = f_11 * di_17[k]
                  + pa_z[k] * dk_17[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_z, pb_z, di_19, di_21, di_23, di_25, \
                         dk_19, dk_21, dk_23, dk_24, fi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * di_19[k]
                  + pa_z[k] * dk_19[k];

        t_83[k] = f_12 * di_21[k]
                  + pa_z[k] * dk_21[k];

        t_84[k] = f_13 * di_23[k]
                  + pa_z[k] * dk_23[k];

        t_85[k] = pa_z[k] * dk_24[k];

        t_86[k] = f_15 * di_25[k]
                  + pb_z[k] * fi_54[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_z, di_26, di_27, di_28, di_29, dk_25, \
                         dk_26, dk_27, dk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_11 * di_26[k]
                  + pa_z[k] * dk_25[k];

        t_88[k] = f_0 * di_27[k]
                  + pa_z[k] * dk_26[k];

        t_89[k] = f_12 * di_28[k]
                  + pa_z[k] * dk_27[k];

        t_90[k] = f_13 * di_29[k]
                  + pa_z[k] * dk_28[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pa_z, pb_y, di_30, di_32, di_34, di_35, \
                         dk_29, dk_31, dk_33, fi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_11 * di_32[k]
                  + pb_y[k] * fi_55[k];

        t_92[k] = f_14 * di_30[k]
                  + pa_z[k] * dk_29[k];

        t_93[k] = f_11 * di_34[k]
                  + pa_y[k] * dk_31[k];

        t_94[k] = f_0 * di_35[k]
                  + pa_y[k] * dk_33[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_y, pb_z, di_31, di_37, di_39, di_42, \
                         dk_35, dk_37, dk_39, fi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_12 * di_37[k]
                  + pa_y[k] * dk_35[k];

        t_96[k] = f_13 * di_39[k]
                  + pa_y[k] * dk_37[k];

        t_97[k] = f_14 * di_42[k]
                  + pa_y[k] * dk_39[k];

        t_98[k] = f_11 * di_31[k]
                  + pb_z[k] * fi_56[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_y, di_43, di_44, di_45, di_46, dk_40, \
                         dk_41, dk_42, dk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * di_43[k]
                  + pa_y[k] * dk_40[k];

        t_100[k] = f_12 * di_44[k]
                   + pa_y[k] * dk_41[k];

        t_101[k] = f_0 * di_45[k]
                   + pa_y[k] * dk_42[k];

        t_102[k] = f_11 * di_46[k]
                   + pa_y[k] * dk_43[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_y, pb_x, pb_y, pb_z, di_33, \
                         di_47, dk_44, fh0_41, fh1_46, fi_61, fi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_15 * di_47[k]
                   + pb_y[k] * fi_61[k];

        t_104[k] = pa_y[k] * dk_44[k];

        t_105[k] = f_1 * fh0_41[k]
                   - f_2 * fh1_46[k]
                   + pb_x[k] * fi_62[k];

        t_106[k] = pb_y[k] * fi_62[k];

        t_107[k] = f_0 * di_33[k]
                   + pb_z[k] * fi_62[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_x, pb_y, fh0_42, fh0_43, fh0_44, \
                         fh1_48, fh1_49, fh1_50, fi_64, fi_65, fi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_9 * fh0_42[k]
                   - f_10 * fh1_48[k]
                   + pb_x[k] * fi_64[k];

        t_109[k] = f_9 * fh0_43[k]
                   - f_10 * fh1_49[k]
                   + pb_x[k] * fi_65[k];

        t_110[k] = f_7 * fh0_44[k]
                   - f_8 * fh1_50[k]
                   + pb_x[k] * fi_66[k];

        t_111[k] = pb_y[k] * fi_65[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_x, pb_y, fh0_45, fh0_46, fh0_47, \
                         fh1_51, fh1_52, fh1_53, fi_67, fi_68, fi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_7 * fh0_45[k]
                   - f_8 * fh1_51[k]
                   + pb_x[k] * fi_67[k];

        t_113[k] = f_5 * fh0_46[k]
                   - f_6 * fh1_52[k]
                   + pb_x[k] * fi_68[k];

        t_114[k] = f_5 * fh0_47[k]
                   - f_6 * fh1_53[k]
                   + pb_x[k] * fi_69[k];

        t_115[k] = pb_y[k] * fi_67[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_x, fh0_48, fh0_49, fh0_50, fh1_54, fh1_55, \
                         fh1_56, fi_70, fi_71, fi_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_5 * fh0_48[k]
                   - f_6 * fh1_54[k]
                   + pb_x[k] * fi_70[k];

        t_117[k] = f_3 * fh0_49[k]
                   - f_4 * fh1_55[k]
                   + pb_x[k] * fi_71[k];

        t_118[k] = f_3 * fh0_50[k]
                   - f_4 * fh1_56[k]
                   + pb_x[k] * fi_72[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, fh0_51, fh0_53, fh1_57, \
                         fh1_59, fi_70, fi_73, fi_74, fi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * fh0_51[k]
                   - f_4 * fh1_57[k]
                   + pb_x[k] * fi_73[k];

        t_120[k] = pb_y[k] * fi_70[k];

        t_121[k] = f_3 * fh0_53[k]
                   - f_4 * fh1_59[k]
                   + pb_x[k] * fi_74[k];

        t_122[k] = pb_x[k] * fi_80[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_y, pb_z, di_42, fh0_49, fh0_50, \
                         fh0_51, fh1_55, fh1_56, fh1_57, fi_75, fi_76, \
                         fi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_1 * fh0_49[k]
                   - f_2 * fh1_55[k]
                   + pb_y[k] * fi_75[k];

        t_124[k] = f_0 * di_42[k]
                   + pb_z[k] * fi_75[k];

        t_125[k] = f_9 * fh0_50[k]
                   - f_10 * fh1_56[k]
                   + pb_y[k] * fi_76[k];

        t_126[k] = f_7 * fh0_51[k]
                   - f_8 * fh1_57[k]
                   + pb_y[k] * fi_77[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pb_y, pb_z, di_47, fh0_52, fh0_53, \
                         fh1_58, fh1_59, fi_78, fi_79, fi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_5 * fh0_52[k]
                   - f_6 * fh1_58[k]
                   + pb_y[k] * fi_78[k];

        t_128[k] = f_3 * fh0_53[k]
                   - f_4 * fh1_59[k]
                   + pb_y[k] * fi_79[k];

        t_129[k] = pb_y[k] * fi_80[k];

        t_130[k] = f_0 * di_47[k]
                   + f_1 * fh0_53[k]
                   - f_2 * fh1_59[k]
                   + pb_z[k] * fi_80[k];
    }
}

auto
compute_prim_fk_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_44 = buffer.data(di + 44);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_13 = buffer.data(fh0 + 13);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_36 = buffer.data(fh0 + 36);
    const auto *fh0_37 = buffer.data(fh0 + 37);
    const auto *fh0_38 = buffer.data(fh0 + 38);
    const auto *fh0_46 = buffer.data(fh0 + 46);
    const auto *fh0_48 = buffer.data(fh0 + 48);
    const auto *fh0_49 = buffer.data(fh0 + 49);
    const auto *fh0_50 = buffer.data(fh0 + 50);
    const auto *fh0_51 = buffer.data(fh0 + 51);
    const auto *fh0_52 = buffer.data(fh0 + 52);
    const auto *fh0_53 = buffer.data(fh0 + 53);
    const auto *fh0_54 = buffer.data(fh0 + 54);
    const auto *fh0_55 = buffer.data(fh0 + 55);
    const auto *fh0_56 = buffer.data(fh0 + 56);
    const auto *fh0_57 = buffer.data(fh0 + 57);
    const auto *fh0_58 = buffer.data(fh0 + 58);
    const auto *fh0_59 = buffer.data(fh0 + 59);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_23 = buffer.data(fh1 + 23);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_34 = buffer.data(fh1 + 34);
    const auto *fh1_35 = buffer.data(fh1 + 35);
    const auto *fh1_41 = buffer.data(fh1 + 41);
    const auto *fh1_42 = buffer.data(fh1 + 42);
    const auto *fh1_43 = buffer.data(fh1 + 43);
    const auto *fh1_44 = buffer.data(fh1 + 44);
    const auto *fh1_45 = buffer.data(fh1 + 45);
    const auto *fh1_46 = buffer.data(fh1 + 46);
    const auto *fh1_47 = buffer.data(fh1 + 47);
    const auto *fh1_48 = buffer.data(fh1 + 48);
    const auto *fh1_49 = buffer.data(fh1 + 49);
    const auto *fh1_50 = buffer.data(fh1 + 50);
    const auto *fh1_51 = buffer.data(fh1 + 51);
    const auto *fh1_52 = buffer.data(fh1 + 52);
    const auto *fh1_53 = buffer.data(fh1 + 53);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, di_0, fh0_0, fh1_0, fi_0, \
                         fi_1, fi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pb_y[k] * fi_0[k];

        t_2[k] = pb_z[k] * fi_0[k];

        t_3[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_4[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, fh0_1, fh0_2, fh0_3, fh1_1, fh1_2, \
                         fh1_3, fi_3, fi_4, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_6[k] = pb_y[k] * fi_4[k];

        t_7[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_8[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, fh0_4, fh0_5, fh1_4, fh1_5, fi_6, \
                         fi_7, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * fh0_4[k]
                 - f_4 * fh1_4[k]
                 + pb_y[k] * fi_6[k];

        t_10[k] = pb_y[k] * fi_7[k];

        t_11[k] = f_7 * fh0_4[k]
                  - f_8 * fh1_4[k]
                  + pb_z[k] * fi_7[k];

        t_12[k] = f_9 * fh0_5[k]
                  - f_10 * fh1_5[k]
                  + pb_y[k] * fi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, fh0_6, fh0_7, fh1_6, fh1_7, fi_9, \
                         fi_10, fi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fh0_6[k]
                  - f_6 * fh1_6[k]
                  + pb_y[k] * fi_9[k];

        t_14[k] = f_3 * fh0_7[k]
                  - f_4 * fh1_7[k]
                  + pb_y[k] * fi_10[k];

        t_15[k] = pb_y[k] * fi_11[k];

        t_16[k] = f_9 * fh0_7[k]
                  - f_10 * fh1_7[k]
                  + pb_z[k] * fi_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, fh0_8, fh0_10, fh0_11, fh1_8, fh1_9, fh1_10, \
                         fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * fh0_8[k]
                  - f_2 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_18[k] = f_9 * fh0_10[k]
                  - f_10 * fh1_9[k]
                  + pb_y[k] * fi_13[k];

        t_19[k] = f_7 * fh0_11[k]
                  - f_8 * fh1_10[k]
                  + pb_y[k] * fi_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, pb_z, dk_0, fh0_12, fh0_13, \
                         fh1_11, fh1_12, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fh0_12[k]
                  - f_6 * fh1_11[k]
                  + pb_y[k] * fi_15[k];

        t_21[k] = f_3 * fh0_13[k]
                  - f_4 * fh1_12[k]
                  + pb_y[k] * fi_16[k];

        t_22[k] = pb_y[k] * fi_17[k];

        t_23[k] = f_1 * fh0_13[k]
                  - f_2 * fh1_12[k]
                  + pb_z[k] * fi_17[k];

        t_24[k] = pa_y[k] * dk_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_z, pb_x, dk_0, dk_1, dk_2, fh0_25, \
                         fh1_23, fi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * dk_0[k];

        t_26[k] = pa_x[k] * dk_1[k];

        t_27[k] = pa_x[k] * dk_2[k];

        t_28[k] = f_1 * fh0_25[k]
                  - f_2 * fh1_23[k]
                  + pb_x[k] * fi_22[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, fh0_27, fh0_28, fh0_29, fh1_24, fh1_25, \
                         fh1_26, fi_23, fi_24, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * fh0_27[k]
                  - f_10 * fh1_24[k]
                  + pb_x[k] * fi_23[k];

        t_30[k] = f_9 * fh0_28[k]
                  - f_10 * fh1_25[k]
                  + pb_x[k] * fi_24[k];

        t_31[k] = f_7 * fh0_29[k]
                  - f_8 * fh1_26[k]
                  + pb_x[k] * fi_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, fh0_30, fh0_31, fh0_32, fh1_27, fh1_28, \
                         fh1_29, fi_26, fi_27, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * fh0_30[k]
                  - f_8 * fh1_27[k]
                  + pb_x[k] * fi_26[k];

        t_33[k] = f_5 * fh0_31[k]
                  - f_6 * fh1_28[k]
                  + pb_x[k] * fi_27[k];

        t_34[k] = f_5 * fh0_32[k]
                  - f_6 * fh1_29[k]
                  + pb_x[k] * fi_28[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, fh0_33, fh0_34, fh0_36, fh1_30, fh1_31, \
                         fh1_33, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * fh0_33[k]
                  - f_6 * fh1_30[k]
                  + pb_x[k] * fi_29[k];

        t_36[k] = f_3 * fh0_34[k]
                  - f_4 * fh1_31[k]
                  + pb_x[k] * fi_30[k];

        t_37[k] = f_3 * fh0_36[k]
                  - f_4 * fh1_33[k]
                  + pb_x[k] * fi_31[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pb_x, fh0_37, fh0_38, fh1_34, fh1_35, \
                         fi_32, fi_33, fi_34, fi_36, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * fh0_37[k]
                  - f_4 * fh1_34[k]
                  + pb_x[k] * fi_32[k];

        t_39[k] = f_3 * fh0_38[k]
                  - f_4 * fh1_35[k]
                  + pb_x[k] * fi_33[k];

        t_40[k] = pb_x[k] * fi_34[k];

        t_41[k] = pb_x[k] * fi_36[k];

        t_42[k] = pb_x[k] * fi_37[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, di_24, fh0_34, \
                         fh1_31, fi_34, fi_35, fi_38, fi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_x[k] * fi_38[k];

        t_44[k] = pb_x[k] * fi_39[k];

        t_45[k] = f_0 * di_24[k]
                  + f_1 * fh0_34[k]
                  - f_2 * fh1_31[k]
                  + pb_y[k] * fi_34[k];

        t_46[k] = pb_z[k] * fi_34[k];

        t_47[k] = f_3 * fh0_34[k]
                  - f_4 * fh1_31[k]
                  + pb_z[k] * fi_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_z, fh0_35, fh0_36, fh0_37, fh1_32, fh1_33, \
                         fh1_34, fi_36, fi_37, fi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * fh0_35[k]
                  - f_6 * fh1_32[k]
                  + pb_z[k] * fi_36[k];

        t_49[k] = f_7 * fh0_36[k]
                  - f_8 * fh1_33[k]
                  + pb_z[k] * fi_37[k];

        t_50[k] = f_9 * fh0_37[k]
                  - f_10 * fh1_34[k]
                  + pb_z[k] * fi_38[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pa_z, pb_x, pb_z, dk_1, dk_2, fh0_38, \
                         fh0_46, fh1_35, fh1_41, fi_39, fi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * fh0_38[k]
                  - f_2 * fh1_35[k]
                  + pb_z[k] * fi_39[k];

        t_52[k] = pa_z[k] * dk_1[k];

        t_53[k] = pa_y[k] * dk_2[k];

        t_54[k] = f_1 * fh0_46[k]
                  - f_2 * fh1_41[k]
                  + pb_x[k] * fi_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, fh0_48, fh0_49, fh0_50, fh1_42, fh1_43, \
                         fh1_44, fi_43, fi_44, fi_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_9 * fh0_48[k]
                  - f_10 * fh1_42[k]
                  + pb_x[k] * fi_43[k];

        t_56[k] = f_9 * fh0_49[k]
                  - f_10 * fh1_43[k]
                  + pb_x[k] * fi_44[k];

        t_57[k] = f_7 * fh0_50[k]
                  - f_8 * fh1_44[k]
                  + pb_x[k] * fi_45[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, fh0_51, fh0_52, fh0_53, fh1_45, fh1_46, \
                         fh1_47, fi_46, fi_47, fi_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * fh0_51[k]
                  - f_8 * fh1_45[k]
                  + pb_x[k] * fi_46[k];

        t_59[k] = f_5 * fh0_52[k]
                  - f_6 * fh1_46[k]
                  + pb_x[k] * fi_47[k];

        t_60[k] = f_5 * fh0_53[k]
                  - f_6 * fh1_47[k]
                  + pb_x[k] * fi_48[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, fh0_54, fh0_55, fh0_56, fh1_48, fh1_49, \
                         fh1_50, fi_49, fi_50, fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * fh0_54[k]
                  - f_6 * fh1_48[k]
                  + pb_x[k] * fi_49[k];

        t_62[k] = f_3 * fh0_55[k]
                  - f_4 * fh1_49[k]
                  + pb_x[k] * fi_50[k];

        t_63[k] = f_3 * fh0_56[k]
                  - f_4 * fh1_50[k]
                  + pb_x[k] * fi_51[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pb_x, fh0_57, fh0_59, fh1_51, fh1_53, \
                         fi_52, fi_53, fi_54, fi_55, fi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * fh0_57[k]
                  - f_4 * fh1_51[k]
                  + pb_x[k] * fi_52[k];

        t_65[k] = f_3 * fh0_59[k]
                  - f_4 * fh1_53[k]
                  + pb_x[k] * fi_53[k];

        t_66[k] = pb_x[k] * fi_54[k];

        t_67[k] = pb_x[k] * fi_55[k];

        t_68[k] = pb_x[k] * fi_56[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, fh0_55, fh0_56, fh1_49, fh1_50, \
                         fi_54, fi_55, fi_57, fi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_x[k] * fi_57[k];

        t_70[k] = pb_x[k] * fi_59[k];

        t_71[k] = f_1 * fh0_55[k]
                  - f_2 * fh1_49[k]
                  + pb_y[k] * fi_54[k];

        t_72[k] = f_9 * fh0_56[k]
                  - f_10 * fh1_50[k]
                  + pb_y[k] * fi_55[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_y, fh0_57, fh0_58, fh0_59, fh1_51, fh1_52, \
                         fh1_53, fi_56, fi_57, fi_58, fi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * fh0_57[k]
                  - f_8 * fh1_51[k]
                  + pb_y[k] * fi_56[k];

        t_74[k] = f_5 * fh0_58[k]
                  - f_6 * fh1_52[k]
                  + pb_y[k] * fi_57[k];

        t_75[k] = f_3 * fh0_59[k]
                  - f_4 * fh1_53[k]
                  + pb_y[k] * fi_58[k];

        t_76[k] = pb_y[k] * fi_59[k];
    }

#pragma omp simd aligned(t_77, pb_z, di_44, fh0_59, fh1_53, fi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * di_44[k]
                  + f_1 * fh0_59[k]
                  - f_2 * fh1_53[k]
                  + pb_z[k] * fi_59[k];
    }
}

auto
compute_prim_fk_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_8 = buffer.data(fh0 + 8);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_17 = buffer.data(fh1 + 17);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_20 = buffer.data(fi + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_5, fh1_11, \
                         fi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_5[k]
                 - f_2 * fh1_11[k]
                 + pb_y[k] * fi_13[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_8, fh1_17, fi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_8[k]
                 - f_2 * fh1_17[k]
                 + pb_z[k] * fi_20[k];
    }
}

auto
compute_prim_fk_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_17 = buffer.data(fh0 + 17);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_56 = buffer.data(fh1 + 56);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_53 = buffer.data(fi + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_11, fh1_33, \
                         fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_11[k]
                 - f_2 * fh1_33[k]
                 + pb_y[k] * fi_31[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_17, fh1_56, fi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_17[k]
                 - f_2 * fh1_56[k]
                 + pb_z[k] * fi_53[k];
    }
}

auto
compute_prim_fk_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_56 = buffer.data(fh0 + 56);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_50 = buffer.data(fh1 + 50);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_33, fh1_29, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_33[k]
                 - f_2 * fh1_29[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_56, fh1_50, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_56[k]
                 - f_2 * fh1_50[k]
                 + pb_z[k] * fi_8[k];
    }
}

auto
compute_prim_fk_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_17 = buffer.data(fh0 + 17);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_17 = buffer.data(fh1 + 17);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, di_0, di_1, dk_0, dk_1, \
                         fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = f_3 * di_1[k]
                 + pa_x[k] * dk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, di_2, di_3, di_4, di_5, dk_2, \
                         dk_3, dk_4, dk_5, fi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * di_2[k]
                 + pa_x[k] * dk_2[k];

        t_5[k] = f_0 * di_3[k]
                 + pa_x[k] * dk_3[k];

        t_6[k] = f_5 * di_4[k]
                 + pa_x[k] * dk_4[k];

        t_7[k] = f_6 * di_5[k]
                 + pb_x[k] * fi_7[k];

        t_8[k] = pa_x[k] * dk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, di_6, di_7, di_8, di_9, dk_6, dk_7, \
                         dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_6[k]
                 + pa_x[k] * dk_6[k];

        t_10[k] = f_4 * di_7[k]
                  + pa_x[k] * dk_7[k];

        t_11[k] = f_0 * di_8[k]
                  + pa_x[k] * dk_8[k];

        t_12[k] = f_5 * di_9[k]
                  + pa_x[k] * dk_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, pb_x, pb_y, di_5, di_14, dk_5, \
                         dk_14, fh0_11, fh1_11, fi_12, fi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * di_14[k]
                  + pb_x[k] * fi_12[k];

        t_14[k] = pa_x[k] * dk_14[k];

        t_15[k] = f_0 * di_5[k]
                  + f_1 * fh0_11[k]
                  - f_2 * fh1_11[k]
                  + pb_y[k] * fi_13[k];

        t_16[k] = pa_z[k] * dk_5[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, di_10, di_11, di_12, di_13, dk_10, \
                         dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * di_10[k]
                  + pa_y[k] * dk_10[k];

        t_18[k] = f_4 * di_11[k]
                  + pa_y[k] * dk_11[k];

        t_19[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_20[k] = f_5 * di_13[k]
                  + pa_y[k] * dk_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_y, pb_z, di_14, dk_14, fh0_17, fh1_17, \
                         fi_19, fi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_6 * di_14[k]
                  + pb_y[k] * fi_19[k];

        t_22[k] = pa_y[k] * dk_14[k];

        t_23[k] = f_0 * di_14[k]
                  + f_1 * fh0_17[k]
                  - f_2 * fh1_17[k]
                  + pb_z[k] * fi_20[k];
    }
}

auto
compute_prim_fk_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_17 = buffer.data(fh0 + 17);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_53 = buffer.data(fh1 + 53);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_68 = buffer.data(fi + 68);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, di_0, di_1, dk_0, dk_1, \
                         fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = f_3 * di_1[k]
                 + pa_x[k] * dk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, di_2, di_3, di_4, di_5, dk_2, \
                         dk_3, dk_4, dk_5, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * di_2[k]
                 + pa_x[k] * dk_2[k];

        t_5[k] = f_0 * di_3[k]
                 + pa_x[k] * dk_3[k];

        t_6[k] = f_5 * di_4[k]
                 + pa_x[k] * dk_4[k];

        t_7[k] = f_6 * di_5[k]
                 + pb_x[k] * fi_24[k];

        t_8[k] = pa_x[k] * dk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, di_6, di_7, di_8, di_9, dk_6, dk_7, \
                         dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_6[k]
                 + pa_x[k] * dk_6[k];

        t_10[k] = f_4 * di_7[k]
                  + pa_x[k] * dk_7[k];

        t_11[k] = f_0 * di_8[k]
                  + pa_x[k] * dk_8[k];

        t_12[k] = f_5 * di_9[k]
                  + pa_x[k] * dk_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, pb_x, pb_y, di_5, di_14, dk_5, \
                         dk_14, fh0_11, fh1_31, fi_30, fi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * di_14[k]
                  + pb_x[k] * fi_30[k];

        t_14[k] = pa_x[k] * dk_14[k];

        t_15[k] = f_0 * di_5[k]
                  + f_1 * fh0_11[k]
                  - f_2 * fh1_31[k]
                  + pb_y[k] * fi_40[k];

        t_16[k] = pa_z[k] * dk_5[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, di_10, di_11, di_12, di_13, dk_10, \
                         dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * di_10[k]
                  + pa_y[k] * dk_10[k];

        t_18[k] = f_4 * di_11[k]
                  + pa_y[k] * dk_11[k];

        t_19[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_20[k] = f_5 * di_13[k]
                  + pa_y[k] * dk_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_y, pb_z, di_14, dk_14, fh0_17, fh1_53, \
                         fi_53, fi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_6 * di_14[k]
                  + pb_y[k] * fi_53[k];

        t_22[k] = pa_y[k] * dk_14[k];

        t_23[k] = f_0 * di_14[k]
                  + f_1 * fh0_17[k]
                  - f_2 * fh1_53[k]
                  + pb_z[k] * fi_68[k];
    }
}

auto
compute_prim_fk_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 1.0 / p;
    const auto f_14 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_4 = buffer.data(fh0 + 4);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_7 = buffer.data(fh0 + 7);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_10 = buffer.data(fh0 + 10);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_28 = buffer.data(fh0 + 28);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);
    const auto *fh0_33 = buffer.data(fh0 + 33);
    const auto *fh0_34 = buffer.data(fh0 + 34);
    const auto *fh0_35 = buffer.data(fh0 + 35);
    const auto *fh0_41 = buffer.data(fh0 + 41);
    const auto *fh0_43 = buffer.data(fh0 + 43);
    const auto *fh0_44 = buffer.data(fh0 + 44);
    const auto *fh0_45 = buffer.data(fh0 + 45);
    const auto *fh0_46 = buffer.data(fh0 + 46);
    const auto *fh0_47 = buffer.data(fh0 + 47);
    const auto *fh0_48 = buffer.data(fh0 + 48);
    const auto *fh0_49 = buffer.data(fh0 + 49);
    const auto *fh0_50 = buffer.data(fh0 + 50);
    const auto *fh0_51 = buffer.data(fh0 + 51);
    const auto *fh0_52 = buffer.data(fh0 + 52);
    const auto *fh0_53 = buffer.data(fh0 + 53);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_4 = buffer.data(fh1 + 4);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_7 = buffer.data(fh1 + 7);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_10 = buffer.data(fh1 + 10);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_22 = buffer.data(fh1 + 22);
    const auto *fh1_23 = buffer.data(fh1 + 23);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_28 = buffer.data(fh1 + 28);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);
    const auto *fh1_33 = buffer.data(fh1 + 33);
    const auto *fh1_39 = buffer.data(fh1 + 39);
    const auto *fh1_40 = buffer.data(fh1 + 40);
    const auto *fh1_41 = buffer.data(fh1 + 41);
    const auto *fh1_42 = buffer.data(fh1 + 42);
    const auto *fh1_43 = buffer.data(fh1 + 43);
    const auto *fh1_44 = buffer.data(fh1 + 44);
    const auto *fh1_45 = buffer.data(fh1 + 45);
    const auto *fh1_46 = buffer.data(fh1 + 46);
    const auto *fh1_47 = buffer.data(fh1 + 47);
    const auto *fh1_48 = buffer.data(fh1 + 48);
    const auto *fh1_49 = buffer.data(fh1 + 49);
    const auto *fh1_50 = buffer.data(fh1 + 50);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, di_0, fh0_0, fh0_1, fh1_0, \
                         fh1_1, fi_0, fi_1, fi_2, fi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_y[k] * fi_1[k];

        t_2[k] = f_3 * fh0_0[k]
                 - f_4 * fh1_0[k]
                 + pb_z[k] * fi_2[k];

        t_3[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, fh0_2, fh0_3, fh0_4, fh1_2, fh1_3, fh1_4, \
                         fi_4, fi_5, fi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_4[k];

        t_5[k] = f_7 * fh0_3[k]
                 - f_8 * fh1_3[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = f_7 * fh0_4[k]
                 - f_8 * fh1_4[k]
                 + pb_z[k] * fi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, fh0_5, fh0_6, fh0_7, fh1_5, fh1_6, fh1_7, \
                         fi_7, fi_8, fi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_9 * fh0_5[k]
                 - f_10 * fh1_5[k]
                 + pb_y[k] * fi_7[k];

        t_8[k] = f_9 * fh0_6[k]
                 - f_10 * fh1_6[k]
                 + pb_z[k] * fi_8[k];

        t_9[k] = f_1 * fh0_7[k]
                 - f_2 * fh1_7[k]
                 + pb_y[k] * fi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, fh0_9, fh0_10, fh0_11, fh1_8, fh1_9, fh1_10, \
                         fi_10, fi_11, fi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_8[k]
                  + pb_y[k] * fi_10[k];

        t_11[k] = f_7 * fh0_10[k]
                  - f_8 * fh1_9[k]
                  + pb_y[k] * fi_11[k];

        t_12[k] = f_5 * fh0_11[k]
                  - f_6 * fh1_10[k]
                  + pb_y[k] * fi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, pb_z, dk_0, fh0_12, fh1_11, \
                         fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * fh0_12[k]
                  - f_4 * fh1_11[k]
                  + pb_y[k] * fi_13[k];

        t_14[k] = f_1 * fh0_12[k]
                  - f_2 * fh1_11[k]
                  + pb_z[k] * fi_14[k];

        t_15[k] = pa_y[k] * dk_0[k];

        t_16[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, di_4, di_5, di_6, di_7, dk_1, dk_2, \
                         dk_3, dk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_11 * di_4[k]
                  + pa_x[k] * dk_1[k];

        t_18[k] = f_12 * di_5[k]
                  + pa_x[k] * dk_2[k];

        t_19[k] = f_0 * di_6[k]
                  + pa_x[k] * dk_3[k];

        t_20[k] = f_13 * di_7[k]
                  + pa_x[k] * dk_4[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pb_x, di_8, di_11, di_12, di_13, \
                         dk_5, dk_6, dk_7, dk_8, fi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_14 * di_8[k]
                  + pb_x[k] * fi_21[k];

        t_22[k] = pa_x[k] * dk_5[k];

        t_23[k] = f_11 * di_11[k]
                  + pa_x[k] * dk_6[k];

        t_24[k] = f_12 * di_12[k]
                  + pa_x[k] * dk_7[k];

        t_25[k] = f_0 * di_13[k]
                  + pa_x[k] * dk_8[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pb_x, di_14, di_20, dk_9, dk_14, \
                         fh0_23, fh1_22, fi_26, fi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_13 * di_14[k]
                  + pa_x[k] * dk_9[k];

        t_27[k] = f_14 * di_20[k]
                  + pb_x[k] * fi_26[k];

        t_28[k] = pa_x[k] * dk_14[k];

        t_29[k] = f_1 * fh0_23[k]
                  - f_2 * fh1_22[k]
                  + pb_x[k] * fi_27[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, fh0_25, fh0_26, fh0_27, fh1_23, fh1_24, \
                         fh1_25, fi_28, fi_29, fi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * fh0_25[k]
                  - f_10 * fh1_23[k]
                  + pb_x[k] * fi_28[k];

        t_31[k] = f_9 * fh0_26[k]
                  - f_10 * fh1_24[k]
                  + pb_x[k] * fi_29[k];

        t_32[k] = f_7 * fh0_27[k]
                  - f_8 * fh1_25[k]
                  + pb_x[k] * fi_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, fh0_28, fh0_29, fh0_30, fh1_26, fh1_27, \
                         fh1_28, fi_31, fi_32, fi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * fh0_28[k]
                  - f_8 * fh1_26[k]
                  + pb_x[k] * fi_31[k];

        t_34[k] = f_5 * fh0_29[k]
                  - f_6 * fh1_27[k]
                  + pb_x[k] * fi_32[k];

        t_35[k] = f_5 * fh0_30[k]
                  - f_6 * fh1_28[k]
                  + pb_x[k] * fi_33[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, pb_y, pb_z, di_8, fh0_31, fh0_35, \
                         fh1_29, fh1_33, fi_34, fi_35, fi_36, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * fh0_31[k]
                  - f_4 * fh1_29[k]
                  + pb_x[k] * fi_34[k];

        t_37[k] = f_3 * fh0_35[k]
                  - f_4 * fh1_33[k]
                  + pb_x[k] * fi_35[k];

        t_38[k] = f_0 * di_8[k]
                  + f_1 * fh0_31[k]
                  - f_2 * fh1_29[k]
                  + pb_y[k] * fi_36[k];

        t_39[k] = f_3 * fh0_31[k]
                  - f_4 * fh1_29[k]
                  + pb_z[k] * fi_37[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_z, fh0_32, fh0_33, fh0_34, fh1_30, fh1_31, \
                         fh1_32, fi_38, fi_39, fi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * fh0_32[k]
                  - f_6 * fh1_30[k]
                  + pb_z[k] * fi_38[k];

        t_41[k] = f_7 * fh0_33[k]
                  - f_8 * fh1_31[k]
                  + pb_z[k] * fi_39[k];

        t_42[k] = f_9 * fh0_34[k]
                  - f_10 * fh1_32[k]
                  + pb_z[k] * fi_40[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pa_z, pb_z, di_16, di_17, dk_5, dk_10, \
                         dk_11, fh0_35, fh1_33, fi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * fh0_35[k]
                  - f_2 * fh1_33[k]
                  + pb_z[k] * fi_41[k];

        t_44[k] = pa_z[k] * dk_5[k];

        t_45[k] = f_11 * di_16[k]
                  + pa_y[k] * dk_10[k];

        t_46[k] = f_12 * di_17[k]
                  + pa_y[k] * dk_11[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_y, di_18, di_19, di_20, dk_12, \
                         dk_13, dk_14, fi_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * di_18[k]
                  + pa_y[k] * dk_12[k];

        t_48[k] = f_13 * di_19[k]
                  + pa_y[k] * dk_13[k];

        t_49[k] = f_14 * di_20[k]
                  + pb_y[k] * fi_47[k];

        t_50[k] = pa_y[k] * dk_14[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, fh0_41, fh0_43, fh0_44, fh1_39, fh1_40, \
                         fh1_41, fi_48, fi_49, fi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * fh0_41[k]
                  - f_2 * fh1_39[k]
                  + pb_x[k] * fi_48[k];

        t_52[k] = f_9 * fh0_43[k]
                  - f_10 * fh1_40[k]
                  + pb_x[k] * fi_49[k];

        t_53[k] = f_9 * fh0_44[k]
                  - f_10 * fh1_41[k]
                  + pb_x[k] * fi_50[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, fh0_45, fh0_46, fh0_47, fh1_42, fh1_43, \
                         fh1_44, fi_51, fi_52, fi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * fh0_45[k]
                  - f_8 * fh1_42[k]
                  + pb_x[k] * fi_51[k];

        t_55[k] = f_7 * fh0_46[k]
                  - f_8 * fh1_43[k]
                  + pb_x[k] * fi_52[k];

        t_56[k] = f_5 * fh0_47[k]
                  - f_6 * fh1_44[k]
                  + pb_x[k] * fi_53[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_x, pb_y, fh0_48, fh0_49, fh0_53, fh1_45, \
                         fh1_46, fh1_50, fi_54, fi_55, fi_56, fi_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * fh0_48[k]
                  - f_6 * fh1_45[k]
                  + pb_x[k] * fi_54[k];

        t_58[k] = f_3 * fh0_49[k]
                  - f_4 * fh1_46[k]
                  + pb_x[k] * fi_55[k];

        t_59[k] = f_3 * fh0_53[k]
                  - f_4 * fh1_50[k]
                  + pb_x[k] * fi_56[k];

        t_60[k] = f_1 * fh0_49[k]
                  - f_2 * fh1_46[k]
                  + pb_y[k] * fi_57[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_y, fh0_50, fh0_51, fh0_52, fh1_47, fh1_48, \
                         fh1_49, fi_58, fi_59, fi_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * fh0_50[k]
                  - f_10 * fh1_47[k]
                  + pb_y[k] * fi_58[k];

        t_62[k] = f_7 * fh0_51[k]
                  - f_8 * fh1_48[k]
                  + pb_y[k] * fi_59[k];

        t_63[k] = f_5 * fh0_52[k]
                  - f_6 * fh1_49[k]
                  + pb_y[k] * fi_60[k];
    }

#pragma omp simd aligned(t_64, t_65, pb_y, pb_z, di_20, fh0_53, fh1_50, fi_61, \
                         fi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * fh0_53[k]
                  - f_4 * fh1_50[k]
                  + pb_y[k] * fi_61[k];

        t_65[k] = f_0 * di_20[k]
                  + f_1 * fh0_53[k]
                  - f_2 * fh1_50[k]
                  + pb_z[k] * fi_62[k];
    }
}

auto
compute_prim_fk_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_14 = buffer.data(di + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_50 = buffer.data(fh0 + 50);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_50 = buffer.data(fh1 + 50);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_5, dk_1, dk_2, fh0_29, fh1_29, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_5[k]
                 + f_1 * fh0_29[k]
                 - f_2 * fh1_29[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_14, fh0_50, fh1_50, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_14[k]
                 + f_1 * fh0_50[k]
                 - f_2 * fh1_50[k]
                 + pb_z[k] * fi_8[k];
    }
}

auto
compute_prim_fk_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_8 = buffer.data(fh0 + 8);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_29 = buffer.data(fi + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_1, dk_1, dk_2, fh0_5, fh1_14, \
                         fi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_1[k]
                 + f_1 * fh0_5[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_16[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_2, fh0_8, fh1_26, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_2[k]
                 + f_1 * fh0_8[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_29[k];
    }
}

auto
compute_prim_fk_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_26 = buffer.data(fh0 + 26);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_29 = buffer.data(fi + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_14, fh1_14, \
                         fi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_14[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_16[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_26, fh1_26, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_26[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_29[k];
    }
}

auto
compute_prim_fk_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_26 = buffer.data(fh0 + 26);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_17 = buffer.data(fh1 + 17);
    const auto *fh1_32 = buffer.data(fh1 + 32);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_29 = buffer.data(fi + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_14, fh1_17, \
                         fi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_14[k]
                 - f_2 * fh1_17[k]
                 + pb_y[k] * fi_16[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_26, fh1_32, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_26[k]
                 - f_2 * fh1_32[k]
                 + pb_z[k] * fi_29[k];
    }
}

auto
compute_prim_fk_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_32 = buffer.data(fh0 + 32);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_17, fh1_14, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_17[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_32, fh1_26, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_32[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_8[k];
    }
}

auto
compute_prim_fk_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_8 = buffer.data(fh0 + 8);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_32 = buffer.data(fi + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_1, dk_1, dk_2, fh0_5, fh1_14, \
                         fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_1[k]
                 + f_1 * fh0_5[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_17[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_2, fh0_8, fh1_26, fi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_2[k]
                 + f_1 * fh0_8[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_32[k];
    }
}

auto
compute_prim_fk_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 2.0 / beta;
    const auto f_8 = 2.0 * alpha / (beta * p);
    const auto f_9 = 1.5 / beta;
    const auto f_10 = 1.5 * alpha / (beta * p);
    const auto f_11 = 1.0 / beta;
    const auto f_12 = alpha / (beta * p);
    const auto f_13 = 0.5 / beta;
    const auto f_14 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_13 = buffer.data(fh0 + 13);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_20 = buffer.data(fh0 + 20);
    const auto *fh0_21 = buffer.data(fh0 + 21);
    const auto *fh0_22 = buffer.data(fh0 + 22);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_11 = buffer.data(fh1 + 11);
    const auto *fh1_12 = buffer.data(fh1 + 12);
    const auto *fh1_13 = buffer.data(fh1 + 13);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_20 = buffer.data(fh1 + 20);
    const auto *fh1_21 = buffer.data(fh1 + 21);
    const auto *fh1_22 = buffer.data(fh1 + 22);
    const auto *fh1_23 = buffer.data(fh1 + 23);
    const auto *fh1_24 = buffer.data(fh1 + 24);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, di_0, di_1, dk_0, dk_1, \
                         fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = f_3 * di_1[k]
                 + pa_x[k] * dk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, di_2, di_3, di_4, di_5, dk_2, \
                         dk_3, dk_4, dk_5, fi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * di_2[k]
                 + pa_x[k] * dk_2[k];

        t_5[k] = f_0 * di_3[k]
                 + pa_x[k] * dk_3[k];

        t_6[k] = f_5 * di_4[k]
                 + pa_x[k] * dk_4[k];

        t_7[k] = f_6 * di_5[k]
                 + pb_x[k] * fi_7[k];

        t_8[k] = pa_x[k] * dk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, di_6, di_7, di_8, di_9, dk_6, dk_7, \
                         dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_6[k]
                 + pa_x[k] * dk_6[k];

        t_10[k] = f_4 * di_7[k]
                  + pa_x[k] * dk_7[k];

        t_11[k] = f_0 * di_8[k]
                  + pa_x[k] * dk_8[k];

        t_12[k] = f_5 * di_9[k]
                  + pa_x[k] * dk_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, di_14, dk_14, fh0_11, fh0_12, \
                         fh1_11, fh1_12, fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * di_14[k]
                  + pb_x[k] * fi_12[k];

        t_14[k] = pa_x[k] * dk_14[k];

        t_15[k] = f_7 * fh0_11[k]
                  - f_8 * fh1_11[k]
                  + pb_x[k] * fi_13[k];

        t_16[k] = f_9 * fh0_12[k]
                  - f_10 * fh1_12[k]
                  + pb_x[k] * fi_14[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, di_5, dk_5, fh0_13, fh0_14, \
                         fh1_13, fh1_14, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_11 * fh0_13[k]
                  - f_12 * fh1_13[k]
                  + pb_x[k] * fi_15[k];

        t_18[k] = f_13 * fh0_14[k]
                  - f_14 * fh1_14[k]
                  + pb_x[k] * fi_16[k];

        t_19[k] = f_0 * di_5[k]
                  + f_1 * fh0_14[k]
                  - f_2 * fh1_14[k]
                  + pb_y[k] * fi_17[k];

        t_20[k] = pa_z[k] * dk_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_y, di_10, di_11, di_12, di_13, dk_10, \
                         dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * di_10[k]
                  + pa_y[k] * dk_10[k];

        t_22[k] = f_4 * di_11[k]
                  + pa_y[k] * dk_11[k];

        t_23[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_24[k] = f_5 * di_13[k]
                  + pa_y[k] * dk_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_x, pb_y, di_14, dk_14, fh0_20, \
                         fh0_21, fh1_20, fh1_21, fi_23, fi_24, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * di_14[k]
                  + pb_y[k] * fi_23[k];

        t_26[k] = pa_y[k] * dk_14[k];

        t_27[k] = f_7 * fh0_20[k]
                  - f_8 * fh1_20[k]
                  + pb_x[k] * fi_24[k];

        t_28[k] = f_9 * fh0_21[k]
                  - f_10 * fh1_21[k]
                  + pb_x[k] * fi_25[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, fh0_22, fh0_23, fh0_26, fh1_22, fh1_23, \
                         fh1_26, fi_26, fi_27, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * fh0_22[k]
                  - f_12 * fh1_22[k]
                  + pb_x[k] * fi_26[k];

        t_30[k] = f_13 * fh0_26[k]
                  - f_14 * fh1_26[k]
                  + pb_x[k] * fi_27[k];

        t_31[k] = f_7 * fh0_23[k]
                  - f_8 * fh1_23[k]
                  + pb_y[k] * fi_28[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_y, fh0_24, fh0_25, fh0_26, fh1_24, fh1_25, \
                         fh1_26, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * fh0_24[k]
                  - f_10 * fh1_24[k]
                  + pb_y[k] * fi_29[k];

        t_33[k] = f_11 * fh0_25[k]
                  - f_12 * fh1_25[k]
                  + pb_y[k] * fi_30[k];

        t_34[k] = f_13 * fh0_26[k]
                  - f_14 * fh1_26[k]
                  + pb_y[k] * fi_31[k];
    }

#pragma omp simd aligned(t_35, pb_z, di_14, fh0_26, fh1_26, fi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * di_14[k]
                  + f_1 * fh0_26[k]
                  - f_2 * fh1_26[k]
                  + pb_z[k] * fi_32[k];
    }
}

auto
compute_prim_fk_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 2.0 / beta;
    const auto f_8 = 2.0 * alpha / (beta * p);
    const auto f_9 = 1.5 / beta;
    const auto f_10 = 1.5 * alpha / (beta * p);
    const auto f_11 = 1.0 / beta;
    const auto f_12 = alpha / (beta * p);
    const auto f_13 = 0.5 / beta;
    const auto f_14 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_12 = buffer.data(fh0 + 12);
    const auto *fh0_13 = buffer.data(fh0 + 13);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_20 = buffer.data(fh0 + 20);
    const auto *fh0_21 = buffer.data(fh0 + 21);
    const auto *fh0_22 = buffer.data(fh0 + 22);
    const auto *fh0_23 = buffer.data(fh0 + 23);
    const auto *fh0_24 = buffer.data(fh0 + 24);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_15 = buffer.data(fh1 + 15);
    const auto *fh1_16 = buffer.data(fh1 + 16);
    const auto *fh1_17 = buffer.data(fh1 + 17);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, di_0, di_1, dk_0, dk_1, \
                         fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = f_3 * di_1[k]
                 + pa_x[k] * dk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, di_2, di_3, di_4, di_5, dk_2, \
                         dk_3, dk_4, dk_5, fi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * di_2[k]
                 + pa_x[k] * dk_2[k];

        t_5[k] = f_0 * di_3[k]
                 + pa_x[k] * dk_3[k];

        t_6[k] = f_5 * di_4[k]
                 + pa_x[k] * dk_4[k];

        t_7[k] = f_6 * di_5[k]
                 + pb_x[k] * fi_9[k];

        t_8[k] = pa_x[k] * dk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, di_6, di_7, di_8, di_9, dk_6, dk_7, \
                         dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_6[k]
                 + pa_x[k] * dk_6[k];

        t_10[k] = f_4 * di_7[k]
                  + pa_x[k] * dk_7[k];

        t_11[k] = f_0 * di_8[k]
                  + pa_x[k] * dk_8[k];

        t_12[k] = f_5 * di_9[k]
                  + pa_x[k] * dk_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, di_14, dk_14, fh0_11, fh0_12, \
                         fh1_14, fh1_15, fi_14, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * di_14[k]
                  + pb_x[k] * fi_14[k];

        t_14[k] = pa_x[k] * dk_14[k];

        t_15[k] = f_7 * fh0_11[k]
                  - f_8 * fh1_14[k]
                  + pb_x[k] * fi_16[k];

        t_16[k] = f_9 * fh0_12[k]
                  - f_10 * fh1_15[k]
                  + pb_x[k] * fi_17[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, di_5, dk_5, fh0_13, fh0_14, \
                         fh1_16, fh1_17, fi_18, fi_19, fi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_11 * fh0_13[k]
                  - f_12 * fh1_16[k]
                  + pb_x[k] * fi_18[k];

        t_18[k] = f_13 * fh0_14[k]
                  - f_14 * fh1_17[k]
                  + pb_x[k] * fi_19[k];

        t_19[k] = f_0 * di_5[k]
                  + f_1 * fh0_14[k]
                  - f_2 * fh1_17[k]
                  + pb_y[k] * fi_20[k];

        t_20[k] = pa_z[k] * dk_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_y, di_10, di_11, di_12, di_13, dk_10, \
                         dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * di_10[k]
                  + pa_y[k] * dk_10[k];

        t_22[k] = f_4 * di_11[k]
                  + pa_y[k] * dk_11[k];

        t_23[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_24[k] = f_5 * di_13[k]
                  + pa_y[k] * dk_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_x, pb_y, di_14, dk_14, fh0_20, \
                         fh0_21, fh1_25, fh1_26, fi_27, fi_29, fi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * di_14[k]
                  + pb_y[k] * fi_27[k];

        t_26[k] = pa_y[k] * dk_14[k];

        t_27[k] = f_7 * fh0_20[k]
                  - f_8 * fh1_25[k]
                  + pb_x[k] * fi_29[k];

        t_28[k] = f_9 * fh0_21[k]
                  - f_10 * fh1_26[k]
                  + pb_x[k] * fi_30[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, fh0_22, fh0_23, fh0_26, fh1_27, fh1_29, \
                         fh1_32, fi_31, fi_32, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * fh0_22[k]
                  - f_12 * fh1_27[k]
                  + pb_x[k] * fi_31[k];

        t_30[k] = f_13 * fh0_26[k]
                  - f_14 * fh1_32[k]
                  + pb_x[k] * fi_32[k];

        t_31[k] = f_7 * fh0_23[k]
                  - f_8 * fh1_29[k]
                  + pb_y[k] * fi_34[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_y, fh0_24, fh0_25, fh0_26, fh1_30, fh1_31, \
                         fh1_32, fi_35, fi_36, fi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * fh0_24[k]
                  - f_10 * fh1_30[k]
                  + pb_y[k] * fi_35[k];

        t_33[k] = f_11 * fh0_25[k]
                  - f_12 * fh1_31[k]
                  + pb_y[k] * fi_36[k];

        t_34[k] = f_13 * fh0_26[k]
                  - f_14 * fh1_32[k]
                  + pb_y[k] * fi_37[k];
    }

#pragma omp simd aligned(t_35, pb_z, di_14, fh0_26, fh1_32, fi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * di_14[k]
                  + f_1 * fh0_26[k]
                  - f_2 * fh1_32[k]
                  + pb_z[k] * fi_38[k];
    }
}

auto
compute_prim_fk_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 2.0 / beta;
    const auto f_8 = 2.0 * alpha / (beta * p);
    const auto f_9 = 1.5 / beta;
    const auto f_10 = 1.5 * alpha / (beta * p);
    const auto f_11 = 1.0 / beta;
    const auto f_12 = alpha / (beta * p);
    const auto f_13 = 0.5 / beta;
    const auto f_14 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_15 = buffer.data(fh0 + 15);
    const auto *fh0_16 = buffer.data(fh0 + 16);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_25 = buffer.data(fh0 + 25);
    const auto *fh0_26 = buffer.data(fh0 + 26);
    const auto *fh0_27 = buffer.data(fh0 + 27);
    const auto *fh0_29 = buffer.data(fh0 + 29);
    const auto *fh0_30 = buffer.data(fh0 + 30);
    const auto *fh0_31 = buffer.data(fh0 + 31);
    const auto *fh0_32 = buffer.data(fh0 + 32);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_15 = buffer.data(fh1 + 15);
    const auto *fh1_16 = buffer.data(fh1 + 16);
    const auto *fh1_17 = buffer.data(fh1 + 17);
    const auto *fh1_25 = buffer.data(fh1 + 25);
    const auto *fh1_26 = buffer.data(fh1 + 26);
    const auto *fh1_27 = buffer.data(fh1 + 27);
    const auto *fh1_29 = buffer.data(fh1 + 29);
    const auto *fh1_30 = buffer.data(fh1 + 30);
    const auto *fh1_31 = buffer.data(fh1 + 31);
    const auto *fh1_32 = buffer.data(fh1 + 32);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, di_0, di_1, dk_0, dk_1, \
                         fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = f_3 * di_1[k]
                 + pa_x[k] * dk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, di_2, di_3, di_4, di_5, dk_2, \
                         dk_3, dk_4, dk_5, fi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * di_2[k]
                 + pa_x[k] * dk_2[k];

        t_5[k] = f_0 * di_3[k]
                 + pa_x[k] * dk_3[k];

        t_6[k] = f_5 * di_4[k]
                 + pa_x[k] * dk_4[k];

        t_7[k] = f_6 * di_5[k]
                 + pb_x[k] * fi_7[k];

        t_8[k] = pa_x[k] * dk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, di_6, di_7, di_8, di_9, dk_6, dk_7, \
                         dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_6[k]
                 + pa_x[k] * dk_6[k];

        t_10[k] = f_4 * di_7[k]
                  + pa_x[k] * dk_7[k];

        t_11[k] = f_0 * di_8[k]
                  + pa_x[k] * dk_8[k];

        t_12[k] = f_5 * di_9[k]
                  + pa_x[k] * dk_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, di_14, dk_14, fh0_14, fh0_15, \
                         fh1_14, fh1_15, fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * di_14[k]
                  + pb_x[k] * fi_12[k];

        t_14[k] = pa_x[k] * dk_14[k];

        t_15[k] = f_7 * fh0_14[k]
                  - f_8 * fh1_14[k]
                  + pb_x[k] * fi_13[k];

        t_16[k] = f_9 * fh0_15[k]
                  - f_10 * fh1_15[k]
                  + pb_x[k] * fi_14[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, di_5, dk_5, fh0_16, fh0_17, \
                         fh1_16, fh1_17, fi_15, fi_16, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_11 * fh0_16[k]
                  - f_12 * fh1_16[k]
                  + pb_x[k] * fi_15[k];

        t_18[k] = f_13 * fh0_17[k]
                  - f_14 * fh1_17[k]
                  + pb_x[k] * fi_16[k];

        t_19[k] = f_0 * di_5[k]
                  + f_1 * fh0_17[k]
                  - f_2 * fh1_17[k]
                  + pb_y[k] * fi_17[k];

        t_20[k] = pa_z[k] * dk_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_y, di_10, di_11, di_12, di_13, dk_10, \
                         dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * di_10[k]
                  + pa_y[k] * dk_10[k];

        t_22[k] = f_4 * di_11[k]
                  + pa_y[k] * dk_11[k];

        t_23[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_24[k] = f_5 * di_13[k]
                  + pa_y[k] * dk_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_x, pb_y, di_14, dk_14, fh0_25, \
                         fh0_26, fh1_25, fh1_26, fi_23, fi_24, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * di_14[k]
                  + pb_y[k] * fi_23[k];

        t_26[k] = pa_y[k] * dk_14[k];

        t_27[k] = f_7 * fh0_25[k]
                  - f_8 * fh1_25[k]
                  + pb_x[k] * fi_24[k];

        t_28[k] = f_9 * fh0_26[k]
                  - f_10 * fh1_26[k]
                  + pb_x[k] * fi_25[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, fh0_27, fh0_29, fh0_32, fh1_27, fh1_29, \
                         fh1_32, fi_26, fi_27, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * fh0_27[k]
                  - f_12 * fh1_27[k]
                  + pb_x[k] * fi_26[k];

        t_30[k] = f_13 * fh0_32[k]
                  - f_14 * fh1_32[k]
                  + pb_x[k] * fi_27[k];

        t_31[k] = f_7 * fh0_29[k]
                  - f_8 * fh1_29[k]
                  + pb_y[k] * fi_28[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_y, fh0_30, fh0_31, fh0_32, fh1_30, fh1_31, \
                         fh1_32, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * fh0_30[k]
                  - f_10 * fh1_30[k]
                  + pb_y[k] * fi_29[k];

        t_33[k] = f_11 * fh0_31[k]
                  - f_12 * fh1_31[k]
                  + pb_y[k] * fi_30[k];

        t_34[k] = f_13 * fh0_32[k]
                  - f_14 * fh1_32[k]
                  + pb_y[k] * fi_31[k];
    }

#pragma omp simd aligned(t_35, pb_z, di_14, fh0_32, fh1_32, fi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * di_14[k]
                  + f_1 * fh0_32[k]
                  - f_2 * fh1_32[k]
                  + pb_z[k] * fi_32[k];
    }
}

auto
compute_prim_fk_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_14 = buffer.data(di + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_32 = buffer.data(fh0 + 32);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_5, dk_1, dk_2, fh0_17, fh1_14, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_5[k]
                 + f_1 * fh0_17[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_14, fh0_32, fh1_26, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_14[k]
                 + f_1 * fh0_32[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_8[k];
    }
}

auto
compute_prim_fk_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_17 = buffer.data(fh0 + 17);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_29 = buffer.data(fi + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_11, fh1_14, \
                         fi_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_11[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_16[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_17, fh1_26, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_17[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_29[k];
    }
}

auto
compute_prim_fk_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_26 = buffer.data(fh0 + 26);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_4, dk_1, dk_2, fh0_14, fh1_14, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_4[k]
                 + f_1 * fh0_14[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_11, fh0_26, fh1_26, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_11[k]
                 + f_1 * fh0_26[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_8[k];
    }
}

auto
compute_prim_fk_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);

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

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_11 = buffer.data(fh0 + 11);
    const auto *fh0_17 = buffer.data(fh0 + 17);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_32 = buffer.data(fi + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, di_0, di_1, dk_0, dk_1, \
                         fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = f_3 * di_1[k]
                 + pa_x[k] * dk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, di_2, di_3, di_4, di_5, dk_2, \
                         dk_3, dk_4, dk_5, fi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * di_2[k]
                 + pa_x[k] * dk_2[k];

        t_5[k] = f_0 * di_3[k]
                 + pa_x[k] * dk_3[k];

        t_6[k] = f_5 * di_4[k]
                 + pa_x[k] * dk_4[k];

        t_7[k] = f_6 * di_5[k]
                 + pb_x[k] * fi_7[k];

        t_8[k] = pa_x[k] * dk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, di_6, di_7, di_8, di_9, dk_6, dk_7, \
                         dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_6[k]
                 + pa_x[k] * dk_6[k];

        t_10[k] = f_4 * di_7[k]
                  + pa_x[k] * dk_7[k];

        t_11[k] = f_0 * di_8[k]
                  + pa_x[k] * dk_8[k];

        t_12[k] = f_5 * di_9[k]
                  + pa_x[k] * dk_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, pb_x, pb_y, di_5, di_14, dk_5, \
                         dk_14, fh0_11, fh1_14, fi_12, fi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * di_14[k]
                  + pb_x[k] * fi_12[k];

        t_14[k] = pa_x[k] * dk_14[k];

        t_15[k] = f_0 * di_5[k]
                  + f_1 * fh0_11[k]
                  - f_2 * fh1_14[k]
                  + pb_y[k] * fi_17[k];

        t_16[k] = pa_z[k] * dk_5[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, di_10, di_11, di_12, di_13, dk_10, \
                         dk_11, dk_12, dk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * di_10[k]
                  + pa_y[k] * dk_10[k];

        t_18[k] = f_4 * di_11[k]
                  + pa_y[k] * dk_11[k];

        t_19[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_12[k];

        t_20[k] = f_5 * di_13[k]
                  + pa_y[k] * dk_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_y, pb_z, di_14, dk_14, fh0_17, fh1_26, \
                         fi_23, fi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_6 * di_14[k]
                  + pb_y[k] * fi_23[k];

        t_22[k] = pa_y[k] * dk_14[k];

        t_23[k] = f_0 * di_14[k]
                  + f_1 * fh0_17[k]
                  - f_2 * fh1_26[k]
                  + pb_z[k] * fi_32[k];
    }
}

auto
compute_prim_fk_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di, const size_t dk,
                                      const size_t fh0, const size_t fh1, const size_t fi,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_14 = buffer.data(di + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_14 = buffer.data(fh0 + 14);
    const auto *fh0_26 = buffer.data(fh0 + 26);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_14 = buffer.data(fh1 + 14);
    const auto *fh1_26 = buffer.data(fh1 + 26);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_8 = buffer.data(fi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, di_0, dk_0, dk_1, \
                         dk_2, fh0_0, fh1_0, fi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_0[k]
                 + f_1 * fh0_0[k]
                 - f_2 * fh1_0[k]
                 + pb_x[k] * fi_0[k];

        t_1[k] = pa_y[k] * dk_0[k];

        t_2[k] = pa_z[k] * dk_0[k];

        t_3[k] = pa_x[k] * dk_1[k];

        t_4[k] = pa_x[k] * dk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, di_5, dk_1, dk_2, fh0_14, fh1_14, \
                         fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * di_5[k]
                 + f_1 * fh0_14[k]
                 - f_2 * fh1_14[k]
                 + pb_y[k] * fi_5[k];

        t_6[k] = pa_z[k] * dk_1[k];

        t_7[k] = pa_y[k] * dk_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, di_14, fh0_26, fh1_26, fi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_14[k]
                 + f_1 * fh0_26[k]
                 - f_2 * fh1_26[k]
                 + pb_z[k] * fi_8[k];
    }
}

}  // namespace simdt2ceri
