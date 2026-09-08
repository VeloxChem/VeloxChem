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
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_143 = buffer.data(dk + 143);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_215 = buffer.data(dk + 215);

    const auto *fh0_0 = buffer.data(fh0 + 0);
    const auto *fh0_1 = buffer.data(fh0 + 1);
    const auto *fh0_2 = buffer.data(fh0 + 2);
    const auto *fh0_3 = buffer.data(fh0 + 3);
    const auto *fh0_5 = buffer.data(fh0 + 5);
    const auto *fh0_6 = buffer.data(fh0 + 6);
    const auto *fh0_8 = buffer.data(fh0 + 8);
    const auto *fh0_9 = buffer.data(fh0 + 9);
    const auto *fh0_15 = buffer.data(fh0 + 15);
    const auto *fh0_17 = buffer.data(fh0 + 17);
    const auto *fh0_18 = buffer.data(fh0 + 18);
    const auto *fh0_19 = buffer.data(fh0 + 19);
    const auto *fh0_20 = buffer.data(fh0 + 20);
    const auto *fh0_126 = buffer.data(fh0 + 126);
    const auto *fh0_129 = buffer.data(fh0 + 129);
    const auto *fh0_131 = buffer.data(fh0 + 131);
    const auto *fh0_132 = buffer.data(fh0 + 132);
    const auto *fh0_135 = buffer.data(fh0 + 135);
    const auto *fh0_136 = buffer.data(fh0 + 136);
    const auto *fh0_138 = buffer.data(fh0 + 138);
    const auto *fh0_140 = buffer.data(fh0 + 140);
    const auto *fh0_141 = buffer.data(fh0 + 141);
    const auto *fh0_142 = buffer.data(fh0 + 142);
    const auto *fh0_143 = buffer.data(fh0 + 143);
    const auto *fh0_144 = buffer.data(fh0 + 144);
    const auto *fh0_146 = buffer.data(fh0 + 146);
    const auto *fh0_189 = buffer.data(fh0 + 189);
    const auto *fh0_192 = buffer.data(fh0 + 192);
    const auto *fh0_194 = buffer.data(fh0 + 194);
    const auto *fh0_195 = buffer.data(fh0 + 195);
    const auto *fh0_198 = buffer.data(fh0 + 198);
    const auto *fh0_199 = buffer.data(fh0 + 199);
    const auto *fh0_201 = buffer.data(fh0 + 201);
    const auto *fh0_203 = buffer.data(fh0 + 203);
    const auto *fh0_204 = buffer.data(fh0 + 204);
    const auto *fh0_206 = buffer.data(fh0 + 206);
    const auto *fh0_207 = buffer.data(fh0 + 207);
    const auto *fh0_208 = buffer.data(fh0 + 208);
    const auto *fh0_209 = buffer.data(fh0 + 209);

    const auto *fh1_0 = buffer.data(fh1 + 0);
    const auto *fh1_1 = buffer.data(fh1 + 1);
    const auto *fh1_2 = buffer.data(fh1 + 2);
    const auto *fh1_3 = buffer.data(fh1 + 3);
    const auto *fh1_5 = buffer.data(fh1 + 5);
    const auto *fh1_6 = buffer.data(fh1 + 6);
    const auto *fh1_8 = buffer.data(fh1 + 8);
    const auto *fh1_9 = buffer.data(fh1 + 9);
    const auto *fh1_15 = buffer.data(fh1 + 15);
    const auto *fh1_17 = buffer.data(fh1 + 17);
    const auto *fh1_18 = buffer.data(fh1 + 18);
    const auto *fh1_19 = buffer.data(fh1 + 19);
    const auto *fh1_20 = buffer.data(fh1 + 20);
    const auto *fh1_126 = buffer.data(fh1 + 126);
    const auto *fh1_129 = buffer.data(fh1 + 129);
    const auto *fh1_131 = buffer.data(fh1 + 131);
    const auto *fh1_132 = buffer.data(fh1 + 132);
    const auto *fh1_135 = buffer.data(fh1 + 135);
    const auto *fh1_136 = buffer.data(fh1 + 136);
    const auto *fh1_138 = buffer.data(fh1 + 138);
    const auto *fh1_140 = buffer.data(fh1 + 140);
    const auto *fh1_141 = buffer.data(fh1 + 141);
    const auto *fh1_142 = buffer.data(fh1 + 142);
    const auto *fh1_143 = buffer.data(fh1 + 143);
    const auto *fh1_144 = buffer.data(fh1 + 144);
    const auto *fh1_146 = buffer.data(fh1 + 146);
    const auto *fh1_189 = buffer.data(fh1 + 189);
    const auto *fh1_192 = buffer.data(fh1 + 192);
    const auto *fh1_194 = buffer.data(fh1 + 194);
    const auto *fh1_195 = buffer.data(fh1 + 195);
    const auto *fh1_198 = buffer.data(fh1 + 198);
    const auto *fh1_199 = buffer.data(fh1 + 199);
    const auto *fh1_201 = buffer.data(fh1 + 201);
    const auto *fh1_203 = buffer.data(fh1 + 203);
    const auto *fh1_204 = buffer.data(fh1 + 204);
    const auto *fh1_206 = buffer.data(fh1 + 206);
    const auto *fh1_207 = buffer.data(fh1 + 207);
    const auto *fh1_208 = buffer.data(fh1 + 208);
    const auto *fh1_209 = buffer.data(fh1 + 209);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

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
                         fh1_2, fh1_3, fi_3, fi_5, fi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * fh0_1[k]
                 - f_6 * fh1_1[k]
                 + pb_y[k] * fi_3[k];

        t_7[k] = pb_z[k] * fi_3[k];

        t_8[k] = pb_y[k] * fi_5[k];

        t_9[k] = f_5 * fh0_2[k]
                 - f_6 * fh1_2[k]
                 + pb_z[k] * fi_5[k];

        t_10[k] = f_7 * fh0_3[k]
                  - f_8 * fh1_3[k]
                  + pb_y[k] * fi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, fh0_5, fh0_6, fh1_5, \
                         fh1_6, fi_6, fi_8, fi_9, fi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * fi_6[k];

        t_12[k] = f_3 * fh0_5[k]
                  - f_4 * fh1_5[k]
                  + pb_y[k] * fi_8[k];

        t_13[k] = pb_y[k] * fi_9[k];

        t_14[k] = f_7 * fh0_5[k]
                  - f_8 * fh1_5[k]
                  + pb_z[k] * fi_9[k];

        t_15[k] = f_9 * fh0_6[k]
                  - f_10 * fh1_6[k]
                  + pb_y[k] * fi_10[k];

        t_16[k] = pb_z[k] * fi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, fh0_8, fh0_9, fh1_8, fh1_9, \
                         fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * fh0_8[k]
                  - f_6 * fh1_8[k]
                  + pb_y[k] * fi_12[k];

        t_18[k] = f_3 * fh0_9[k]
                  - f_4 * fh1_9[k]
                  + pb_y[k] * fi_13[k];

        t_19[k] = pb_y[k] * fi_14[k];

        t_20[k] = f_9 * fh0_9[k]
                  - f_10 * fh1_9[k]
                  + pb_z[k] * fi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, di_21, di_23, di_24, di_25, \
                         fi_15, fi_21, fi_23, fi_24, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * di_21[k]
                  + pb_x[k] * fi_21[k];

        t_22[k] = pb_z[k] * fi_15[k];

        t_23[k] = f_0 * di_23[k]
                  + pb_x[k] * fi_23[k];

        t_24[k] = f_0 * di_24[k]
                  + pb_x[k] * fi_24[k];

        t_25[k] = f_0 * di_25[k]
                  + pb_x[k] * fi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, di_27, fh0_15, fh1_15, \
                         fi_20, fi_21, fi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * fi_20[k];

        t_27[k] = f_0 * di_27[k]
                  + pb_x[k] * fi_27[k];

        t_28[k] = f_1 * fh0_15[k]
                  - f_2 * fh1_15[k]
                  + pb_y[k] * fi_21[k];

        t_29[k] = pb_z[k] * fi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, fh0_17, fh0_18, fh0_19, fh1_17, fh1_18, \
                         fh1_19, fi_23, fi_24, fi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * fh0_17[k]
                  - f_10 * fh1_17[k]
                  + pb_y[k] * fi_23[k];

        t_31[k] = f_7 * fh0_18[k]
                  - f_8 * fh1_18[k]
                  + pb_y[k] * fi_24[k];

        t_32[k] = f_5 * fh0_19[k]
                  - f_6 * fh1_19[k]
                  + pb_y[k] * fi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, di_0, dk_0, \
                         fh0_20, fh1_20, fi_26, fi_27, fi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * fh0_20[k]
                  - f_4 * fh1_20[k]
                  + pb_y[k] * fi_26[k];

        t_34[k] = pb_y[k] * fi_27[k];

        t_35[k] = f_1 * fh0_20[k]
                  - f_2 * fh1_20[k]
                  + pb_z[k] * fi_27[k];

        t_36[k] = pa_y[k] * dk_0[k];

        t_37[k] = f_11 * di_0[k]
                  + pb_y[k] * fi_28[k];

        t_38[k] = pb_z[k] * fi_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, di_1, di_3, dk_3, dk_5, \
                         dk_6, fi_29, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * di_1[k]
                  + pa_y[k] * dk_3[k];

        t_40[k] = pb_z[k] * fi_29[k];

        t_41[k] = pa_y[k] * dk_5[k];

        t_42[k] = f_0 * di_3[k]
                  + pa_y[k] * dk_6[k];

        t_43[k] = pb_z[k] * fi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, di_5, di_6, di_8, \
                         dk_9, dk_10, dk_12, fi_33, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * di_5[k]
                  + pb_y[k] * fi_33[k];

        t_45[k] = pa_y[k] * dk_9[k];

        t_46[k] = f_13 * di_6[k]
                  + pa_y[k] * dk_10[k];

        t_47[k] = pb_z[k] * fi_34[k];

        t_48[k] = f_12 * di_8[k]
                  + pa_y[k] * dk_12[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, di_9, di_10, di_12, \
                         dk_14, dk_15, dk_17, fi_37, fi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * di_9[k]
                  + pb_y[k] * fi_37[k];

        t_50[k] = pa_y[k] * dk_14[k];

        t_51[k] = f_14 * di_10[k]
                  + pa_y[k] * dk_15[k];

        t_52[k] = pb_z[k] * fi_38[k];

        t_53[k] = f_0 * di_12[k]
                  + pa_y[k] * dk_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, di_13, di_14, di_49, dk_18, \
                         dk_20, fi_42, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * di_13[k]
                  + pa_y[k] * dk_18[k];

        t_55[k] = f_11 * di_14[k]
                  + pb_y[k] * fi_42[k];

        t_56[k] = pa_y[k] * dk_20[k];

        t_57[k] = f_12 * di_49[k]
                  + pb_x[k] * fi_49[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, di_51, di_52, di_53, di_54, \
                         fi_43, fi_51, fi_52, fi_53, fi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * fi_43[k];

        t_59[k] = f_12 * di_51[k]
                  + pb_x[k] * fi_51[k];

        t_60[k] = f_12 * di_52[k]
                  + pb_x[k] * fi_52[k];

        t_61[k] = f_12 * di_53[k]
                  + pb_x[k] * fi_53[k];

        t_62[k] = f_12 * di_54[k]
                  + pb_x[k] * fi_54[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, di_21, di_23, di_24, dk_27, \
                         dk_28, dk_30, dk_31, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * dk_27[k];

        t_64[k] = f_15 * di_21[k]
                  + pa_y[k] * dk_28[k];

        t_65[k] = pb_z[k] * fi_49[k];

        t_66[k] = f_14 * di_23[k]
                  + pa_y[k] * dk_30[k];

        t_67[k] = f_13 * di_24[k]
                  + pa_y[k] * dk_31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, di_25, di_26, di_27, \
                         dk_0, dk_32, dk_33, dk_35, fi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * di_25[k]
                  + pa_y[k] * dk_32[k];

        t_69[k] = f_12 * di_26[k]
                  + pa_y[k] * dk_33[k];

        t_70[k] = f_11 * di_27[k]
                  + pb_y[k] * fi_55[k];

        t_71[k] = pa_y[k] * dk_35[k];

        t_72[k] = pa_z[k] * dk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, di_0, di_2, \
                         dk_3, dk_5, dk_6, fi_56, fi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * fi_56[k];

        t_74[k] = f_11 * di_0[k]
                  + pb_z[k] * fi_56[k];

        t_75[k] = pa_z[k] * dk_3[k];

        t_76[k] = pb_y[k] * fi_58[k];

        t_77[k] = f_12 * di_2[k]
                  + pa_z[k] * dk_5[k];

        t_78[k] = pa_z[k] * dk_6[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, di_3, di_5, di_6, \
                         dk_9, dk_10, fi_59, fi_61, fi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * di_3[k]
                  + pb_z[k] * fi_59[k];

        t_80[k] = pb_y[k] * fi_61[k];

        t_81[k] = f_0 * di_5[k]
                  + pa_z[k] * dk_9[k];

        t_82[k] = pa_z[k] * dk_10[k];

        t_83[k] = f_11 * di_6[k]
                  + pb_z[k] * fi_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, di_7, di_9, di_10, \
                         dk_12, dk_14, dk_15, fi_65, fi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * di_7[k]
                  + pa_z[k] * dk_12[k];

        t_85[k] = pb_y[k] * fi_65[k];

        t_86[k] = f_13 * di_9[k]
                  + pa_z[k] * dk_14[k];

        t_87[k] = pa_z[k] * dk_15[k];

        t_88[k] = f_11 * di_10[k]
                  + pb_z[k] * fi_66[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, di_11, di_12, di_14, dk_17, \
                         dk_18, dk_20, dk_21, fi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * di_11[k]
                  + pa_z[k] * dk_17[k];

        t_90[k] = f_0 * di_12[k]
                  + pa_z[k] * dk_18[k];

        t_91[k] = pb_y[k] * fi_70[k];

        t_92[k] = f_14 * di_14[k]
                  + pa_z[k] * dk_20[k];

        t_93[k] = pa_z[k] * dk_21[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, di_78, di_79, di_80, di_81, \
                         fi_76, fi_78, fi_79, fi_80, fi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_12 * di_78[k]
                  + pb_x[k] * fi_78[k];

        t_95[k] = f_12 * di_79[k]
                  + pb_x[k] * fi_79[k];

        t_96[k] = f_12 * di_80[k]
                  + pb_x[k] * fi_80[k];

        t_97[k] = f_12 * di_81[k]
                  + pb_x[k] * fi_81[k];

        t_98[k] = pb_y[k] * fi_76[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, di_21, di_22, di_83, \
                         dk_28, dk_30, fi_77, fi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_12 * di_83[k]
                  + pb_x[k] * fi_83[k];

        t_100[k] = pa_z[k] * dk_28[k];

        t_101[k] = f_11 * di_21[k]
                   + pb_z[k] * fi_77[k];

        t_102[k] = f_12 * di_22[k]
                   + pa_z[k] * dk_30[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, di_23, di_24, di_25, \
                         di_27, dk_31, dk_32, dk_33, dk_35, fi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * di_23[k]
                   + pa_z[k] * dk_31[k];

        t_104[k] = f_13 * di_24[k]
                   + pa_z[k] * dk_32[k];

        t_105[k] = f_14 * di_25[k]
                   + pa_z[k] * dk_33[k];

        t_106[k] = pb_y[k] * fi_83[k];

        t_107[k] = f_15 * di_27[k]
                   + pa_z[k] * dk_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_x, pb_y, pb_z, di_28, di_84, \
                         di_87, dk_108, dk_111, fi_84, fi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * di_84[k]
                   + pa_x[k] * dk_108[k];

        t_109[k] = f_12 * di_28[k]
                   + pb_y[k] * fi_84[k];

        t_110[k] = pb_z[k] * fi_84[k];

        t_111[k] = f_14 * di_87[k]
                   + pa_x[k] * dk_111[k];

        t_112[k] = pb_z[k] * fi_85[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pb_y, pb_z, di_33, di_89, di_90, \
                         dk_113, dk_114, fi_87, fi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * di_89[k]
                   + pa_x[k] * dk_113[k];

        t_114[k] = f_13 * di_90[k]
                   + pa_x[k] * dk_114[k];

        t_115[k] = pb_z[k] * fi_87[k];

        t_116[k] = f_12 * di_33[k]
                   + pb_y[k] * fi_89[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pb_z, di_93, di_94, di_96, dk_117, \
                         dk_118, dk_120, fi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_13 * di_93[k]
                   + pa_x[k] * dk_117[k];

        t_118[k] = f_0 * di_94[k]
                   + pa_x[k] * dk_118[k];

        t_119[k] = pb_z[k] * fi_90[k];

        t_120[k] = f_0 * di_96[k]
                   + pa_x[k] * dk_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pb_y, pb_z, di_37, di_98, di_99, \
                         dk_122, dk_123, fi_93, fi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * di_37[k]
                   + pb_y[k] * fi_93[k];

        t_122[k] = f_0 * di_98[k]
                   + pa_x[k] * dk_122[k];

        t_123[k] = f_12 * di_99[k]
                   + pa_x[k] * dk_123[k];

        t_124[k] = pb_z[k] * fi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_x, pb_y, di_42, di_101, di_102, \
                         di_104, dk_125, dk_126, dk_128, fi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * di_101[k]
                   + pa_x[k] * dk_125[k];

        t_126[k] = f_12 * di_102[k]
                   + pa_x[k] * dk_126[k];

        t_127[k] = f_12 * di_42[k]
                   + pb_y[k] * fi_98[k];

        t_128[k] = f_12 * di_104[k]
                   + pa_x[k] * dk_128[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, di_105, di_107, \
                         di_108, di_109, fi_99, fi_105, fi_107, fi_108, \
                         fi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * di_105[k]
                   + pb_x[k] * fi_105[k];

        t_130[k] = pb_z[k] * fi_99[k];

        t_131[k] = f_11 * di_107[k]
                   + pb_x[k] * fi_107[k];

        t_132[k] = f_11 * di_108[k]
                   + pb_x[k] * fi_108[k];

        t_133[k] = f_11 * di_109[k]
                   + pb_x[k] * fi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_x, pb_z, di_110, di_111, \
                         dk_136, dk_138, fi_105, fi_110, fi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * di_110[k]
                   + pb_x[k] * fi_110[k];

        t_135[k] = f_11 * di_111[k]
                   + pb_x[k] * fi_111[k];

        t_136[k] = pa_x[k] * dk_136[k];

        t_137[k] = pb_z[k] * fi_105[k];

        t_138[k] = pa_x[k] * dk_138[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, pa_x, pa_y, dk_72, dk_139, \
                         dk_140, dk_141, dk_142, dk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_x[k] * dk_139[k];

        t_140[k] = pa_x[k] * dk_140[k];

        t_141[k] = pa_x[k] * dk_141[k];

        t_142[k] = pa_x[k] * dk_142[k];

        t_143[k] = pa_x[k] * dk_143[k];

        t_144[k] = pa_y[k] * dk_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, di_58, \
                         dk_37, dk_39, dk_42, dk_74, dk_77, fi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * dk_37[k];

        t_146[k] = pa_y[k] * dk_74[k];

        t_147[k] = pa_z[k] * dk_39[k];

        t_148[k] = f_11 * di_58[k]
                   + pb_y[k] * fi_114[k];

        t_149[k] = pa_y[k] * dk_77[k];

        t_150[k] = pa_z[k] * dk_42[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, di_31, di_61, \
                         dk_46, dk_81, fi_115, fi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * di_31[k]
                   + pb_z[k] * fi_115[k];

        t_152[k] = f_11 * di_61[k]
                   + pb_y[k] * fi_117[k];

        t_153[k] = pa_y[k] * dk_81[k];

        t_154[k] = pa_z[k] * dk_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_x, pa_y, pb_y, pb_z, di_34, di_65, \
                         di_124, dk_86, dk_156, fi_118, fi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * di_34[k]
                   + pb_z[k] * fi_118[k];

        t_156[k] = f_0 * di_124[k]
                   + pa_x[k] * dk_156[k];

        t_157[k] = f_11 * di_65[k]
                   + pb_y[k] * fi_121[k];

        t_158[k] = pa_y[k] * dk_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pa_z, pb_z, di_38, di_129, di_130, \
                         dk_51, dk_161, dk_162, fi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * dk_51[k];

        t_160[k] = f_11 * di_38[k]
                   + pb_z[k] * fi_122[k];

        t_161[k] = f_12 * di_129[k]
                   + pa_x[k] * dk_161[k];

        t_162[k] = f_12 * di_130[k]
                   + pa_x[k] * dk_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, di_70, di_134, \
                         dk_57, dk_92, fi_126, fi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * di_70[k]
                   + pb_y[k] * fi_126[k];

        t_164[k] = pa_y[k] * dk_92[k];

        t_165[k] = pa_z[k] * dk_57[k];

        t_166[k] = f_11 * di_134[k]
                   + pb_x[k] * fi_134[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, di_135, di_136, \
                         di_137, di_138, dk_99, fi_135, fi_136, fi_137, \
                         fi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_11 * di_135[k]
                   + pb_x[k] * fi_135[k];

        t_168[k] = f_11 * di_136[k]
                   + pb_x[k] * fi_136[k];

        t_169[k] = f_11 * di_137[k]
                   + pb_x[k] * fi_137[k];

        t_170[k] = f_11 * di_138[k]
                   + pb_x[k] * fi_138[k];

        t_171[k] = pa_y[k] * dk_99[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, t_177, t_178, pa_x, dk_172, \
                         dk_173, dk_174, dk_175, dk_176, dk_177, \
                         dk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_x[k] * dk_172[k];

        t_173[k] = pa_x[k] * dk_173[k];

        t_174[k] = pa_x[k] * dk_174[k];

        t_175[k] = pa_x[k] * dk_175[k];

        t_176[k] = pa_x[k] * dk_176[k];

        t_177[k] = pa_x[k] * dk_177[k];

        t_178[k] = pa_x[k] * dk_178[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pa_x, pb_y, pb_z, di_56, di_140, \
                         di_143, dk_179, dk_180, dk_183, fi_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_x[k] * dk_179[k];

        t_180[k] = f_15 * di_140[k]
                   + pa_x[k] * dk_180[k];

        t_181[k] = pb_y[k] * fi_140[k];

        t_182[k] = f_12 * di_56[k]
                   + pb_z[k] * fi_140[k];

        t_183[k] = f_14 * di_143[k]
                   + pa_x[k] * dk_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pa_x, pb_y, pb_z, di_59, di_145, \
                         di_146, dk_185, dk_186, fi_142, fi_143, \
                         fi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * fi_142[k];

        t_185[k] = f_14 * di_145[k]
                   + pa_x[k] * dk_185[k];

        t_186[k] = f_13 * di_146[k]
                   + pa_x[k] * dk_186[k];

        t_187[k] = f_12 * di_59[k]
                   + pb_z[k] * fi_143[k];

        t_188[k] = pb_y[k] * fi_145[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pb_z, di_62, di_149, di_150, \
                         di_152, dk_189, dk_190, dk_192, fi_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * di_149[k]
                   + pa_x[k] * dk_189[k];

        t_190[k] = f_0 * di_150[k]
                   + pa_x[k] * dk_190[k];

        t_191[k] = f_12 * di_62[k]
                   + pb_z[k] * fi_146[k];

        t_192[k] = f_0 * di_152[k]
                   + pa_x[k] * dk_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_x, pb_y, pb_z, di_66, di_154, di_155, \
                         dk_194, dk_195, fi_149, fi_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_y[k] * fi_149[k];

        t_194[k] = f_0 * di_154[k]
                   + pa_x[k] * dk_194[k];

        t_195[k] = f_12 * di_155[k]
                   + pa_x[k] * dk_195[k];

        t_196[k] = f_12 * di_66[k]
                   + pb_z[k] * fi_150[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_x, pb_y, di_157, di_158, di_160, \
                         dk_197, dk_198, dk_200, fi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * di_157[k]
                   + pa_x[k] * dk_197[k];

        t_198[k] = f_12 * di_158[k]
                   + pa_x[k] * dk_198[k];

        t_199[k] = pb_y[k] * fi_154[k];

        t_200[k] = f_12 * di_160[k]
                   + pa_x[k] * dk_200[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pb_x, di_161, di_162, di_163, \
                         di_164, di_165, fi_161, fi_162, fi_163, fi_164, \
                         fi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_11 * di_161[k]
                   + pb_x[k] * fi_161[k];

        t_202[k] = f_11 * di_162[k]
                   + pb_x[k] * fi_162[k];

        t_203[k] = f_11 * di_163[k]
                   + pb_x[k] * fi_163[k];

        t_204[k] = f_11 * di_164[k]
                   + pb_x[k] * fi_164[k];

        t_205[k] = f_11 * di_165[k]
                   + pb_x[k] * fi_165[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, t_211, pa_x, pb_x, pb_y, di_167, \
                         dk_208, dk_209, dk_210, dk_211, fi_160, \
                         fi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_y[k] * fi_160[k];

        t_207[k] = f_11 * di_167[k]
                   + pb_x[k] * fi_167[k];

        t_208[k] = pa_x[k] * dk_208[k];

        t_209[k] = pa_x[k] * dk_209[k];

        t_210[k] = pa_x[k] * dk_210[k];

        t_211[k] = pa_x[k] * dk_211[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pb_x, pb_y, dk_212, dk_213, \
                         dk_215, fh0_126, fh1_126, fi_167, fi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pa_x[k] * dk_212[k];

        t_213[k] = pa_x[k] * dk_213[k];

        t_214[k] = pb_y[k] * fi_167[k];

        t_215[k] = pa_x[k] * dk_215[k];

        t_216[k] = f_1 * fh0_126[k]
                   - f_2 * fh1_126[k]
                   + pb_x[k] * fi_168[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_x, pb_y, pb_z, di_84, fh0_129, \
                         fh1_129, fi_168, fi_169, fi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_0 * di_84[k]
                   + pb_y[k] * fi_168[k];

        t_218[k] = pb_z[k] * fi_168[k];

        t_219[k] = f_9 * fh0_129[k]
                   - f_10 * fh1_129[k]
                   + pb_x[k] * fi_171[k];

        t_220[k] = pb_z[k] * fi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_x, pb_y, pb_z, di_89, fh0_131, \
                         fh0_132, fh1_131, fh1_132, fi_171, fi_173, \
                         fi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * fh0_131[k]
                   - f_10 * fh1_131[k]
                   + pb_x[k] * fi_173[k];

        t_222[k] = f_7 * fh0_132[k]
                   - f_8 * fh1_132[k]
                   + pb_x[k] * fi_174[k];

        t_223[k] = pb_z[k] * fi_171[k];

        t_224[k] = f_0 * di_89[k]
                   + pb_y[k] * fi_173[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_x, pb_z, fh0_135, fh0_136, fh0_138, \
                         fh1_135, fh1_136, fh1_138, fi_174, fi_177, fi_178, \
                         fi_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_7 * fh0_135[k]
                   - f_8 * fh1_135[k]
                   + pb_x[k] * fi_177[k];

        t_226[k] = f_5 * fh0_136[k]
                   - f_6 * fh1_136[k]
                   + pb_x[k] * fi_178[k];

        t_227[k] = pb_z[k] * fi_174[k];

        t_228[k] = f_5 * fh0_138[k]
                   - f_6 * fh1_138[k]
                   + pb_x[k] * fi_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, di_93, fh0_140, \
                         fh0_141, fh1_140, fh1_141, fi_177, fi_178, fi_182, \
                         fi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * di_93[k]
                   + pb_y[k] * fi_177[k];

        t_230[k] = f_5 * fh0_140[k]
                   - f_6 * fh1_140[k]
                   + pb_x[k] * fi_182[k];

        t_231[k] = f_3 * fh0_141[k]
                   - f_4 * fh1_141[k]
                   + pb_x[k] * fi_183[k];

        t_232[k] = pb_z[k] * fi_178[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_x, pb_y, di_98, fh0_143, fh0_144, fh1_143, \
                         fh1_144, fi_182, fi_185, fi_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * fh0_143[k]
                   - f_4 * fh1_143[k]
                   + pb_x[k] * fi_185[k];

        t_234[k] = f_3 * fh0_144[k]
                   - f_4 * fh1_144[k]
                   + pb_x[k] * fi_186[k];

        t_235[k] = f_0 * di_98[k]
                   + pb_y[k] * fi_182[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, pb_x, fh0_146, fh1_146, \
                         fi_188, fi_189, fi_190, fi_191, fi_192, \
                         fi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_3 * fh0_146[k]
                   - f_4 * fh1_146[k]
                   + pb_x[k] * fi_188[k];

        t_237[k] = pb_x[k] * fi_189[k];

        t_238[k] = pb_x[k] * fi_190[k];

        t_239[k] = pb_x[k] * fi_191[k];

        t_240[k] = pb_x[k] * fi_192[k];

        t_241[k] = pb_x[k] * fi_193[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, pb_y, pb_z, di_105, fh0_141, \
                         fh1_141, fi_189, fi_190, fi_194, fi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = pb_x[k] * fi_194[k];

        t_243[k] = pb_x[k] * fi_195[k];

        t_244[k] = f_0 * di_105[k]
                   + f_1 * fh0_141[k]
                   - f_2 * fh1_141[k]
                   + pb_y[k] * fi_189[k];

        t_245[k] = pb_z[k] * fi_189[k];

        t_246[k] = f_3 * fh0_141[k]
                   - f_4 * fh1_141[k]
                   + pb_z[k] * fi_190[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pb_z, fh0_142, fh0_143, fh0_144, fh1_142, \
                         fh1_143, fh1_144, fi_191, fi_192, fi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_5 * fh0_142[k]
                   - f_6 * fh1_142[k]
                   + pb_z[k] * fi_191[k];

        t_248[k] = f_7 * fh0_143[k]
                   - f_8 * fh1_143[k]
                   + pb_z[k] * fi_192[k];

        t_249[k] = f_9 * fh0_144[k]
                   - f_10 * fh1_144[k]
                   + pb_z[k] * fi_193[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, pa_z, pb_y, pb_z, di_84, di_111, \
                         dk_108, dk_109, fh0_146, fh1_146, fi_195, \
                         fi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * di_111[k]
                   + pb_y[k] * fi_195[k];

        t_251[k] = f_1 * fh0_146[k]
                   - f_2 * fh1_146[k]
                   + pb_z[k] * fi_195[k];

        t_252[k] = pa_z[k] * dk_108[k];

        t_253[k] = pa_z[k] * dk_109[k];

        t_254[k] = f_11 * di_84[k]
                   + pb_z[k] * fi_196[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pa_z, pb_y, pb_z, di_86, di_87, \
                         di_114, dk_111, dk_113, dk_114, fi_198, \
                         fi_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pa_z[k] * dk_111[k];

        t_256[k] = f_12 * di_114[k]
                   + pb_y[k] * fi_198[k];

        t_257[k] = f_12 * di_86[k]
                   + pa_z[k] * dk_113[k];

        t_258[k] = pa_z[k] * dk_114[k];

        t_259[k] = f_11 * di_87[k]
                   + pb_z[k] * fi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_z, pb_y, pb_z, di_89, di_90, di_117, \
                         dk_117, dk_118, fi_201, fi_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_12 * di_117[k]
                   + pb_y[k] * fi_201[k];

        t_261[k] = f_0 * di_89[k]
                   + pa_z[k] * dk_117[k];

        t_262[k] = pa_z[k] * dk_118[k];

        t_263[k] = f_11 * di_90[k]
                   + pb_z[k] * fi_202[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_z, pb_y, di_91, di_93, di_121, dk_120, \
                         dk_122, dk_123, fi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_12 * di_91[k]
                   + pa_z[k] * dk_120[k];

        t_265[k] = f_12 * di_121[k]
                   + pb_y[k] * fi_205[k];

        t_266[k] = f_13 * di_93[k]
                   + pa_z[k] * dk_122[k];

        t_267[k] = pa_z[k] * dk_123[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_z, pb_y, pb_z, di_94, di_95, di_96, \
                         di_126, dk_125, dk_126, fi_206, fi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_11 * di_94[k]
                   + pb_z[k] * fi_206[k];

        t_269[k] = f_12 * di_95[k]
                   + pa_z[k] * dk_125[k];

        t_270[k] = f_0 * di_96[k]
                   + pa_z[k] * dk_126[k];

        t_271[k] = f_12 * di_126[k]
                   + pb_y[k] * fi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, t_277, pa_z, pb_x, di_98, dk_128, \
                         fi_217, fi_218, fi_219, fi_220, fi_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_14 * di_98[k]
                   + pa_z[k] * dk_128[k];

        t_273[k] = pb_x[k] * fi_217[k];

        t_274[k] = pb_x[k] * fi_218[k];

        t_275[k] = pb_x[k] * fi_219[k];

        t_276[k] = pb_x[k] * fi_220[k];

        t_277[k] = pb_x[k] * fi_221[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, pa_z, pb_x, pb_z, di_105, di_106, \
                         dk_136, dk_138, fi_217, fi_222, fi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pb_x[k] * fi_222[k];

        t_279[k] = pb_x[k] * fi_223[k];

        t_280[k] = pa_z[k] * dk_136[k];

        t_281[k] = f_11 * di_105[k]
                   + pb_z[k] * fi_217[k];

        t_282[k] = f_12 * di_106[k]
                   + pa_z[k] * dk_138[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_z, pb_y, di_107, di_108, di_109, \
                         di_139, dk_139, dk_140, dk_141, fi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_0 * di_107[k]
                   + pa_z[k] * dk_139[k];

        t_284[k] = f_13 * di_108[k]
                   + pa_z[k] * dk_140[k];

        t_285[k] = f_14 * di_109[k]
                   + pa_z[k] * dk_141[k];

        t_286[k] = f_12 * di_139[k]
                   + pb_y[k] * fi_223[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_y, pa_z, pb_y, di_111, di_140, \
                         di_141, dk_143, dk_180, dk_182, dk_183, \
                         fi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_15 * di_111[k]
                   + pa_z[k] * dk_143[k];

        t_288[k] = pa_y[k] * dk_180[k];

        t_289[k] = f_11 * di_140[k]
                   + pb_y[k] * fi_224[k];

        t_290[k] = pa_y[k] * dk_182[k];

        t_291[k] = f_12 * di_141[k]
                   + pa_y[k] * dk_183[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_y, pb_z, di_115, di_142, di_143, \
                         dk_185, dk_186, fi_226, fi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_11 * di_142[k]
                   + pb_y[k] * fi_226[k];

        t_293[k] = pa_y[k] * dk_185[k];

        t_294[k] = f_0 * di_143[k]
                   + pa_y[k] * dk_186[k];

        t_295[k] = f_12 * di_115[k]
                   + pb_z[k] * fi_227[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_y, pb_z, di_118, di_145, di_146, \
                         dk_189, dk_190, fi_229, fi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_11 * di_145[k]
                   + pb_y[k] * fi_229[k];

        t_297[k] = pa_y[k] * dk_189[k];

        t_298[k] = f_13 * di_146[k]
                   + pa_y[k] * dk_190[k];

        t_299[k] = f_12 * di_118[k]
                   + pb_z[k] * fi_230[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_y, pb_y, di_148, di_149, di_150, \
                         dk_192, dk_194, dk_195, fi_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_12 * di_148[k]
                   + pa_y[k] * dk_192[k];

        t_301[k] = f_11 * di_149[k]
                   + pb_y[k] * fi_233[k];

        t_302[k] = pa_y[k] * dk_194[k];

        t_303[k] = f_14 * di_150[k]
                   + pa_y[k] * dk_195[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_y, pb_y, pb_z, di_122, di_152, di_153, \
                         di_154, dk_197, dk_198, fi_234, fi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_12 * di_122[k]
                   + pb_z[k] * fi_234[k];

        t_305[k] = f_0 * di_152[k]
                   + pa_y[k] * dk_197[k];

        t_306[k] = f_12 * di_153[k]
                   + pa_y[k] * dk_198[k];

        t_307[k] = f_11 * di_154[k]
                   + pb_y[k] * fi_238[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, t_314, pa_y, pb_x, dk_200, \
                         fi_245, fi_246, fi_247, fi_248, fi_249, \
                         fi_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_y[k] * dk_200[k];

        t_309[k] = pb_x[k] * fi_245[k];

        t_310[k] = pb_x[k] * fi_246[k];

        t_311[k] = pb_x[k] * fi_247[k];

        t_312[k] = pb_x[k] * fi_248[k];

        t_313[k] = pb_x[k] * fi_249[k];

        t_314[k] = pb_x[k] * fi_250[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_y, pb_x, pb_z, di_133, di_161, di_163, \
                         dk_208, dk_210, fi_245, fi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_x[k] * fi_251[k];

        t_316[k] = f_15 * di_161[k]
                   + pa_y[k] * dk_208[k];

        t_317[k] = f_12 * di_133[k]
                   + pb_z[k] * fi_245[k];

        t_318[k] = f_14 * di_163[k]
                   + pa_y[k] * dk_210[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, pa_y, pb_y, di_164, di_165, \
                         di_166, di_167, dk_211, dk_212, dk_213, dk_215, \
                         fi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_13 * di_164[k]
                   + pa_y[k] * dk_211[k];

        t_320[k] = f_0 * di_165[k]
                   + pa_y[k] * dk_212[k];

        t_321[k] = f_12 * di_166[k]
                   + pa_y[k] * dk_213[k];

        t_322[k] = f_11 * di_167[k]
                   + pb_y[k] * fi_251[k];

        t_323[k] = pa_y[k] * dk_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pb_x, pb_y, pb_z, di_140, fh0_189, \
                         fh0_192, fh1_189, fh1_192, fi_252, fi_254, \
                         fi_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * fh0_189[k]
                   - f_2 * fh1_189[k]
                   + pb_x[k] * fi_252[k];

        t_325[k] = pb_y[k] * fi_252[k];

        t_326[k] = f_0 * di_140[k]
                   + pb_z[k] * fi_252[k];

        t_327[k] = f_9 * fh0_192[k]
                   - f_10 * fh1_192[k]
                   + pb_x[k] * fi_255[k];

        t_328[k] = pb_y[k] * fi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pb_x, pb_y, pb_z, di_143, fh0_194, \
                         fh0_195, fh1_194, fh1_195, fi_255, fi_257, \
                         fi_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_9 * fh0_194[k]
                   - f_10 * fh1_194[k]
                   + pb_x[k] * fi_257[k];

        t_330[k] = f_7 * fh0_195[k]
                   - f_8 * fh1_195[k]
                   + pb_x[k] * fi_258[k];

        t_331[k] = f_0 * di_143[k]
                   + pb_z[k] * fi_255[k];

        t_332[k] = pb_y[k] * fi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pb_x, pb_z, di_146, fh0_198, fh0_199, fh1_198, \
                         fh1_199, fi_258, fi_261, fi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_7 * fh0_198[k]
                   - f_8 * fh1_198[k]
                   + pb_x[k] * fi_261[k];

        t_334[k] = f_5 * fh0_199[k]
                   - f_6 * fh1_199[k]
                   + pb_x[k] * fi_262[k];

        t_335[k] = f_0 * di_146[k]
                   + pb_z[k] * fi_258[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pb_x, pb_y, fh0_201, fh0_203, fh0_204, \
                         fh1_201, fh1_203, fh1_204, fi_261, fi_264, fi_266, \
                         fi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_5 * fh0_201[k]
                   - f_6 * fh1_201[k]
                   + pb_x[k] * fi_264[k];

        t_337[k] = pb_y[k] * fi_261[k];

        t_338[k] = f_5 * fh0_203[k]
                   - f_6 * fh1_203[k]
                   + pb_x[k] * fi_266[k];

        t_339[k] = f_3 * fh0_204[k]
                   - f_4 * fh1_204[k]
                   + pb_x[k] * fi_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pb_x, pb_y, pb_z, di_150, fh0_206, \
                         fh0_207, fh1_206, fh1_207, fi_262, fi_266, fi_269, \
                         fi_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_0 * di_150[k]
                   + pb_z[k] * fi_262[k];

        t_341[k] = f_3 * fh0_206[k]
                   - f_4 * fh1_206[k]
                   + pb_x[k] * fi_269[k];

        t_342[k] = f_3 * fh0_207[k]
                   - f_4 * fh1_207[k]
                   + pb_x[k] * fi_270[k];

        t_343[k] = pb_y[k] * fi_266[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, t_349, pb_x, fh0_209, fh1_209, \
                         fi_272, fi_273, fi_274, fi_275, fi_276, \
                         fi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_3 * fh0_209[k]
                   - f_4 * fh1_209[k]
                   + pb_x[k] * fi_272[k];

        t_345[k] = pb_x[k] * fi_273[k];

        t_346[k] = pb_x[k] * fi_274[k];

        t_347[k] = pb_x[k] * fi_275[k];

        t_348[k] = pb_x[k] * fi_276[k];

        t_349[k] = pb_x[k] * fi_277[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_x, pb_y, pb_z, di_161, fh0_204, \
                         fh1_204, fi_273, fi_278, fi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_x[k] * fi_278[k];

        t_351[k] = pb_x[k] * fi_279[k];

        t_352[k] = f_1 * fh0_204[k]
                   - f_2 * fh1_204[k]
                   + pb_y[k] * fi_273[k];

        t_353[k] = f_0 * di_161[k]
                   + pb_z[k] * fi_273[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pb_y, fh0_206, fh0_207, fh0_208, fh1_206, \
                         fh1_207, fh1_208, fi_275, fi_276, fi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_9 * fh0_206[k]
                   - f_10 * fh1_206[k]
                   + pb_y[k] * fi_275[k];

        t_355[k] = f_7 * fh0_207[k]
                   - f_8 * fh1_207[k]
                   + pb_y[k] * fi_276[k];

        t_356[k] = f_5 * fh0_208[k]
                   - f_6 * fh1_208[k]
                   + pb_y[k] * fi_277[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pb_y, pb_z, di_167, fh0_209, fh1_209, fi_278, \
                         fi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_3 * fh0_209[k]
                   - f_4 * fh1_209[k]
                   + pb_y[k] * fi_278[k];

        t_358[k] = pb_y[k] * fi_279[k];

        t_359[k] = f_0 * di_167[k]
                   + f_1 * fh0_209[k]
                   - f_2 * fh1_209[k]
                   + pb_z[k] * fi_279[k];
    }
}

}  // namespace simdt2ceri
