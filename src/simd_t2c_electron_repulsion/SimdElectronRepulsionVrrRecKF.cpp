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


#include "SimdElectronRepulsionVrrRecKF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_kf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 2.0 / alpha;
    const auto f_11 = 2.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.0 / p;
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);

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

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_30 = buffer.data(hf0 + 30);
    const auto *hf0_36 = buffer.data(hf0 + 36);
    const auto *hf0_50 = buffer.data(hf0 + 50);
    const auto *hf0_59 = buffer.data(hf0 + 59);
    const auto *hf0_60 = buffer.data(hf0 + 60);
    const auto *hf0_66 = buffer.data(hf0 + 66);
    const auto *hf0_80 = buffer.data(hf0 + 80);
    const auto *hf0_90 = buffer.data(hf0 + 90);
    const auto *hf0_99 = buffer.data(hf0 + 99);
    const auto *hf0_106 = buffer.data(hf0 + 106);
    const auto *hf0_126 = buffer.data(hf0 + 126);
    const auto *hf0_129 = buffer.data(hf0 + 129);
    const auto *hf0_149 = buffer.data(hf0 + 149);
    const auto *hf0_156 = buffer.data(hf0 + 156);
    const auto *hf0_166 = buffer.data(hf0 + 166);
    const auto *hf0_176 = buffer.data(hf0 + 176);
    const auto *hf0_179 = buffer.data(hf0 + 179);
    const auto *hf0_186 = buffer.data(hf0 + 186);
    const auto *hf0_189 = buffer.data(hf0 + 189);
    const auto *hf0_199 = buffer.data(hf0 + 199);
    const auto *hf0_209 = buffer.data(hf0 + 209);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_30 = buffer.data(hf1 + 30);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_50 = buffer.data(hf1 + 50);
    const auto *hf1_59 = buffer.data(hf1 + 59);
    const auto *hf1_60 = buffer.data(hf1 + 60);
    const auto *hf1_66 = buffer.data(hf1 + 66);
    const auto *hf1_80 = buffer.data(hf1 + 80);
    const auto *hf1_90 = buffer.data(hf1 + 90);
    const auto *hf1_99 = buffer.data(hf1 + 99);
    const auto *hf1_106 = buffer.data(hf1 + 106);
    const auto *hf1_126 = buffer.data(hf1 + 126);
    const auto *hf1_129 = buffer.data(hf1 + 129);
    const auto *hf1_149 = buffer.data(hf1 + 149);
    const auto *hf1_156 = buffer.data(hf1 + 156);
    const auto *hf1_166 = buffer.data(hf1 + 166);
    const auto *hf1_176 = buffer.data(hf1 + 176);
    const auto *hf1_179 = buffer.data(hf1 + 179);
    const auto *hf1_186 = buffer.data(hf1 + 186);
    const auto *hf1_189 = buffer.data(hf1 + 189);
    const auto *hf1_199 = buffer.data(hf1 + 199);
    const auto *hf1_209 = buffer.data(hf1 + 209);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_167 = buffer.data(id + 167);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__211 = buffer.data(if_ + 211);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__272 = buffer.data(if_ + 272);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__277 = buffer.data(if_ + 277);
    const auto *if__279 = buffer.data(if_ + 279);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_20 = buffer.data(kp0 + 20);
    const auto *kp0_28 = buffer.data(kp0 + 28);
    const auto *kp0_32 = buffer.data(kp0 + 32);
    const auto *kp0_43 = buffer.data(kp0 + 43);
    const auto *kp0_47 = buffer.data(kp0 + 47);
    const auto *kp0_61 = buffer.data(kp0 + 61);
    const auto *kp0_84 = buffer.data(kp0 + 84);
    const auto *kp0_85 = buffer.data(kp0 + 85);
    const auto *kp0_86 = buffer.data(kp0 + 86);
    const auto *kp0_90 = buffer.data(kp0 + 90);
    const auto *kp0_93 = buffer.data(kp0 + 93);
    const auto *kp0_96 = buffer.data(kp0 + 96);
    const auto *kp0_99 = buffer.data(kp0 + 99);
    const auto *kp0_105 = buffer.data(kp0 + 105);
    const auto *kp0_106 = buffer.data(kp0 + 106);
    const auto *kp0_107 = buffer.data(kp0 + 107);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_20 = buffer.data(kp1 + 20);
    const auto *kp1_28 = buffer.data(kp1 + 28);
    const auto *kp1_32 = buffer.data(kp1 + 32);
    const auto *kp1_43 = buffer.data(kp1 + 43);
    const auto *kp1_47 = buffer.data(kp1 + 47);
    const auto *kp1_61 = buffer.data(kp1 + 61);
    const auto *kp1_84 = buffer.data(kp1 + 84);
    const auto *kp1_85 = buffer.data(kp1 + 85);
    const auto *kp1_86 = buffer.data(kp1 + 86);
    const auto *kp1_90 = buffer.data(kp1 + 90);
    const auto *kp1_93 = buffer.data(kp1 + 93);
    const auto *kp1_96 = buffer.data(kp1 + 96);
    const auto *kp1_99 = buffer.data(kp1 + 99);
    const auto *kp1_105 = buffer.data(kp1 + 105);
    const auto *kp1_106 = buffer.data(kp1 + 106);
    const auto *kp1_107 = buffer.data(kp1 + 107);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, id_3, kp0_0, kp1_0, \
                         kd_0, kd_2, kd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_0 * id_3[k]
                 + pb_x[k] * kd_3[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, id_5, kp0_1, kp0_2, kp1_1, \
                         kp1_2, kd_3, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * id_5[k]
                 + pb_x[k] * kd_5[k];

        t_6[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_3[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = pb_y[k] * kd_5[k];

        t_9[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, id_0, id_9, \
                         if__0, kd_6, kd_7, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * if__0[k];

        t_11[k] = f_3 * id_0[k]
                  + pb_y[k] * kd_6[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_4 * id_9[k]
                  + pb_x[k] * kd_9[k];

        t_14[k] = pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, id_3, id_5, if__5, \
                         if__6, if__9, kd_9, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * if__5[k];

        t_16[k] = f_5 * id_3[k]
                  + pa_y[k] * if__6[k];

        t_17[k] = pb_z[k] * kd_9[k];

        t_18[k] = f_3 * id_5[k]
                  + pb_y[k] * kd_11[k];

        t_19[k] = pa_y[k] * if__9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, id_0, if__0, if__3, \
                         kd_12, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * if__0[k];

        t_21[k] = pb_y[k] * kd_12[k];

        t_22[k] = f_3 * id_0[k]
                  + pb_z[k] * kd_12[k];

        t_23[k] = pa_z[k] * if__3[k];

        t_24[k] = pb_y[k] * kd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, id_3, id_5, \
                         id_17, if__6, if__9, kd_15, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * id_17[k]
                  + pb_x[k] * kd_17[k];

        t_26[k] = pa_z[k] * if__6[k];

        t_27[k] = f_3 * id_3[k]
                  + pb_z[k] * kd_15[k];

        t_28[k] = pb_y[k] * kd_17[k];

        t_29[k] = f_5 * id_5[k]
                  + pa_z[k] * if__9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, hf0_0, hf1_0, id_6, \
                         id_21, if__10, kd_18, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_y[k] * if__10[k];

        t_31[k] = f_8 * id_6[k]
                  + pb_y[k] * kd_18[k];

        t_32[k] = pb_z[k] * kd_18[k];

        t_33[k] = f_9 * id_21[k]
                  + pb_x[k] * kd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_36, hf1_36, id_23, \
                         if__36, kd_19, kd_21, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * kd_19[k];

        t_35[k] = f_9 * id_23[k]
                  + pb_x[k] * kd_23[k];

        t_36[k] = f_10 * hf0_36[k]
                  - f_11 * hf1_36[k]
                  + pa_x[k] * if__36[k];

        t_37[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, id_11, if__11, \
                         if__20, if__22, kp0_11, kp1_11, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * id_11[k]
                  + pb_y[k] * kd_23[k];

        t_39[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_z[k] * kd_23[k];

        t_40[k] = pa_y[k] * if__20[k];

        t_41[k] = pa_z[k] * if__11[k];

        t_42[k] = pa_y[k] * if__22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, id_9, id_28, \
                         if__13, if__16, if__25, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * if__13[k];

        t_44[k] = f_9 * id_28[k]
                  + pb_x[k] * kd_28[k];

        t_45[k] = pa_y[k] * if__25[k];

        t_46[k] = pa_z[k] * if__16[k];

        t_47[k] = f_3 * id_9[k]
                  + pb_z[k] * kd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, hf0_0, hf1_0, id_17, \
                         if__20, if__29, kd_29, kd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * id_17[k]
                  + pb_y[k] * kd_29[k];

        t_49[k] = pa_y[k] * if__29[k];

        t_50[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_z[k] * if__20[k];

        t_51[k] = pb_y[k] * kd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, id_12, id_33, id_35, kd_30, \
                         kd_32, kd_33, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * id_12[k]
                  + pb_z[k] * kd_30[k];

        t_53[k] = f_9 * id_33[k]
                  + pb_x[k] * kd_33[k];

        t_54[k] = pb_y[k] * kd_32[k];

        t_55[k] = f_9 * id_35[k]
                  + pb_x[k] * kd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, hf0_59, hf1_59, id_15, \
                         if__59, kp0_16, kp1_16, kd_33, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_y[k] * kd_33[k];

        t_57[k] = f_8 * id_15[k]
                  + pb_z[k] * kd_33[k];

        t_58[k] = pb_y[k] * kd_35[k];

        t_59[k] = f_10 * hf0_59[k]
                  - f_11 * hf1_59[k]
                  + pa_x[k] * if__59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, hf0_10, hf1_10, \
                         id_18, id_39, if__30, kd_36, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * hf0_10[k]
                  - f_13 * hf1_10[k]
                  + pa_y[k] * if__30[k];

        t_61[k] = f_5 * id_18[k]
                  + pb_y[k] * kd_36[k];

        t_62[k] = pb_z[k] * kd_36[k];

        t_63[k] = f_14 * id_39[k]
                  + pb_x[k] * kd_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, hf0_66, hf1_66, id_41, \
                         if__66, kd_37, kd_39, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * kd_37[k];

        t_65[k] = f_14 * id_41[k]
                  + pb_x[k] * kd_41[k];

        t_66[k] = f_15 * hf0_66[k]
                  - f_16 * hf1_66[k]
                  + pa_x[k] * if__66[k];

        t_67[k] = pb_z[k] * kd_39[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, id_18, id_23, if__30, \
                         if__31, kp0_20, kp1_20, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * id_23[k]
                  + pb_y[k] * kd_41[k];

        t_69[k] = f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_41[k];

        t_70[k] = pa_z[k] * if__30[k];

        t_71[k] = pa_z[k] * if__31[k];

        t_72[k] = f_3 * id_18[k]
                  + pb_z[k] * kd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, id_21, id_46, id_47, \
                         if__33, if__36, kd_45, kd_46, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * if__33[k];

        t_74[k] = f_14 * id_46[k]
                  + pb_x[k] * kd_46[k];

        t_75[k] = f_14 * id_47[k]
                  + pb_x[k] * kd_47[k];

        t_76[k] = pa_z[k] * if__36[k];

        t_77[k] = f_3 * id_21[k]
                  + pb_z[k] * kd_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, id_23, id_29, id_30, \
                         if__39, if__50, if__52, kd_47, kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * id_29[k]
                  + pb_y[k] * kd_47[k];

        t_79[k] = f_5 * id_23[k]
                  + pa_z[k] * if__39[k];

        t_80[k] = pa_y[k] * if__50[k];

        t_81[k] = f_3 * id_30[k]
                  + pb_y[k] * kd_48[k];

        t_82[k] = pa_y[k] * if__52[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, id_27, id_33, id_51, \
                         id_52, if__55, if__56, kd_51, kd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_14 * id_51[k]
                  + pb_x[k] * kd_51[k];

        t_84[k] = f_14 * id_52[k]
                  + pb_x[k] * kd_52[k];

        t_85[k] = pa_y[k] * if__55[k];

        t_86[k] = f_5 * id_33[k]
                  + pa_y[k] * if__56[k];

        t_87[k] = f_8 * id_27[k]
                  + pb_z[k] * kd_51[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, hf0_20, hf1_20, id_35, \
                         if__50, if__59, kd_53, kd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * id_35[k]
                  + pb_y[k] * kd_53[k];

        t_89[k] = pa_y[k] * if__59[k];

        t_90[k] = f_12 * hf0_20[k]
                  - f_13 * hf1_20[k]
                  + pa_z[k] * if__50[k];

        t_91[k] = pb_y[k] * kd_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, id_30, id_57, id_59, kd_54, \
                         kd_56, kd_57, kd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * id_30[k]
                  + pb_z[k] * kd_54[k];

        t_93[k] = f_14 * id_57[k]
                  + pb_x[k] * kd_57[k];

        t_94[k] = pb_y[k] * kd_56[k];

        t_95[k] = f_14 * id_59[k]
                  + pb_x[k] * kd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, hf0_99, hf1_99, id_33, \
                         if__99, kp0_28, kp1_28, kd_57, kd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * kp0_28[k]
                  - f_2 * kp1_28[k]
                  + pb_y[k] * kd_57[k];

        t_97[k] = f_5 * id_33[k]
                  + pb_z[k] * kd_57[k];

        t_98[k] = pb_y[k] * kd_59[k];

        t_99[k] = f_15 * hf0_99[k]
                  - f_16 * hf1_99[k]
                  + pa_x[k] * if__99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, pb_z, hf0_30, hf1_30, \
                         id_36, id_63, if__60, kd_60, kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_15 * hf0_30[k]
                   - f_16 * hf1_30[k]
                   + pa_y[k] * if__60[k];

        t_101[k] = f_14 * id_36[k]
                   + pb_y[k] * kd_60[k];

        t_102[k] = pb_z[k] * kd_60[k];

        t_103[k] = f_5 * id_63[k]
                   + pb_x[k] * kd_63[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, pb_z, hf0_106, hf1_106, \
                         id_65, if__106, kd_61, kd_63, kd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * kd_61[k];

        t_105[k] = f_5 * id_65[k]
                   + pb_x[k] * kd_65[k];

        t_106[k] = f_12 * hf0_106[k]
                   - f_13 * hf1_106[k]
                   + pa_x[k] * if__106[k];

        t_107[k] = pb_z[k] * kd_63[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, id_36, id_41, \
                         if__60, if__61, kp0_32, kp1_32, kd_65, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_14 * id_41[k]
                   + pb_y[k] * kd_65[k];

        t_109[k] = f_1 * kp0_32[k]
                   - f_2 * kp1_32[k]
                   + pb_z[k] * kd_65[k];

        t_110[k] = pa_z[k] * if__60[k];

        t_111[k] = pa_z[k] * if__61[k];

        t_112[k] = f_3 * id_36[k]
                   + pb_z[k] * kd_66[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, id_39, id_70, \
                         id_71, if__63, if__66, kd_69, kd_70, kd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * if__63[k];

        t_114[k] = f_5 * id_70[k]
                   + pb_x[k] * kd_70[k];

        t_115[k] = f_5 * id_71[k]
                   + pb_x[k] * kd_71[k];

        t_116[k] = pa_z[k] * if__66[k];

        t_117[k] = f_3 * id_39[k]
                   + pb_z[k] * kd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pa_z, pb_y, hf0_50, hf1_50, id_41, \
                         id_47, id_48, if__69, if__80, kd_71, kd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * id_47[k]
                   + pb_y[k] * kd_71[k];

        t_119[k] = f_5 * id_41[k]
                   + pa_z[k] * if__69[k];

        t_120[k] = f_6 * hf0_50[k]
                   - f_7 * hf1_50[k]
                   + pa_y[k] * if__80[k];

        t_121[k] = f_8 * id_48[k]
                   + pb_y[k] * kd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pb_x, pb_z, id_42, id_75, id_76, id_77, \
                         kd_72, kd_75, kd_76, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * id_42[k]
                   + pb_z[k] * kd_72[k];

        t_123[k] = f_5 * id_75[k]
                   + pb_x[k] * kd_75[k];

        t_124[k] = f_5 * id_76[k]
                   + pb_x[k] * kd_76[k];

        t_125[k] = f_5 * id_77[k]
                   + pb_x[k] * kd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, hf0_126, hf1_126, id_45, \
                         id_53, if__126, kd_75, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_12 * hf0_126[k]
                   - f_13 * hf1_126[k]
                   + pa_x[k] * if__126[k];

        t_127[k] = f_8 * id_45[k]
                   + pb_z[k] * kd_75[k];

        t_128[k] = f_8 * id_53[k]
                   + pb_y[k] * kd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, hf0_129, hf1_129, \
                         id_54, if__90, if__92, if__129, kd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * hf0_129[k]
                   - f_13 * hf1_129[k]
                   + pa_x[k] * if__129[k];

        t_130[k] = pa_y[k] * if__90[k];

        t_131[k] = f_3 * id_54[k]
                   + pb_y[k] * kd_78[k];

        t_132[k] = pa_y[k] * if__92[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pa_y, pb_x, pb_z, id_51, id_57, \
                         id_81, id_82, if__95, if__96, kd_81, kd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * id_81[k]
                   + pb_x[k] * kd_81[k];

        t_134[k] = f_5 * id_82[k]
                   + pb_x[k] * kd_82[k];

        t_135[k] = pa_y[k] * if__95[k];

        t_136[k] = f_5 * id_57[k]
                   + pa_y[k] * if__96[k];

        t_137[k] = f_5 * id_51[k]
                   + pb_z[k] * kd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pa_z, pb_y, hf0_50, hf1_50, id_59, \
                         if__90, if__99, kd_83, kd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * id_59[k]
                   + pb_y[k] * kd_83[k];

        t_139[k] = pa_y[k] * if__99[k];

        t_140[k] = f_15 * hf0_50[k]
                   - f_16 * hf1_50[k]
                   + pa_z[k] * if__90[k];

        t_141[k] = pb_y[k] * kd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, id_54, id_87, id_89, \
                         kd_84, kd_86, kd_87, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * id_54[k]
                   + pb_z[k] * kd_84[k];

        t_143[k] = f_5 * id_87[k]
                   + pb_x[k] * kd_87[k];

        t_144[k] = pb_y[k] * kd_86[k];

        t_145[k] = f_5 * id_89[k]
                   + pb_x[k] * kd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, pb_z, hf0_149, hf1_149, \
                         id_57, if__149, kp0_43, kp1_43, kd_87, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * kp0_43[k]
                   - f_2 * kp1_43[k]
                   + pb_y[k] * kd_87[k];

        t_147[k] = f_14 * id_57[k]
                   + pb_z[k] * kd_87[k];

        t_148[k] = pb_y[k] * kd_89[k];

        t_149[k] = f_12 * hf0_149[k]
                   - f_13 * hf1_149[k]
                   + pa_x[k] * if__149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pb_x, pb_y, pb_z, hf0_60, hf1_60, \
                         id_60, id_93, if__100, kd_90, kd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_10 * hf0_60[k]
                   - f_11 * hf1_60[k]
                   + pa_y[k] * if__100[k];

        t_151[k] = f_9 * id_60[k]
                   + pb_y[k] * kd_90[k];

        t_152[k] = pb_z[k] * kd_90[k];

        t_153[k] = f_8 * id_93[k]
                   + pb_x[k] * kd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, hf0_156, hf1_156, \
                         id_95, if__156, kd_91, kd_93, kd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * kd_91[k];

        t_155[k] = f_8 * id_95[k]
                   + pb_x[k] * kd_95[k];

        t_156[k] = f_6 * hf0_156[k]
                   - f_7 * hf1_156[k]
                   + pa_x[k] * if__156[k];

        t_157[k] = pb_z[k] * kd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_y, pb_z, id_60, id_65, \
                         if__100, if__101, kp0_47, kp1_47, kd_95, \
                         kd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_9 * id_65[k]
                   + pb_y[k] * kd_95[k];

        t_159[k] = f_1 * kp0_47[k]
                   - f_2 * kp1_47[k]
                   + pb_z[k] * kd_95[k];

        t_160[k] = pa_z[k] * if__100[k];

        t_161[k] = pa_z[k] * if__101[k];

        t_162[k] = f_3 * id_60[k]
                   + pb_z[k] * kd_96[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_z, pb_x, pb_z, id_63, id_100, \
                         id_101, if__103, if__106, kd_99, kd_100, \
                         kd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * if__103[k];

        t_164[k] = f_8 * id_100[k]
                   + pb_x[k] * kd_100[k];

        t_165[k] = f_8 * id_101[k]
                   + pb_x[k] * kd_101[k];

        t_166[k] = pa_z[k] * if__106[k];

        t_167[k] = f_3 * id_63[k]
                   + pb_z[k] * kd_99[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, hf0_80, hf1_80, id_65, \
                         id_71, id_72, if__109, if__120, kd_101, \
                         kd_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_14 * id_71[k]
                   + pb_y[k] * kd_101[k];

        t_169[k] = f_5 * id_65[k]
                   + pa_z[k] * if__109[k];

        t_170[k] = f_12 * hf0_80[k]
                   - f_13 * hf1_80[k]
                   + pa_y[k] * if__120[k];

        t_171[k] = f_5 * id_72[k]
                   + pb_y[k] * kd_102[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, pb_z, id_66, id_105, id_106, \
                         id_107, kd_102, kd_105, kd_106, kd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_8 * id_66[k]
                   + pb_z[k] * kd_102[k];

        t_173[k] = f_8 * id_105[k]
                   + pb_x[k] * kd_105[k];

        t_174[k] = f_8 * id_106[k]
                   + pb_x[k] * kd_106[k];

        t_175[k] = f_8 * id_107[k]
                   + pb_x[k] * kd_107[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, pb_y, pb_z, hf0_176, hf1_176, id_69, \
                         id_77, if__176, kd_105, kd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_6 * hf0_176[k]
                   - f_7 * hf1_176[k]
                   + pa_x[k] * if__176[k];

        t_177[k] = f_8 * id_69[k]
                   + pb_z[k] * kd_105[k];

        t_178[k] = f_5 * id_77[k]
                   + pb_y[k] * kd_107[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pb_y, hf0_90, hf0_179, hf1_90, \
                         hf1_179, id_78, if__130, if__179, kd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * hf0_179[k]
                   - f_7 * hf1_179[k]
                   + pa_x[k] * if__179[k];

        t_180[k] = f_6 * hf0_90[k]
                   - f_7 * hf1_90[k]
                   + pa_y[k] * if__130[k];

        t_181[k] = f_8 * id_78[k]
                   + pb_y[k] * kd_108[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_x, pb_z, id_72, id_111, id_112, \
                         id_113, kd_108, kd_111, kd_112, kd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_5 * id_72[k]
                   + pb_z[k] * kd_108[k];

        t_183[k] = f_8 * id_111[k]
                   + pb_x[k] * kd_111[k];

        t_184[k] = f_8 * id_112[k]
                   + pb_x[k] * kd_112[k];

        t_185[k] = f_8 * id_113[k]
                   + pb_x[k] * kd_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_x, pb_y, pb_z, hf0_186, hf1_186, id_75, \
                         id_83, if__186, kd_111, kd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_6 * hf0_186[k]
                   - f_7 * hf1_186[k]
                   + pa_x[k] * if__186[k];

        t_187[k] = f_5 * id_75[k]
                   + pb_z[k] * kd_111[k];

        t_188[k] = f_8 * id_83[k]
                   + pb_y[k] * kd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pa_y, pb_y, hf0_189, hf1_189, \
                         id_84, if__140, if__142, if__189, kd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_6 * hf0_189[k]
                   - f_7 * hf1_189[k]
                   + pa_x[k] * if__189[k];

        t_190[k] = pa_y[k] * if__140[k];

        t_191[k] = f_3 * id_84[k]
                   + pb_y[k] * kd_114[k];

        t_192[k] = pa_y[k] * if__142[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_z, id_81, id_87, \
                         id_117, id_118, if__145, if__146, kd_117, \
                         kd_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * id_117[k]
                   + pb_x[k] * kd_117[k];

        t_194[k] = f_8 * id_118[k]
                   + pb_x[k] * kd_118[k];

        t_195[k] = pa_y[k] * if__145[k];

        t_196[k] = f_5 * id_87[k]
                   + pa_y[k] * if__146[k];

        t_197[k] = f_14 * id_81[k]
                   + pb_z[k] * kd_117[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_y, pa_z, pb_y, hf0_90, hf1_90, id_89, \
                         if__140, if__149, kd_119, kd_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * id_89[k]
                   + pb_y[k] * kd_119[k];

        t_199[k] = pa_y[k] * if__149[k];

        t_200[k] = f_10 * hf0_90[k]
                   - f_11 * hf1_90[k]
                   + pa_z[k] * if__140[k];

        t_201[k] = pb_y[k] * kd_120[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, id_84, id_123, id_125, \
                         kd_120, kd_122, kd_123, kd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_9 * id_84[k]
                   + pb_z[k] * kd_120[k];

        t_203[k] = f_8 * id_123[k]
                   + pb_x[k] * kd_123[k];

        t_204[k] = pb_y[k] * kd_122[k];

        t_205[k] = f_8 * id_125[k]
                   + pb_x[k] * kd_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pb_y, pb_z, hf0_209, hf1_209, \
                         id_87, if__209, kp0_61, kp1_61, kd_123, \
                         kd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * kp0_61[k]
                   - f_2 * kp1_61[k]
                   + pb_y[k] * kd_123[k];

        t_207[k] = f_9 * id_87[k]
                   + pb_z[k] * kd_123[k];

        t_208[k] = pb_y[k] * kd_125[k];

        t_209[k] = f_6 * hf0_209[k]
                   - f_7 * hf1_209[k]
                   + pa_x[k] * if__209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pa_x, pb_x, pb_y, pb_z, id_90, \
                         id_126, id_129, if__210, kd_126, kd_127, \
                         kd_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_5 * id_126[k]
                   + pa_x[k] * if__210[k];

        t_211[k] = f_4 * id_90[k]
                   + pb_y[k] * kd_126[k];

        t_212[k] = pb_z[k] * kd_126[k];

        t_213[k] = f_3 * id_129[k]
                   + pb_x[k] * kd_129[k];

        t_214[k] = pb_z[k] * kd_127[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pa_x, pb_x, pb_z, id_131, if__216, \
                         if__218, if__219, kd_129, kd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_3 * id_131[k]
                   + pb_x[k] * kd_131[k];

        t_216[k] = pa_x[k] * if__216[k];

        t_217[k] = pb_z[k] * kd_129[k];

        t_218[k] = pa_x[k] * if__218[k];

        t_219[k] = pa_x[k] * if__219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pa_z, pb_x, pb_z, id_90, id_136, \
                         if__150, if__151, if__153, kd_132, kd_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_z[k] * if__150[k];

        t_221[k] = pa_z[k] * if__151[k];

        t_222[k] = f_3 * id_90[k]
                   + pb_z[k] * kd_132[k];

        t_223[k] = pa_z[k] * if__153[k];

        t_224[k] = f_3 * id_136[k]
                   + pb_x[k] * kd_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, t_230, pa_x, pb_x, id_137, id_138, \
                         if__226, if__227, if__228, if__229, if__230, \
                         kd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_3 * id_137[k]
                   + pb_x[k] * kd_137[k];

        t_226[k] = pa_x[k] * if__226[k];

        t_227[k] = pa_x[k] * if__227[k];

        t_228[k] = pa_x[k] * if__228[k];

        t_229[k] = pa_x[k] * if__229[k];

        t_230[k] = f_5 * id_138[k]
                   + pa_x[k] * if__230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pb_y, pb_z, id_96, id_102, id_141, \
                         id_142, kd_138, kd_141, kd_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_14 * id_102[k]
                   + pb_y[k] * kd_138[k];

        t_232[k] = f_8 * id_96[k]
                   + pb_z[k] * kd_138[k];

        t_233[k] = f_3 * id_141[k]
                   + pb_x[k] * kd_141[k];

        t_234[k] = f_3 * id_142[k]
                   + pb_x[k] * kd_142[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, t_240, pa_x, pb_x, id_143, id_144, \
                         if__236, if__237, if__238, if__239, if__240, \
                         kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_3 * id_143[k]
                   + pb_x[k] * kd_143[k];

        t_236[k] = pa_x[k] * if__236[k];

        t_237[k] = pa_x[k] * if__237[k];

        t_238[k] = pa_x[k] * if__238[k];

        t_239[k] = pa_x[k] * if__239[k];

        t_240[k] = f_5 * id_144[k]
                   + pa_x[k] * if__240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_x, pb_y, pb_z, id_102, id_108, id_147, \
                         id_148, kd_144, kd_147, kd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_5 * id_108[k]
                   + pb_y[k] * kd_144[k];

        t_242[k] = f_5 * id_102[k]
                   + pb_z[k] * kd_144[k];

        t_243[k] = f_3 * id_147[k]
                   + pb_x[k] * kd_147[k];

        t_244[k] = f_3 * id_148[k]
                   + pb_x[k] * kd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, t_250, pa_x, pb_x, id_149, id_150, \
                         if__246, if__247, if__248, if__249, if__250, \
                         kd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_3 * id_149[k]
                   + pb_x[k] * kd_149[k];

        t_246[k] = pa_x[k] * if__246[k];

        t_247[k] = pa_x[k] * if__247[k];

        t_248[k] = pa_x[k] * if__248[k];

        t_249[k] = pa_x[k] * if__249[k];

        t_250[k] = f_5 * id_150[k]
                   + pa_x[k] * if__250[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_x, pb_y, pb_z, id_108, id_114, id_153, \
                         id_154, kd_150, kd_153, kd_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_8 * id_114[k]
                   + pb_y[k] * kd_150[k];

        t_252[k] = f_14 * id_108[k]
                   + pb_z[k] * kd_150[k];

        t_253[k] = f_3 * id_153[k]
                   + pb_x[k] * kd_153[k];

        t_254[k] = f_3 * id_154[k]
                   + pb_x[k] * kd_154[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, pa_x, pa_y, pb_x, id_155, \
                         if__200, if__256, if__257, if__258, if__259, \
                         kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_3 * id_155[k]
                   + pb_x[k] * kd_155[k];

        t_256[k] = pa_x[k] * if__256[k];

        t_257[k] = pa_x[k] * if__257[k];

        t_258[k] = pa_x[k] * if__258[k];

        t_259[k] = pa_x[k] * if__259[k];

        t_260[k] = pa_y[k] * if__200[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, pa_y, pb_x, pb_y, id_120, id_159, \
                         id_160, if__202, if__205, kd_156, kd_159, \
                         kd_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * id_120[k]
                   + pb_y[k] * kd_156[k];

        t_262[k] = pa_y[k] * if__202[k];

        t_263[k] = f_3 * id_159[k]
                   + pb_x[k] * kd_159[k];

        t_264[k] = f_3 * id_160[k]
                   + pb_x[k] * kd_160[k];

        t_265[k] = pa_y[k] * if__205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, t_271, pa_x, pb_y, id_162, \
                         if__266, if__267, if__268, if__269, if__270, \
                         kd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_x[k] * if__266[k];

        t_267[k] = pa_x[k] * if__267[k];

        t_268[k] = pa_x[k] * if__268[k];

        t_269[k] = pa_x[k] * if__269[k];

        t_270[k] = f_5 * id_162[k]
                   + pa_x[k] * if__270[k];

        t_271[k] = pb_y[k] * kd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, id_120, id_165, id_167, \
                         kd_162, kd_164, kd_165, kd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_4 * id_120[k]
                   + pb_z[k] * kd_162[k];

        t_273[k] = f_3 * id_165[k]
                   + pb_x[k] * kd_165[k];

        t_274[k] = pb_y[k] * kd_164[k];

        t_275[k] = f_3 * id_167[k]
                   + pb_x[k] * kd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, pa_x, pb_x, pb_y, if__276, \
                         if__277, if__279, kp0_84, kp1_84, kd_167, \
                         kd_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_x[k] * if__276[k];

        t_277[k] = pa_x[k] * if__277[k];

        t_278[k] = pb_y[k] * kd_167[k];

        t_279[k] = pa_x[k] * if__279[k];

        t_280[k] = f_1 * kp0_84[k]
                   - f_2 * kp1_84[k]
                   + pb_x[k] * kd_168[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, pb_x, pb_y, pb_z, id_126, kd_168, \
                         kd_171, kd_172, kd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_0 * id_126[k]
                   + pb_y[k] * kd_168[k];

        t_282[k] = pb_z[k] * kd_168[k];

        t_283[k] = pb_x[k] * kd_171[k];

        t_284[k] = pb_x[k] * kd_172[k];

        t_285[k] = pb_x[k] * kd_173[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pb_z, id_129, id_131, kp0_85, \
                         kp0_86, kp1_85, kp1_86, kd_171, kd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * id_129[k]
                   + f_1 * kp0_85[k]
                   - f_2 * kp1_85[k]
                   + pb_y[k] * kd_171[k];

        t_287[k] = pb_z[k] * kd_171[k];

        t_288[k] = f_0 * id_131[k]
                   + pb_y[k] * kd_173[k];

        t_289[k] = f_1 * kp0_86[k]
                   - f_2 * kp1_86[k]
                   + pb_z[k] * kd_173[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, pa_z, pb_x, pb_z, id_126, \
                         if__210, if__211, kd_174, kd_177, kd_178, \
                         kd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_z[k] * if__210[k];

        t_291[k] = pa_z[k] * if__211[k];

        t_292[k] = f_3 * id_126[k]
                   + pb_z[k] * kd_174[k];

        t_293[k] = pb_x[k] * kd_177[k];

        t_294[k] = pb_x[k] * kd_178[k];

        t_295[k] = pb_x[k] * kd_179[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_z, pb_y, pb_z, id_129, id_131, id_137, \
                         if__216, if__219, kd_177, kd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = pa_z[k] * if__216[k];

        t_297[k] = f_3 * id_129[k]
                   + pb_z[k] * kd_177[k];

        t_298[k] = f_4 * id_137[k]
                   + pb_y[k] * kd_179[k];

        t_299[k] = f_5 * id_131[k]
                   + pa_z[k] * if__219[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pb_x, pb_y, pb_z, id_132, id_138, \
                         kp0_90, kp1_90, kd_180, kd_183, kd_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * kp0_90[k]
                   - f_2 * kp1_90[k]
                   + pb_x[k] * kd_180[k];

        t_301[k] = f_9 * id_138[k]
                   + pb_y[k] * kd_180[k];

        t_302[k] = f_8 * id_132[k]
                   + pb_z[k] * kd_180[k];

        t_303[k] = pb_x[k] * kd_183[k];

        t_304[k] = pb_x[k] * kd_184[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pb_x, pb_y, pb_z, hf0_156, hf1_156, \
                         id_135, id_143, if__226, kd_183, kd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_x[k] * kd_185[k];

        t_306[k] = f_6 * hf0_156[k]
                   - f_7 * hf1_156[k]
                   + pa_z[k] * if__226[k];

        t_307[k] = f_8 * id_135[k]
                   + pb_z[k] * kd_183[k];

        t_308[k] = f_9 * id_143[k]
                   + pb_y[k] * kd_185[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_y, pb_x, pb_y, pb_z, hf0_179, hf1_179, \
                         id_138, id_144, if__239, kp0_93, kp1_93, \
                         kd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_10 * hf0_179[k]
                   - f_11 * hf1_179[k]
                   + pa_y[k] * if__239[k];

        t_310[k] = f_1 * kp0_93[k]
                   - f_2 * kp1_93[k]
                   + pb_x[k] * kd_186[k];

        t_311[k] = f_14 * id_144[k]
                   + pb_y[k] * kd_186[k];

        t_312[k] = f_5 * id_138[k]
                   + pb_z[k] * kd_186[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, pa_z, pb_x, pb_z, hf0_166, \
                         hf1_166, id_141, if__236, kd_189, kd_190, \
                         kd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_x[k] * kd_189[k];

        t_314[k] = pb_x[k] * kd_190[k];

        t_315[k] = pb_x[k] * kd_191[k];

        t_316[k] = f_12 * hf0_166[k]
                   - f_13 * hf1_166[k]
                   + pa_z[k] * if__236[k];

        t_317[k] = f_5 * id_141[k]
                   + pb_z[k] * kd_189[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pa_y, pb_x, pb_y, hf0_189, hf1_189, \
                         id_149, id_150, if__249, kp0_96, kp1_96, kd_191, \
                         kd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_14 * id_149[k]
                   + pb_y[k] * kd_191[k];

        t_319[k] = f_15 * hf0_189[k]
                   - f_16 * hf1_189[k]
                   + pa_y[k] * if__249[k];

        t_320[k] = f_1 * kp0_96[k]
                   - f_2 * kp1_96[k]
                   + pb_x[k] * kd_192[k];

        t_321[k] = f_5 * id_150[k]
                   + pb_y[k] * kd_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, pa_z, pb_x, pb_z, hf0_176, \
                         hf1_176, id_144, if__246, kd_192, kd_195, kd_196, \
                         kd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * id_144[k]
                   + pb_z[k] * kd_192[k];

        t_323[k] = pb_x[k] * kd_195[k];

        t_324[k] = pb_x[k] * kd_196[k];

        t_325[k] = pb_x[k] * kd_197[k];

        t_326[k] = f_15 * hf0_176[k]
                   - f_16 * hf1_176[k]
                   + pa_z[k] * if__246[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, pa_y, pb_y, pb_z, hf0_199, hf1_199, id_147, \
                         id_155, if__259, kd_195, kd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_14 * id_147[k]
                   + pb_z[k] * kd_195[k];

        t_328[k] = f_5 * id_155[k]
                   + pb_y[k] * kd_197[k];

        t_329[k] = f_12 * hf0_199[k]
                   - f_13 * hf1_199[k]
                   + pa_y[k] * if__259[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, pb_x, pb_y, pb_z, id_150, id_156, \
                         kp0_99, kp1_99, kd_198, kd_201, kd_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_1 * kp0_99[k]
                   - f_2 * kp1_99[k]
                   + pb_x[k] * kd_198[k];

        t_331[k] = f_8 * id_156[k]
                   + pb_y[k] * kd_198[k];

        t_332[k] = f_9 * id_150[k]
                   + pb_z[k] * kd_198[k];

        t_333[k] = pb_x[k] * kd_201[k];

        t_334[k] = pb_x[k] * kd_202[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pb_x, pb_y, pb_z, hf0_186, hf1_186, \
                         id_153, id_161, if__256, kd_201, kd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pb_x[k] * kd_203[k];

        t_336[k] = f_10 * hf0_186[k]
                   - f_11 * hf1_186[k]
                   + pa_z[k] * if__256[k];

        t_337[k] = f_9 * id_153[k]
                   + pb_z[k] * kd_201[k];

        t_338[k] = f_8 * id_161[k]
                   + pb_y[k] * kd_203[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, pa_y, pb_x, pb_y, hf0_209, \
                         hf1_209, id_162, if__269, if__270, if__272, kd_204, \
                         kd_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_6 * hf0_209[k]
                   - f_7 * hf1_209[k]
                   + pa_y[k] * if__269[k];

        t_340[k] = pa_y[k] * if__270[k];

        t_341[k] = f_3 * id_162[k]
                   + pb_y[k] * kd_204[k];

        t_342[k] = pa_y[k] * if__272[k];

        t_343[k] = pb_x[k] * kd_207[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, pa_y, pb_x, pb_y, pb_z, id_159, \
                         id_165, id_167, if__276, kd_207, kd_208, \
                         kd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pb_x[k] * kd_208[k];

        t_345[k] = pb_x[k] * kd_209[k];

        t_346[k] = f_5 * id_165[k]
                   + pa_y[k] * if__276[k];

        t_347[k] = f_4 * id_159[k]
                   + pb_z[k] * kd_207[k];

        t_348[k] = f_3 * id_167[k]
                   + pb_y[k] * kd_209[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, pa_y, pb_x, pb_y, pb_z, id_162, \
                         if__279, kp0_105, kp1_105, kd_210, kd_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = pa_y[k] * if__279[k];

        t_350[k] = f_1 * kp0_105[k]
                   - f_2 * kp1_105[k]
                   + pb_x[k] * kd_210[k];

        t_351[k] = pb_y[k] * kd_210[k];

        t_352[k] = f_0 * id_162[k]
                   + pb_z[k] * kd_210[k];

        t_353[k] = pb_x[k] * kd_213[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pb_x, pb_y, pb_z, id_165, kp0_106, \
                         kp1_106, kd_213, kd_214, kd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * kd_214[k];

        t_355[k] = pb_x[k] * kd_215[k];

        t_356[k] = f_1 * kp0_106[k]
                   - f_2 * kp1_106[k]
                   + pb_y[k] * kd_213[k];

        t_357[k] = f_0 * id_165[k]
                   + pb_z[k] * kd_213[k];

        t_358[k] = pb_y[k] * kd_215[k];
    }

#pragma omp simd aligned(t_359, pb_z, id_167, kp0_107, kp1_107, \
                         kd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_0 * id_167[k]
                   + f_1 * kp0_107[k]
                   - f_2 * kp1_107[k]
                   + pb_z[k] * kd_215[k];
    }
}

}  // namespace simdt2ceri
