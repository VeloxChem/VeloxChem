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


#include "SimdTransferMD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_md(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t mp, const size_t np, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *mp_0 = buffer.data(mp + 0);
    const auto *mp_1 = buffer.data(mp + 1);
    const auto *mp_2 = buffer.data(mp + 2);
    const auto *mp_3 = buffer.data(mp + 3);
    const auto *mp_4 = buffer.data(mp + 4);
    const auto *mp_5 = buffer.data(mp + 5);
    const auto *mp_6 = buffer.data(mp + 6);
    const auto *mp_7 = buffer.data(mp + 7);
    const auto *mp_8 = buffer.data(mp + 8);
    const auto *mp_9 = buffer.data(mp + 9);
    const auto *mp_10 = buffer.data(mp + 10);
    const auto *mp_11 = buffer.data(mp + 11);
    const auto *mp_12 = buffer.data(mp + 12);
    const auto *mp_13 = buffer.data(mp + 13);
    const auto *mp_14 = buffer.data(mp + 14);
    const auto *mp_15 = buffer.data(mp + 15);
    const auto *mp_16 = buffer.data(mp + 16);
    const auto *mp_17 = buffer.data(mp + 17);
    const auto *mp_18 = buffer.data(mp + 18);
    const auto *mp_19 = buffer.data(mp + 19);
    const auto *mp_20 = buffer.data(mp + 20);
    const auto *mp_21 = buffer.data(mp + 21);
    const auto *mp_22 = buffer.data(mp + 22);
    const auto *mp_23 = buffer.data(mp + 23);
    const auto *mp_24 = buffer.data(mp + 24);
    const auto *mp_25 = buffer.data(mp + 25);
    const auto *mp_26 = buffer.data(mp + 26);
    const auto *mp_27 = buffer.data(mp + 27);
    const auto *mp_28 = buffer.data(mp + 28);
    const auto *mp_29 = buffer.data(mp + 29);
    const auto *mp_30 = buffer.data(mp + 30);
    const auto *mp_31 = buffer.data(mp + 31);
    const auto *mp_32 = buffer.data(mp + 32);
    const auto *mp_33 = buffer.data(mp + 33);
    const auto *mp_34 = buffer.data(mp + 34);
    const auto *mp_35 = buffer.data(mp + 35);
    const auto *mp_36 = buffer.data(mp + 36);
    const auto *mp_37 = buffer.data(mp + 37);
    const auto *mp_38 = buffer.data(mp + 38);
    const auto *mp_39 = buffer.data(mp + 39);
    const auto *mp_40 = buffer.data(mp + 40);
    const auto *mp_41 = buffer.data(mp + 41);
    const auto *mp_42 = buffer.data(mp + 42);
    const auto *mp_43 = buffer.data(mp + 43);
    const auto *mp_44 = buffer.data(mp + 44);
    const auto *mp_45 = buffer.data(mp + 45);
    const auto *mp_46 = buffer.data(mp + 46);
    const auto *mp_47 = buffer.data(mp + 47);
    const auto *mp_48 = buffer.data(mp + 48);
    const auto *mp_49 = buffer.data(mp + 49);
    const auto *mp_50 = buffer.data(mp + 50);
    const auto *mp_51 = buffer.data(mp + 51);
    const auto *mp_52 = buffer.data(mp + 52);
    const auto *mp_53 = buffer.data(mp + 53);
    const auto *mp_54 = buffer.data(mp + 54);
    const auto *mp_55 = buffer.data(mp + 55);
    const auto *mp_56 = buffer.data(mp + 56);
    const auto *mp_57 = buffer.data(mp + 57);
    const auto *mp_58 = buffer.data(mp + 58);
    const auto *mp_59 = buffer.data(mp + 59);
    const auto *mp_60 = buffer.data(mp + 60);
    const auto *mp_61 = buffer.data(mp + 61);
    const auto *mp_62 = buffer.data(mp + 62);
    const auto *mp_63 = buffer.data(mp + 63);
    const auto *mp_64 = buffer.data(mp + 64);
    const auto *mp_65 = buffer.data(mp + 65);
    const auto *mp_66 = buffer.data(mp + 66);
    const auto *mp_67 = buffer.data(mp + 67);
    const auto *mp_68 = buffer.data(mp + 68);
    const auto *mp_69 = buffer.data(mp + 69);
    const auto *mp_70 = buffer.data(mp + 70);
    const auto *mp_71 = buffer.data(mp + 71);
    const auto *mp_72 = buffer.data(mp + 72);
    const auto *mp_73 = buffer.data(mp + 73);
    const auto *mp_74 = buffer.data(mp + 74);
    const auto *mp_75 = buffer.data(mp + 75);
    const auto *mp_76 = buffer.data(mp + 76);
    const auto *mp_77 = buffer.data(mp + 77);
    const auto *mp_78 = buffer.data(mp + 78);
    const auto *mp_79 = buffer.data(mp + 79);
    const auto *mp_80 = buffer.data(mp + 80);
    const auto *mp_81 = buffer.data(mp + 81);
    const auto *mp_82 = buffer.data(mp + 82);
    const auto *mp_83 = buffer.data(mp + 83);
    const auto *mp_84 = buffer.data(mp + 84);
    const auto *mp_85 = buffer.data(mp + 85);
    const auto *mp_86 = buffer.data(mp + 86);
    const auto *mp_87 = buffer.data(mp + 87);
    const auto *mp_88 = buffer.data(mp + 88);
    const auto *mp_89 = buffer.data(mp + 89);
    const auto *mp_90 = buffer.data(mp + 90);
    const auto *mp_91 = buffer.data(mp + 91);
    const auto *mp_92 = buffer.data(mp + 92);
    const auto *mp_93 = buffer.data(mp + 93);
    const auto *mp_94 = buffer.data(mp + 94);
    const auto *mp_95 = buffer.data(mp + 95);
    const auto *mp_96 = buffer.data(mp + 96);
    const auto *mp_97 = buffer.data(mp + 97);
    const auto *mp_98 = buffer.data(mp + 98);
    const auto *mp_99 = buffer.data(mp + 99);
    const auto *mp_100 = buffer.data(mp + 100);
    const auto *mp_101 = buffer.data(mp + 101);
    const auto *mp_102 = buffer.data(mp + 102);
    const auto *mp_103 = buffer.data(mp + 103);
    const auto *mp_104 = buffer.data(mp + 104);
    const auto *mp_105 = buffer.data(mp + 105);
    const auto *mp_106 = buffer.data(mp + 106);
    const auto *mp_107 = buffer.data(mp + 107);
    const auto *mp_108 = buffer.data(mp + 108);
    const auto *mp_109 = buffer.data(mp + 109);
    const auto *mp_110 = buffer.data(mp + 110);
    const auto *mp_111 = buffer.data(mp + 111);
    const auto *mp_112 = buffer.data(mp + 112);
    const auto *mp_113 = buffer.data(mp + 113);
    const auto *mp_114 = buffer.data(mp + 114);
    const auto *mp_115 = buffer.data(mp + 115);
    const auto *mp_116 = buffer.data(mp + 116);
    const auto *mp_117 = buffer.data(mp + 117);
    const auto *mp_118 = buffer.data(mp + 118);
    const auto *mp_119 = buffer.data(mp + 119);
    const auto *mp_120 = buffer.data(mp + 120);
    const auto *mp_121 = buffer.data(mp + 121);
    const auto *mp_122 = buffer.data(mp + 122);
    const auto *mp_123 = buffer.data(mp + 123);
    const auto *mp_124 = buffer.data(mp + 124);
    const auto *mp_125 = buffer.data(mp + 125);
    const auto *mp_126 = buffer.data(mp + 126);
    const auto *mp_127 = buffer.data(mp + 127);
    const auto *mp_128 = buffer.data(mp + 128);
    const auto *mp_129 = buffer.data(mp + 129);
    const auto *mp_130 = buffer.data(mp + 130);
    const auto *mp_131 = buffer.data(mp + 131);
    const auto *mp_132 = buffer.data(mp + 132);
    const auto *mp_133 = buffer.data(mp + 133);
    const auto *mp_134 = buffer.data(mp + 134);
    const auto *mp_135 = buffer.data(mp + 135);
    const auto *mp_136 = buffer.data(mp + 136);
    const auto *mp_137 = buffer.data(mp + 137);
    const auto *mp_138 = buffer.data(mp + 138);
    const auto *mp_139 = buffer.data(mp + 139);
    const auto *mp_140 = buffer.data(mp + 140);
    const auto *mp_141 = buffer.data(mp + 141);
    const auto *mp_142 = buffer.data(mp + 142);
    const auto *mp_143 = buffer.data(mp + 143);
    const auto *mp_144 = buffer.data(mp + 144);
    const auto *mp_145 = buffer.data(mp + 145);
    const auto *mp_146 = buffer.data(mp + 146);
    const auto *mp_147 = buffer.data(mp + 147);
    const auto *mp_148 = buffer.data(mp + 148);
    const auto *mp_149 = buffer.data(mp + 149);
    const auto *mp_150 = buffer.data(mp + 150);
    const auto *mp_151 = buffer.data(mp + 151);
    const auto *mp_152 = buffer.data(mp + 152);
    const auto *mp_153 = buffer.data(mp + 153);
    const auto *mp_154 = buffer.data(mp + 154);
    const auto *mp_155 = buffer.data(mp + 155);
    const auto *mp_156 = buffer.data(mp + 156);
    const auto *mp_157 = buffer.data(mp + 157);
    const auto *mp_158 = buffer.data(mp + 158);
    const auto *mp_159 = buffer.data(mp + 159);
    const auto *mp_160 = buffer.data(mp + 160);
    const auto *mp_161 = buffer.data(mp + 161);
    const auto *mp_162 = buffer.data(mp + 162);
    const auto *mp_163 = buffer.data(mp + 163);
    const auto *mp_164 = buffer.data(mp + 164);

    const auto *np_0 = buffer.data(np + 0);
    const auto *np_1 = buffer.data(np + 1);
    const auto *np_2 = buffer.data(np + 2);
    const auto *np_3 = buffer.data(np + 3);
    const auto *np_4 = buffer.data(np + 4);
    const auto *np_5 = buffer.data(np + 5);
    const auto *np_6 = buffer.data(np + 6);
    const auto *np_7 = buffer.data(np + 7);
    const auto *np_8 = buffer.data(np + 8);
    const auto *np_9 = buffer.data(np + 9);
    const auto *np_10 = buffer.data(np + 10);
    const auto *np_11 = buffer.data(np + 11);
    const auto *np_12 = buffer.data(np + 12);
    const auto *np_13 = buffer.data(np + 13);
    const auto *np_14 = buffer.data(np + 14);
    const auto *np_15 = buffer.data(np + 15);
    const auto *np_16 = buffer.data(np + 16);
    const auto *np_17 = buffer.data(np + 17);
    const auto *np_18 = buffer.data(np + 18);
    const auto *np_19 = buffer.data(np + 19);
    const auto *np_20 = buffer.data(np + 20);
    const auto *np_21 = buffer.data(np + 21);
    const auto *np_22 = buffer.data(np + 22);
    const auto *np_23 = buffer.data(np + 23);
    const auto *np_24 = buffer.data(np + 24);
    const auto *np_25 = buffer.data(np + 25);
    const auto *np_26 = buffer.data(np + 26);
    const auto *np_27 = buffer.data(np + 27);
    const auto *np_28 = buffer.data(np + 28);
    const auto *np_29 = buffer.data(np + 29);
    const auto *np_30 = buffer.data(np + 30);
    const auto *np_31 = buffer.data(np + 31);
    const auto *np_32 = buffer.data(np + 32);
    const auto *np_33 = buffer.data(np + 33);
    const auto *np_34 = buffer.data(np + 34);
    const auto *np_35 = buffer.data(np + 35);
    const auto *np_36 = buffer.data(np + 36);
    const auto *np_37 = buffer.data(np + 37);
    const auto *np_38 = buffer.data(np + 38);
    const auto *np_39 = buffer.data(np + 39);
    const auto *np_40 = buffer.data(np + 40);
    const auto *np_41 = buffer.data(np + 41);
    const auto *np_42 = buffer.data(np + 42);
    const auto *np_43 = buffer.data(np + 43);
    const auto *np_44 = buffer.data(np + 44);
    const auto *np_45 = buffer.data(np + 45);
    const auto *np_46 = buffer.data(np + 46);
    const auto *np_47 = buffer.data(np + 47);
    const auto *np_48 = buffer.data(np + 48);
    const auto *np_49 = buffer.data(np + 49);
    const auto *np_50 = buffer.data(np + 50);
    const auto *np_51 = buffer.data(np + 51);
    const auto *np_52 = buffer.data(np + 52);
    const auto *np_53 = buffer.data(np + 53);
    const auto *np_54 = buffer.data(np + 54);
    const auto *np_55 = buffer.data(np + 55);
    const auto *np_56 = buffer.data(np + 56);
    const auto *np_57 = buffer.data(np + 57);
    const auto *np_58 = buffer.data(np + 58);
    const auto *np_59 = buffer.data(np + 59);
    const auto *np_60 = buffer.data(np + 60);
    const auto *np_61 = buffer.data(np + 61);
    const auto *np_62 = buffer.data(np + 62);
    const auto *np_63 = buffer.data(np + 63);
    const auto *np_64 = buffer.data(np + 64);
    const auto *np_65 = buffer.data(np + 65);
    const auto *np_66 = buffer.data(np + 66);
    const auto *np_67 = buffer.data(np + 67);
    const auto *np_68 = buffer.data(np + 68);
    const auto *np_69 = buffer.data(np + 69);
    const auto *np_70 = buffer.data(np + 70);
    const auto *np_71 = buffer.data(np + 71);
    const auto *np_72 = buffer.data(np + 72);
    const auto *np_73 = buffer.data(np + 73);
    const auto *np_74 = buffer.data(np + 74);
    const auto *np_75 = buffer.data(np + 75);
    const auto *np_76 = buffer.data(np + 76);
    const auto *np_77 = buffer.data(np + 77);
    const auto *np_78 = buffer.data(np + 78);
    const auto *np_79 = buffer.data(np + 79);
    const auto *np_80 = buffer.data(np + 80);
    const auto *np_81 = buffer.data(np + 81);
    const auto *np_82 = buffer.data(np + 82);
    const auto *np_83 = buffer.data(np + 83);
    const auto *np_84 = buffer.data(np + 84);
    const auto *np_85 = buffer.data(np + 85);
    const auto *np_86 = buffer.data(np + 86);
    const auto *np_87 = buffer.data(np + 87);
    const auto *np_88 = buffer.data(np + 88);
    const auto *np_89 = buffer.data(np + 89);
    const auto *np_90 = buffer.data(np + 90);
    const auto *np_91 = buffer.data(np + 91);
    const auto *np_92 = buffer.data(np + 92);
    const auto *np_93 = buffer.data(np + 93);
    const auto *np_94 = buffer.data(np + 94);
    const auto *np_95 = buffer.data(np + 95);
    const auto *np_96 = buffer.data(np + 96);
    const auto *np_97 = buffer.data(np + 97);
    const auto *np_98 = buffer.data(np + 98);
    const auto *np_99 = buffer.data(np + 99);
    const auto *np_100 = buffer.data(np + 100);
    const auto *np_101 = buffer.data(np + 101);
    const auto *np_102 = buffer.data(np + 102);
    const auto *np_103 = buffer.data(np + 103);
    const auto *np_104 = buffer.data(np + 104);
    const auto *np_105 = buffer.data(np + 105);
    const auto *np_106 = buffer.data(np + 106);
    const auto *np_107 = buffer.data(np + 107);
    const auto *np_108 = buffer.data(np + 108);
    const auto *np_109 = buffer.data(np + 109);
    const auto *np_110 = buffer.data(np + 110);
    const auto *np_111 = buffer.data(np + 111);
    const auto *np_112 = buffer.data(np + 112);
    const auto *np_113 = buffer.data(np + 113);
    const auto *np_114 = buffer.data(np + 114);
    const auto *np_115 = buffer.data(np + 115);
    const auto *np_116 = buffer.data(np + 116);
    const auto *np_117 = buffer.data(np + 117);
    const auto *np_118 = buffer.data(np + 118);
    const auto *np_119 = buffer.data(np + 119);
    const auto *np_120 = buffer.data(np + 120);
    const auto *np_121 = buffer.data(np + 121);
    const auto *np_122 = buffer.data(np + 122);
    const auto *np_123 = buffer.data(np + 123);
    const auto *np_124 = buffer.data(np + 124);
    const auto *np_125 = buffer.data(np + 125);
    const auto *np_126 = buffer.data(np + 126);
    const auto *np_127 = buffer.data(np + 127);
    const auto *np_128 = buffer.data(np + 128);
    const auto *np_129 = buffer.data(np + 129);
    const auto *np_130 = buffer.data(np + 130);
    const auto *np_131 = buffer.data(np + 131);
    const auto *np_132 = buffer.data(np + 132);
    const auto *np_133 = buffer.data(np + 133);
    const auto *np_134 = buffer.data(np + 134);
    const auto *np_135 = buffer.data(np + 135);
    const auto *np_136 = buffer.data(np + 136);
    const auto *np_137 = buffer.data(np + 137);
    const auto *np_138 = buffer.data(np + 138);
    const auto *np_139 = buffer.data(np + 139);
    const auto *np_140 = buffer.data(np + 140);
    const auto *np_141 = buffer.data(np + 141);
    const auto *np_142 = buffer.data(np + 142);
    const auto *np_143 = buffer.data(np + 143);
    const auto *np_144 = buffer.data(np + 144);
    const auto *np_145 = buffer.data(np + 145);
    const auto *np_146 = buffer.data(np + 146);
    const auto *np_147 = buffer.data(np + 147);
    const auto *np_148 = buffer.data(np + 148);
    const auto *np_149 = buffer.data(np + 149);
    const auto *np_150 = buffer.data(np + 150);
    const auto *np_151 = buffer.data(np + 151);
    const auto *np_152 = buffer.data(np + 152);
    const auto *np_153 = buffer.data(np + 153);
    const auto *np_154 = buffer.data(np + 154);
    const auto *np_155 = buffer.data(np + 155);
    const auto *np_156 = buffer.data(np + 156);
    const auto *np_157 = buffer.data(np + 157);
    const auto *np_158 = buffer.data(np + 158);
    const auto *np_159 = buffer.data(np + 159);
    const auto *np_160 = buffer.data(np + 160);
    const auto *np_161 = buffer.data(np + 161);
    const auto *np_162 = buffer.data(np + 162);
    const auto *np_163 = buffer.data(np + 163);
    const auto *np_164 = buffer.data(np + 164);
    const auto *np_166 = buffer.data(np + 166);
    const auto *np_167 = buffer.data(np + 167);
    const auto *np_169 = buffer.data(np + 169);
    const auto *np_170 = buffer.data(np + 170);
    const auto *np_172 = buffer.data(np + 172);
    const auto *np_173 = buffer.data(np + 173);
    const auto *np_175 = buffer.data(np + 175);
    const auto *np_176 = buffer.data(np + 176);
    const auto *np_178 = buffer.data(np + 178);
    const auto *np_179 = buffer.data(np + 179);
    const auto *np_181 = buffer.data(np + 181);
    const auto *np_182 = buffer.data(np + 182);
    const auto *np_184 = buffer.data(np + 184);
    const auto *np_185 = buffer.data(np + 185);
    const auto *np_187 = buffer.data(np + 187);
    const auto *np_188 = buffer.data(np + 188);
    const auto *np_190 = buffer.data(np + 190);
    const auto *np_191 = buffer.data(np + 191);
    const auto *np_193 = buffer.data(np + 193);
    const auto *np_194 = buffer.data(np + 194);
    const auto *np_197 = buffer.data(np + 197);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, mp_0, mp_1, mp_2, np_0, np_1, \
                         np_2, np_4, np_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * mp_0[k]
                 + np_0[k];

        t_1[k] = ab_x[k] * mp_1[k]
                 + np_1[k];

        t_2[k] = ab_x[k] * mp_2[k]
                 + np_2[k];

        t_3[k] = ab_y[k] * mp_1[k]
                 + np_4[k];

        t_4[k] = ab_y[k] * mp_2[k]
                 + np_5[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, mp_2, mp_3, mp_4, mp_5, np_3, np_4, \
                         np_5, np_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_z[k] * mp_2[k]
                 + np_8[k];

        t_6[k] = ab_x[k] * mp_3[k]
                 + np_3[k];

        t_7[k] = ab_x[k] * mp_4[k]
                 + np_4[k];

        t_8[k] = ab_x[k] * mp_5[k]
                 + np_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, mp_4, mp_5, mp_6, np_6, \
                         np_10, np_11, np_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = ab_y[k] * mp_4[k]
                 + np_10[k];

        t_10[k] = ab_y[k] * mp_5[k]
                  + np_11[k];

        t_11[k] = ab_z[k] * mp_5[k]
                  + np_14[k];

        t_12[k] = ab_x[k] * mp_6[k]
                  + np_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, mp_7, mp_8, np_7, \
                         np_8, np_13, np_14, np_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = ab_x[k] * mp_7[k]
                  + np_7[k];

        t_14[k] = ab_x[k] * mp_8[k]
                  + np_8[k];

        t_15[k] = ab_y[k] * mp_7[k]
                  + np_13[k];

        t_16[k] = ab_y[k] * mp_8[k]
                  + np_14[k];

        t_17[k] = ab_z[k] * mp_8[k]
                  + np_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, mp_9, mp_10, mp_11, np_9, \
                         np_10, np_11, np_19, np_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_18[k] = ab_x[k] * mp_9[k]
                  + np_9[k];

        t_19[k] = ab_x[k] * mp_10[k]
                  + np_10[k];

        t_20[k] = ab_x[k] * mp_11[k]
                  + np_11[k];

        t_21[k] = ab_y[k] * mp_10[k]
                  + np_19[k];

        t_22[k] = ab_y[k] * mp_11[k]
                  + np_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, mp_11, mp_12, mp_13, mp_14, \
                         np_12, np_13, np_14, np_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_23[k] = ab_z[k] * mp_11[k]
                  + np_23[k];

        t_24[k] = ab_x[k] * mp_12[k]
                  + np_12[k];

        t_25[k] = ab_x[k] * mp_13[k]
                  + np_13[k];

        t_26[k] = ab_x[k] * mp_14[k]
                  + np_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, mp_13, mp_14, mp_15, np_15, \
                         np_22, np_23, np_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_y[k] * mp_13[k]
                  + np_22[k];

        t_28[k] = ab_y[k] * mp_14[k]
                  + np_23[k];

        t_29[k] = ab_z[k] * mp_14[k]
                  + np_26[k];

        t_30[k] = ab_x[k] * mp_15[k]
                  + np_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, mp_16, mp_17, np_16, \
                         np_17, np_25, np_26, np_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_31[k] = ab_x[k] * mp_16[k]
                  + np_16[k];

        t_32[k] = ab_x[k] * mp_17[k]
                  + np_17[k];

        t_33[k] = ab_y[k] * mp_16[k]
                  + np_25[k];

        t_34[k] = ab_y[k] * mp_17[k]
                  + np_26[k];

        t_35[k] = ab_z[k] * mp_17[k]
                  + np_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, mp_18, mp_19, mp_20, np_18, \
                         np_19, np_20, np_31, np_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_36[k] = ab_x[k] * mp_18[k]
                  + np_18[k];

        t_37[k] = ab_x[k] * mp_19[k]
                  + np_19[k];

        t_38[k] = ab_x[k] * mp_20[k]
                  + np_20[k];

        t_39[k] = ab_y[k] * mp_19[k]
                  + np_31[k];

        t_40[k] = ab_y[k] * mp_20[k]
                  + np_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, mp_20, mp_21, mp_22, mp_23, \
                         np_21, np_22, np_23, np_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_41[k] = ab_z[k] * mp_20[k]
                  + np_35[k];

        t_42[k] = ab_x[k] * mp_21[k]
                  + np_21[k];

        t_43[k] = ab_x[k] * mp_22[k]
                  + np_22[k];

        t_44[k] = ab_x[k] * mp_23[k]
                  + np_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, mp_22, mp_23, mp_24, np_24, \
                         np_34, np_35, np_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_y[k] * mp_22[k]
                  + np_34[k];

        t_46[k] = ab_y[k] * mp_23[k]
                  + np_35[k];

        t_47[k] = ab_z[k] * mp_23[k]
                  + np_38[k];

        t_48[k] = ab_x[k] * mp_24[k]
                  + np_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, mp_25, mp_26, np_25, \
                         np_26, np_37, np_38, np_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_x[k] * mp_25[k]
                  + np_25[k];

        t_50[k] = ab_x[k] * mp_26[k]
                  + np_26[k];

        t_51[k] = ab_y[k] * mp_25[k]
                  + np_37[k];

        t_52[k] = ab_y[k] * mp_26[k]
                  + np_38[k];

        t_53[k] = ab_z[k] * mp_26[k]
                  + np_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, mp_27, mp_28, mp_29, np_27, \
                         np_28, np_29, np_40, np_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * mp_27[k]
                  + np_27[k];

        t_55[k] = ab_x[k] * mp_28[k]
                  + np_28[k];

        t_56[k] = ab_x[k] * mp_29[k]
                  + np_29[k];

        t_57[k] = ab_y[k] * mp_28[k]
                  + np_40[k];

        t_58[k] = ab_y[k] * mp_29[k]
                  + np_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, mp_29, mp_30, mp_31, mp_32, \
                         np_30, np_31, np_32, np_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_59[k] = ab_z[k] * mp_29[k]
                  + np_44[k];

        t_60[k] = ab_x[k] * mp_30[k]
                  + np_30[k];

        t_61[k] = ab_x[k] * mp_31[k]
                  + np_31[k];

        t_62[k] = ab_x[k] * mp_32[k]
                  + np_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, mp_31, mp_32, mp_33, np_33, \
                         np_46, np_47, np_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_63[k] = ab_y[k] * mp_31[k]
                  + np_46[k];

        t_64[k] = ab_y[k] * mp_32[k]
                  + np_47[k];

        t_65[k] = ab_z[k] * mp_32[k]
                  + np_50[k];

        t_66[k] = ab_x[k] * mp_33[k]
                  + np_33[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, mp_34, mp_35, np_34, \
                         np_35, np_49, np_50, np_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_67[k] = ab_x[k] * mp_34[k]
                  + np_34[k];

        t_68[k] = ab_x[k] * mp_35[k]
                  + np_35[k];

        t_69[k] = ab_y[k] * mp_34[k]
                  + np_49[k];

        t_70[k] = ab_y[k] * mp_35[k]
                  + np_50[k];

        t_71[k] = ab_z[k] * mp_35[k]
                  + np_53[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, mp_36, mp_37, mp_38, np_36, \
                         np_37, np_38, np_52, np_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_72[k] = ab_x[k] * mp_36[k]
                  + np_36[k];

        t_73[k] = ab_x[k] * mp_37[k]
                  + np_37[k];

        t_74[k] = ab_x[k] * mp_38[k]
                  + np_38[k];

        t_75[k] = ab_y[k] * mp_37[k]
                  + np_52[k];

        t_76[k] = ab_y[k] * mp_38[k]
                  + np_53[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, mp_38, mp_39, mp_40, mp_41, \
                         np_39, np_40, np_41, np_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * mp_38[k]
                  + np_56[k];

        t_78[k] = ab_x[k] * mp_39[k]
                  + np_39[k];

        t_79[k] = ab_x[k] * mp_40[k]
                  + np_40[k];

        t_80[k] = ab_x[k] * mp_41[k]
                  + np_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, mp_40, mp_41, mp_42, np_42, \
                         np_55, np_56, np_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_81[k] = ab_y[k] * mp_40[k]
                  + np_55[k];

        t_82[k] = ab_y[k] * mp_41[k]
                  + np_56[k];

        t_83[k] = ab_z[k] * mp_41[k]
                  + np_59[k];

        t_84[k] = ab_x[k] * mp_42[k]
                  + np_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, mp_43, mp_44, np_43, \
                         np_44, np_58, np_59, np_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * mp_43[k]
                  + np_43[k];

        t_86[k] = ab_x[k] * mp_44[k]
                  + np_44[k];

        t_87[k] = ab_y[k] * mp_43[k]
                  + np_58[k];

        t_88[k] = ab_y[k] * mp_44[k]
                  + np_59[k];

        t_89[k] = ab_z[k] * mp_44[k]
                  + np_62[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, mp_45, mp_46, mp_47, np_45, \
                         np_46, np_47, np_64, np_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * mp_45[k]
                  + np_45[k];

        t_91[k] = ab_x[k] * mp_46[k]
                  + np_46[k];

        t_92[k] = ab_x[k] * mp_47[k]
                  + np_47[k];

        t_93[k] = ab_y[k] * mp_46[k]
                  + np_64[k];

        t_94[k] = ab_y[k] * mp_47[k]
                  + np_65[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, mp_47, mp_48, mp_49, mp_50, \
                         np_48, np_49, np_50, np_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_z[k] * mp_47[k]
                  + np_68[k];

        t_96[k] = ab_x[k] * mp_48[k]
                  + np_48[k];

        t_97[k] = ab_x[k] * mp_49[k]
                  + np_49[k];

        t_98[k] = ab_x[k] * mp_50[k]
                  + np_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, mp_49, mp_50, mp_51, \
                         np_51, np_67, np_68, np_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = ab_y[k] * mp_49[k]
                  + np_67[k];

        t_100[k] = ab_y[k] * mp_50[k]
                   + np_68[k];

        t_101[k] = ab_z[k] * mp_50[k]
                   + np_71[k];

        t_102[k] = ab_x[k] * mp_51[k]
                   + np_51[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, mp_52, mp_53, \
                         np_52, np_53, np_70, np_71, np_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_103[k] = ab_x[k] * mp_52[k]
                   + np_52[k];

        t_104[k] = ab_x[k] * mp_53[k]
                   + np_53[k];

        t_105[k] = ab_y[k] * mp_52[k]
                   + np_70[k];

        t_106[k] = ab_y[k] * mp_53[k]
                   + np_71[k];

        t_107[k] = ab_z[k] * mp_53[k]
                   + np_74[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, mp_54, mp_55, mp_56, \
                         np_54, np_55, np_56, np_73, np_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_108[k] = ab_x[k] * mp_54[k]
                   + np_54[k];

        t_109[k] = ab_x[k] * mp_55[k]
                   + np_55[k];

        t_110[k] = ab_x[k] * mp_56[k]
                   + np_56[k];

        t_111[k] = ab_y[k] * mp_55[k]
                   + np_73[k];

        t_112[k] = ab_y[k] * mp_56[k]
                   + np_74[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, mp_56, mp_57, mp_58, mp_59, \
                         np_57, np_58, np_59, np_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_113[k] = ab_z[k] * mp_56[k]
                   + np_77[k];

        t_114[k] = ab_x[k] * mp_57[k]
                   + np_57[k];

        t_115[k] = ab_x[k] * mp_58[k]
                   + np_58[k];

        t_116[k] = ab_x[k] * mp_59[k]
                   + np_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, mp_58, mp_59, mp_60, \
                         np_60, np_76, np_77, np_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_117[k] = ab_y[k] * mp_58[k]
                   + np_76[k];

        t_118[k] = ab_y[k] * mp_59[k]
                   + np_77[k];

        t_119[k] = ab_z[k] * mp_59[k]
                   + np_80[k];

        t_120[k] = ab_x[k] * mp_60[k]
                   + np_60[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, mp_61, mp_62, \
                         np_61, np_62, np_79, np_80, np_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_121[k] = ab_x[k] * mp_61[k]
                   + np_61[k];

        t_122[k] = ab_x[k] * mp_62[k]
                   + np_62[k];

        t_123[k] = ab_y[k] * mp_61[k]
                   + np_79[k];

        t_124[k] = ab_y[k] * mp_62[k]
                   + np_80[k];

        t_125[k] = ab_z[k] * mp_62[k]
                   + np_83[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ab_y, mp_63, mp_64, mp_65, \
                         np_63, np_64, np_65, np_85, np_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_126[k] = ab_x[k] * mp_63[k]
                   + np_63[k];

        t_127[k] = ab_x[k] * mp_64[k]
                   + np_64[k];

        t_128[k] = ab_x[k] * mp_65[k]
                   + np_65[k];

        t_129[k] = ab_y[k] * mp_64[k]
                   + np_85[k];

        t_130[k] = ab_y[k] * mp_65[k]
                   + np_86[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, ab_x, ab_z, mp_65, mp_66, mp_67, mp_68, \
                         np_66, np_67, np_68, np_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_131[k] = ab_z[k] * mp_65[k]
                   + np_89[k];

        t_132[k] = ab_x[k] * mp_66[k]
                   + np_66[k];

        t_133[k] = ab_x[k] * mp_67[k]
                   + np_67[k];

        t_134[k] = ab_x[k] * mp_68[k]
                   + np_68[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, mp_67, mp_68, mp_69, \
                         np_69, np_88, np_89, np_92 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_y[k] * mp_67[k]
                   + np_88[k];

        t_136[k] = ab_y[k] * mp_68[k]
                   + np_89[k];

        t_137[k] = ab_z[k] * mp_68[k]
                   + np_92[k];

        t_138[k] = ab_x[k] * mp_69[k]
                   + np_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, mp_70, mp_71, \
                         np_70, np_71, np_91, np_92, np_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_139[k] = ab_x[k] * mp_70[k]
                   + np_70[k];

        t_140[k] = ab_x[k] * mp_71[k]
                   + np_71[k];

        t_141[k] = ab_y[k] * mp_70[k]
                   + np_91[k];

        t_142[k] = ab_y[k] * mp_71[k]
                   + np_92[k];

        t_143[k] = ab_z[k] * mp_71[k]
                   + np_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_x, ab_y, mp_72, mp_73, mp_74, \
                         np_72, np_73, np_74, np_94, np_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_144[k] = ab_x[k] * mp_72[k]
                   + np_72[k];

        t_145[k] = ab_x[k] * mp_73[k]
                   + np_73[k];

        t_146[k] = ab_x[k] * mp_74[k]
                   + np_74[k];

        t_147[k] = ab_y[k] * mp_73[k]
                   + np_94[k];

        t_148[k] = ab_y[k] * mp_74[k]
                   + np_95[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, ab_x, ab_z, mp_74, mp_75, mp_76, mp_77, \
                         np_75, np_76, np_77, np_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_149[k] = ab_z[k] * mp_74[k]
                   + np_98[k];

        t_150[k] = ab_x[k] * mp_75[k]
                   + np_75[k];

        t_151[k] = ab_x[k] * mp_76[k]
                   + np_76[k];

        t_152[k] = ab_x[k] * mp_77[k]
                   + np_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, ab_x, ab_y, ab_z, mp_76, mp_77, mp_78, \
                         np_78, np_97, np_98, np_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_153[k] = ab_y[k] * mp_76[k]
                   + np_97[k];

        t_154[k] = ab_y[k] * mp_77[k]
                   + np_98[k];

        t_155[k] = ab_z[k] * mp_77[k]
                   + np_101[k];

        t_156[k] = ab_x[k] * mp_78[k]
                   + np_78[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, mp_79, mp_80, \
                         np_79, np_80, np_100, np_101, np_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_157[k] = ab_x[k] * mp_79[k]
                   + np_79[k];

        t_158[k] = ab_x[k] * mp_80[k]
                   + np_80[k];

        t_159[k] = ab_y[k] * mp_79[k]
                   + np_100[k];

        t_160[k] = ab_y[k] * mp_80[k]
                   + np_101[k];

        t_161[k] = ab_z[k] * mp_80[k]
                   + np_104[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, mp_81, mp_82, mp_83, \
                         np_81, np_82, np_83, np_103, np_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_162[k] = ab_x[k] * mp_81[k]
                   + np_81[k];

        t_163[k] = ab_x[k] * mp_82[k]
                   + np_82[k];

        t_164[k] = ab_x[k] * mp_83[k]
                   + np_83[k];

        t_165[k] = ab_y[k] * mp_82[k]
                   + np_103[k];

        t_166[k] = ab_y[k] * mp_83[k]
                   + np_104[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, ab_x, ab_z, mp_83, mp_84, mp_85, mp_86, \
                         np_84, np_85, np_86, np_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_167[k] = ab_z[k] * mp_83[k]
                   + np_107[k];

        t_168[k] = ab_x[k] * mp_84[k]
                   + np_84[k];

        t_169[k] = ab_x[k] * mp_85[k]
                   + np_85[k];

        t_170[k] = ab_x[k] * mp_86[k]
                   + np_86[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, ab_x, ab_y, ab_z, mp_85, mp_86, mp_87, \
                         np_87, np_109, np_110, np_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_171[k] = ab_y[k] * mp_85[k]
                   + np_109[k];

        t_172[k] = ab_y[k] * mp_86[k]
                   + np_110[k];

        t_173[k] = ab_z[k] * mp_86[k]
                   + np_113[k];

        t_174[k] = ab_x[k] * mp_87[k]
                   + np_87[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, mp_88, mp_89, \
                         np_88, np_89, np_112, np_113, np_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = ab_x[k] * mp_88[k]
                   + np_88[k];

        t_176[k] = ab_x[k] * mp_89[k]
                   + np_89[k];

        t_177[k] = ab_y[k] * mp_88[k]
                   + np_112[k];

        t_178[k] = ab_y[k] * mp_89[k]
                   + np_113[k];

        t_179[k] = ab_z[k] * mp_89[k]
                   + np_116[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, mp_90, mp_91, mp_92, \
                         np_90, np_91, np_92, np_115, np_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = ab_x[k] * mp_90[k]
                   + np_90[k];

        t_181[k] = ab_x[k] * mp_91[k]
                   + np_91[k];

        t_182[k] = ab_x[k] * mp_92[k]
                   + np_92[k];

        t_183[k] = ab_y[k] * mp_91[k]
                   + np_115[k];

        t_184[k] = ab_y[k] * mp_92[k]
                   + np_116[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, ab_x, ab_z, mp_92, mp_93, mp_94, mp_95, \
                         np_93, np_94, np_95, np_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = ab_z[k] * mp_92[k]
                   + np_119[k];

        t_186[k] = ab_x[k] * mp_93[k]
                   + np_93[k];

        t_187[k] = ab_x[k] * mp_94[k]
                   + np_94[k];

        t_188[k] = ab_x[k] * mp_95[k]
                   + np_95[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, ab_x, ab_y, ab_z, mp_94, mp_95, mp_96, \
                         np_96, np_118, np_119, np_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_189[k] = ab_y[k] * mp_94[k]
                   + np_118[k];

        t_190[k] = ab_y[k] * mp_95[k]
                   + np_119[k];

        t_191[k] = ab_z[k] * mp_95[k]
                   + np_122[k];

        t_192[k] = ab_x[k] * mp_96[k]
                   + np_96[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_x, ab_y, ab_z, mp_97, mp_98, \
                         np_97, np_98, np_121, np_122, np_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_193[k] = ab_x[k] * mp_97[k]
                   + np_97[k];

        t_194[k] = ab_x[k] * mp_98[k]
                   + np_98[k];

        t_195[k] = ab_y[k] * mp_97[k]
                   + np_121[k];

        t_196[k] = ab_y[k] * mp_98[k]
                   + np_122[k];

        t_197[k] = ab_z[k] * mp_98[k]
                   + np_125[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, ab_x, ab_y, mp_99, mp_100, mp_101, \
                         np_99, np_100, np_101, np_124, np_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_198[k] = ab_x[k] * mp_99[k]
                   + np_99[k];

        t_199[k] = ab_x[k] * mp_100[k]
                   + np_100[k];

        t_200[k] = ab_x[k] * mp_101[k]
                   + np_101[k];

        t_201[k] = ab_y[k] * mp_100[k]
                   + np_124[k];

        t_202[k] = ab_y[k] * mp_101[k]
                   + np_125[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, ab_x, ab_z, mp_101, mp_102, mp_103, \
                         mp_104, np_102, np_103, np_104, np_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_203[k] = ab_z[k] * mp_101[k]
                   + np_128[k];

        t_204[k] = ab_x[k] * mp_102[k]
                   + np_102[k];

        t_205[k] = ab_x[k] * mp_103[k]
                   + np_103[k];

        t_206[k] = ab_x[k] * mp_104[k]
                   + np_104[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, ab_x, ab_y, ab_z, mp_103, mp_104, mp_105, \
                         np_105, np_127, np_128, np_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_207[k] = ab_y[k] * mp_103[k]
                   + np_127[k];

        t_208[k] = ab_y[k] * mp_104[k]
                   + np_128[k];

        t_209[k] = ab_z[k] * mp_104[k]
                   + np_131[k];

        t_210[k] = ab_x[k] * mp_105[k]
                   + np_105[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ab_x, ab_y, ab_z, mp_106, mp_107, \
                         np_106, np_107, np_130, np_131, np_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_211[k] = ab_x[k] * mp_106[k]
                   + np_106[k];

        t_212[k] = ab_x[k] * mp_107[k]
                   + np_107[k];

        t_213[k] = ab_y[k] * mp_106[k]
                   + np_130[k];

        t_214[k] = ab_y[k] * mp_107[k]
                   + np_131[k];

        t_215[k] = ab_z[k] * mp_107[k]
                   + np_134[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, ab_x, ab_y, mp_108, mp_109, \
                         mp_110, np_108, np_109, np_110, np_136, \
                         np_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_216[k] = ab_x[k] * mp_108[k]
                   + np_108[k];

        t_217[k] = ab_x[k] * mp_109[k]
                   + np_109[k];

        t_218[k] = ab_x[k] * mp_110[k]
                   + np_110[k];

        t_219[k] = ab_y[k] * mp_109[k]
                   + np_136[k];

        t_220[k] = ab_y[k] * mp_110[k]
                   + np_137[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, ab_x, ab_z, mp_110, mp_111, mp_112, \
                         mp_113, np_111, np_112, np_113, np_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_221[k] = ab_z[k] * mp_110[k]
                   + np_140[k];

        t_222[k] = ab_x[k] * mp_111[k]
                   + np_111[k];

        t_223[k] = ab_x[k] * mp_112[k]
                   + np_112[k];

        t_224[k] = ab_x[k] * mp_113[k]
                   + np_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, ab_x, ab_y, ab_z, mp_112, mp_113, mp_114, \
                         np_114, np_139, np_140, np_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = ab_y[k] * mp_112[k]
                   + np_139[k];

        t_226[k] = ab_y[k] * mp_113[k]
                   + np_140[k];

        t_227[k] = ab_z[k] * mp_113[k]
                   + np_143[k];

        t_228[k] = ab_x[k] * mp_114[k]
                   + np_114[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, ab_x, ab_y, ab_z, mp_115, mp_116, \
                         np_115, np_116, np_142, np_143, np_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_229[k] = ab_x[k] * mp_115[k]
                   + np_115[k];

        t_230[k] = ab_x[k] * mp_116[k]
                   + np_116[k];

        t_231[k] = ab_y[k] * mp_115[k]
                   + np_142[k];

        t_232[k] = ab_y[k] * mp_116[k]
                   + np_143[k];

        t_233[k] = ab_z[k] * mp_116[k]
                   + np_146[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, ab_y, mp_117, mp_118, \
                         mp_119, np_117, np_118, np_119, np_145, \
                         np_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_234[k] = ab_x[k] * mp_117[k]
                   + np_117[k];

        t_235[k] = ab_x[k] * mp_118[k]
                   + np_118[k];

        t_236[k] = ab_x[k] * mp_119[k]
                   + np_119[k];

        t_237[k] = ab_y[k] * mp_118[k]
                   + np_145[k];

        t_238[k] = ab_y[k] * mp_119[k]
                   + np_146[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, ab_x, ab_z, mp_119, mp_120, mp_121, \
                         mp_122, np_120, np_121, np_122, np_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_239[k] = ab_z[k] * mp_119[k]
                   + np_149[k];

        t_240[k] = ab_x[k] * mp_120[k]
                   + np_120[k];

        t_241[k] = ab_x[k] * mp_121[k]
                   + np_121[k];

        t_242[k] = ab_x[k] * mp_122[k]
                   + np_122[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, ab_x, ab_y, ab_z, mp_121, mp_122, mp_123, \
                         np_123, np_148, np_149, np_152 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_243[k] = ab_y[k] * mp_121[k]
                   + np_148[k];

        t_244[k] = ab_y[k] * mp_122[k]
                   + np_149[k];

        t_245[k] = ab_z[k] * mp_122[k]
                   + np_152[k];

        t_246[k] = ab_x[k] * mp_123[k]
                   + np_123[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, ab_x, ab_y, ab_z, mp_124, mp_125, \
                         np_124, np_125, np_151, np_152, np_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_247[k] = ab_x[k] * mp_124[k]
                   + np_124[k];

        t_248[k] = ab_x[k] * mp_125[k]
                   + np_125[k];

        t_249[k] = ab_y[k] * mp_124[k]
                   + np_151[k];

        t_250[k] = ab_y[k] * mp_125[k]
                   + np_152[k];

        t_251[k] = ab_z[k] * mp_125[k]
                   + np_155[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, ab_y, mp_126, mp_127, \
                         mp_128, np_126, np_127, np_128, np_154, \
                         np_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_252[k] = ab_x[k] * mp_126[k]
                   + np_126[k];

        t_253[k] = ab_x[k] * mp_127[k]
                   + np_127[k];

        t_254[k] = ab_x[k] * mp_128[k]
                   + np_128[k];

        t_255[k] = ab_y[k] * mp_127[k]
                   + np_154[k];

        t_256[k] = ab_y[k] * mp_128[k]
                   + np_155[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, ab_x, ab_z, mp_128, mp_129, mp_130, \
                         mp_131, np_129, np_130, np_131, np_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_257[k] = ab_z[k] * mp_128[k]
                   + np_158[k];

        t_258[k] = ab_x[k] * mp_129[k]
                   + np_129[k];

        t_259[k] = ab_x[k] * mp_130[k]
                   + np_130[k];

        t_260[k] = ab_x[k] * mp_131[k]
                   + np_131[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, ab_x, ab_y, ab_z, mp_130, mp_131, mp_132, \
                         np_132, np_157, np_158, np_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_261[k] = ab_y[k] * mp_130[k]
                   + np_157[k];

        t_262[k] = ab_y[k] * mp_131[k]
                   + np_158[k];

        t_263[k] = ab_z[k] * mp_131[k]
                   + np_161[k];

        t_264[k] = ab_x[k] * mp_132[k]
                   + np_132[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, mp_133, mp_134, \
                         np_133, np_134, np_160, np_161, np_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = ab_x[k] * mp_133[k]
                   + np_133[k];

        t_266[k] = ab_x[k] * mp_134[k]
                   + np_134[k];

        t_267[k] = ab_y[k] * mp_133[k]
                   + np_160[k];

        t_268[k] = ab_y[k] * mp_134[k]
                   + np_161[k];

        t_269[k] = ab_z[k] * mp_134[k]
                   + np_164[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, ab_y, mp_135, mp_136, \
                         mp_137, np_135, np_136, np_137, np_166, \
                         np_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = ab_x[k] * mp_135[k]
                   + np_135[k];

        t_271[k] = ab_x[k] * mp_136[k]
                   + np_136[k];

        t_272[k] = ab_x[k] * mp_137[k]
                   + np_137[k];

        t_273[k] = ab_y[k] * mp_136[k]
                   + np_166[k];

        t_274[k] = ab_y[k] * mp_137[k]
                   + np_167[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, ab_x, ab_z, mp_137, mp_138, mp_139, \
                         mp_140, np_138, np_139, np_140, np_170 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = ab_z[k] * mp_137[k]
                   + np_170[k];

        t_276[k] = ab_x[k] * mp_138[k]
                   + np_138[k];

        t_277[k] = ab_x[k] * mp_139[k]
                   + np_139[k];

        t_278[k] = ab_x[k] * mp_140[k]
                   + np_140[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, ab_x, ab_y, ab_z, mp_139, mp_140, mp_141, \
                         np_141, np_169, np_170, np_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_279[k] = ab_y[k] * mp_139[k]
                   + np_169[k];

        t_280[k] = ab_y[k] * mp_140[k]
                   + np_170[k];

        t_281[k] = ab_z[k] * mp_140[k]
                   + np_173[k];

        t_282[k] = ab_x[k] * mp_141[k]
                   + np_141[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, ab_x, ab_y, ab_z, mp_142, mp_143, \
                         np_142, np_143, np_172, np_173, np_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_283[k] = ab_x[k] * mp_142[k]
                   + np_142[k];

        t_284[k] = ab_x[k] * mp_143[k]
                   + np_143[k];

        t_285[k] = ab_y[k] * mp_142[k]
                   + np_172[k];

        t_286[k] = ab_y[k] * mp_143[k]
                   + np_173[k];

        t_287[k] = ab_z[k] * mp_143[k]
                   + np_176[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, ab_x, ab_y, mp_144, mp_145, \
                         mp_146, np_144, np_145, np_146, np_175, \
                         np_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_288[k] = ab_x[k] * mp_144[k]
                   + np_144[k];

        t_289[k] = ab_x[k] * mp_145[k]
                   + np_145[k];

        t_290[k] = ab_x[k] * mp_146[k]
                   + np_146[k];

        t_291[k] = ab_y[k] * mp_145[k]
                   + np_175[k];

        t_292[k] = ab_y[k] * mp_146[k]
                   + np_176[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, ab_x, ab_z, mp_146, mp_147, mp_148, \
                         mp_149, np_147, np_148, np_149, np_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_293[k] = ab_z[k] * mp_146[k]
                   + np_179[k];

        t_294[k] = ab_x[k] * mp_147[k]
                   + np_147[k];

        t_295[k] = ab_x[k] * mp_148[k]
                   + np_148[k];

        t_296[k] = ab_x[k] * mp_149[k]
                   + np_149[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, ab_x, ab_y, ab_z, mp_148, mp_149, mp_150, \
                         np_150, np_178, np_179, np_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_297[k] = ab_y[k] * mp_148[k]
                   + np_178[k];

        t_298[k] = ab_y[k] * mp_149[k]
                   + np_179[k];

        t_299[k] = ab_z[k] * mp_149[k]
                   + np_182[k];

        t_300[k] = ab_x[k] * mp_150[k]
                   + np_150[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, t_305, ab_x, ab_y, ab_z, mp_151, mp_152, \
                         np_151, np_152, np_181, np_182, np_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_301[k] = ab_x[k] * mp_151[k]
                   + np_151[k];

        t_302[k] = ab_x[k] * mp_152[k]
                   + np_152[k];

        t_303[k] = ab_y[k] * mp_151[k]
                   + np_181[k];

        t_304[k] = ab_y[k] * mp_152[k]
                   + np_182[k];

        t_305[k] = ab_z[k] * mp_152[k]
                   + np_185[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, ab_x, ab_y, mp_153, mp_154, \
                         mp_155, np_153, np_154, np_155, np_184, \
                         np_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_306[k] = ab_x[k] * mp_153[k]
                   + np_153[k];

        t_307[k] = ab_x[k] * mp_154[k]
                   + np_154[k];

        t_308[k] = ab_x[k] * mp_155[k]
                   + np_155[k];

        t_309[k] = ab_y[k] * mp_154[k]
                   + np_184[k];

        t_310[k] = ab_y[k] * mp_155[k]
                   + np_185[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, ab_x, ab_z, mp_155, mp_156, mp_157, \
                         mp_158, np_156, np_157, np_158, np_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_311[k] = ab_z[k] * mp_155[k]
                   + np_188[k];

        t_312[k] = ab_x[k] * mp_156[k]
                   + np_156[k];

        t_313[k] = ab_x[k] * mp_157[k]
                   + np_157[k];

        t_314[k] = ab_x[k] * mp_158[k]
                   + np_158[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, ab_x, ab_y, ab_z, mp_157, mp_158, mp_159, \
                         np_159, np_187, np_188, np_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_315[k] = ab_y[k] * mp_157[k]
                   + np_187[k];

        t_316[k] = ab_y[k] * mp_158[k]
                   + np_188[k];

        t_317[k] = ab_z[k] * mp_158[k]
                   + np_191[k];

        t_318[k] = ab_x[k] * mp_159[k]
                   + np_159[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, ab_x, ab_y, ab_z, mp_160, mp_161, \
                         np_160, np_161, np_190, np_191, np_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_319[k] = ab_x[k] * mp_160[k]
                   + np_160[k];

        t_320[k] = ab_x[k] * mp_161[k]
                   + np_161[k];

        t_321[k] = ab_y[k] * mp_160[k]
                   + np_190[k];

        t_322[k] = ab_y[k] * mp_161[k]
                   + np_191[k];

        t_323[k] = ab_z[k] * mp_161[k]
                   + np_194[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, ab_x, ab_y, mp_162, mp_163, \
                         mp_164, np_162, np_163, np_164, np_193, \
                         np_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_324[k] = ab_x[k] * mp_162[k]
                   + np_162[k];

        t_325[k] = ab_x[k] * mp_163[k]
                   + np_163[k];

        t_326[k] = ab_x[k] * mp_164[k]
                   + np_164[k];

        t_327[k] = ab_y[k] * mp_163[k]
                   + np_193[k];

        t_328[k] = ab_y[k] * mp_164[k]
                   + np_194[k];
    }

#pragma omp simd aligned(t_329, ab_z, mp_164, np_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_329[k] = ab_z[k] * mp_164[k]
                   + np_197[k];
    }
}

}  // namespace simdtrf
