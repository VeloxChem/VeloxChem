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


#include "SimdElectronRepulsionVrrRecGH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_gh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_35 = buffer.data(gf0 + 35);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_31 = buffer.data(gf1 + 31);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_35 = buffer.data(gf1 + 35);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fg_0, gf0_0, gf1_0, \
                         gg_0, gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pb_y[k] * gg_0[k];

        t_2[k] = pb_z[k] * gg_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_4[k] = pb_y[k] * gg_2[k];

        t_5[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, fg_5, gf0_1, gf0_2, \
                         gf1_1, gf1_2, gg_3, gg_4, gg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];

        t_7[k] = pb_z[k] * gg_3[k];

        t_8[k] = pb_y[k] * gg_4[k];

        t_9[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_10[k] = f_0 * fg_5[k]
                  + pb_x[k] * gg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, fg_7, fg_9, gg_5, gg_6, \
                         gg_8, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gg_5[k];

        t_12[k] = f_0 * fg_7[k]
                  + pb_x[k] * gg_8[k];

        t_13[k] = pb_y[k] * gg_6[k];

        t_14[k] = f_0 * fg_9[k]
                  + pb_x[k] * gg_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, gf0_3, gf0_4, gf0_5, gf1_3, \
                         gf1_4, gf1_5, gg_7, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * gf0_3[k]
                  - f_2 * gf1_3[k]
                  + pb_y[k] * gg_7[k];

        t_16[k] = pb_z[k] * gg_7[k];

        t_17[k] = f_5 * gf0_4[k]
                  - f_6 * gf1_4[k]
                  + pb_y[k] * gg_8[k];

        t_18[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pb_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, fg_0, fh_0, gf0_5, \
                         gf1_5, gg_10, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * gg_10[k];

        t_20[k] = f_1 * gf0_5[k]
                  - f_2 * gf1_5[k]
                  + pb_z[k] * gg_10[k];

        t_21[k] = pa_y[k] * fh_0[k];

        t_22[k] = f_7 * fg_0[k]
                  + pb_y[k] * gg_11[k];

        t_23[k] = pb_z[k] * gg_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_z, fg_1, fg_3, fh_1, fh_2, \
                         fh_3, gg_12, gg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * fg_1[k]
                  + pa_y[k] * fh_1[k];

        t_25[k] = pb_z[k] * gg_12[k];

        t_26[k] = pa_y[k] * fh_2[k];

        t_27[k] = f_9 * fg_3[k]
                  + pa_y[k] * fh_3[k];

        t_28[k] = pb_z[k] * gg_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, fg_4, fg_13, fh_4, \
                         gg_14, gg_15, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * fg_4[k]
                  + pb_y[k] * gg_14[k];

        t_30[k] = pa_y[k] * fh_4[k];

        t_31[k] = f_9 * fg_13[k]
                  + pb_x[k] * gg_16[k];

        t_32[k] = pb_z[k] * gg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, fg_5, fg_14, fg_15, \
                         fh_6, fh_7, gg_16, gg_17, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_9 * fg_14[k]
                  + pb_x[k] * gg_17[k];

        t_34[k] = f_9 * fg_15[k]
                  + pb_x[k] * gg_18[k];

        t_35[k] = pa_y[k] * fh_6[k];

        t_36[k] = f_10 * fg_5[k]
                  + pa_y[k] * fh_7[k];

        t_37[k] = pb_z[k] * gg_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, fg_7, fg_8, fg_9, \
                         fh_0, fh_8, fh_9, fh_10, gg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * fg_7[k]
                  + pa_y[k] * fh_8[k];

        t_39[k] = f_8 * fg_8[k]
                  + pa_y[k] * fh_9[k];

        t_40[k] = f_7 * fg_9[k]
                  + pb_y[k] * gg_19[k];

        t_41[k] = pa_y[k] * fh_10[k];

        t_42[k] = pa_z[k] * fh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, fg_0, fg_2, \
                         fh_1, fh_2, fh_3, gg_20, gg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * gg_20[k];

        t_44[k] = f_7 * fg_0[k]
                  + pb_z[k] * gg_20[k];

        t_45[k] = pa_z[k] * fh_1[k];

        t_46[k] = pb_y[k] * gg_21[k];

        t_47[k] = f_8 * fg_2[k]
                  + pa_z[k] * fh_2[k];

        t_48[k] = pa_z[k] * fh_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, fg_3, fg_4, fh_4, fh_5, \
                         gg_22, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * fg_3[k]
                  + pb_z[k] * gg_22[k];

        t_50[k] = pb_y[k] * gg_23[k];

        t_51[k] = f_9 * fg_4[k]
                  + pa_z[k] * fh_4[k];

        t_52[k] = pa_z[k] * fh_5[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, fg_22, fg_23, fg_25, \
                         fh_7, gg_24, gg_26, gg_27, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_9 * fg_22[k]
                  + pb_x[k] * gg_26[k];

        t_54[k] = f_9 * fg_23[k]
                  + pb_x[k] * gg_27[k];

        t_55[k] = pb_y[k] * gg_24[k];

        t_56[k] = f_9 * fg_25[k]
                  + pb_x[k] * gg_28[k];

        t_57[k] = pa_z[k] * fh_7[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_y, pb_z, fg_5, fg_6, fg_7, fh_8, \
                         fh_9, gg_25, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * fg_5[k]
                  + pb_z[k] * gg_25[k];

        t_59[k] = f_8 * fg_6[k]
                  + pa_z[k] * fh_8[k];

        t_60[k] = f_9 * fg_7[k]
                  + pa_z[k] * fh_9[k];

        t_61[k] = pb_y[k] * gg_28[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, pb_z, dh0_0, dh1_0, fg_9, \
                         fg_10, fh_10, fh_11, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_10 * fg_9[k]
                  + pa_z[k] * fh_10[k];

        t_63[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_y[k] * fh_11[k];

        t_64[k] = f_8 * fg_10[k]
                  + pb_y[k] * gg_29[k];

        t_65[k] = pb_z[k] * gg_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, fg_27, gf0_6, gf0_8, gf1_6, gf1_8, \
                         gg_30, gg_31, gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_8 * fg_27[k]
                  + f_5 * gf0_8[k]
                  - f_6 * gf1_8[k]
                  + pb_x[k] * gg_32[k];

        t_67[k] = pb_z[k] * gg_30[k];

        t_68[k] = f_3 * gf0_6[k]
                  - f_4 * gf1_6[k]
                  + pb_z[k] * gg_31[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, fg_12, fg_29, gf0_7, gf0_9, \
                         gf1_7, gf1_9, gg_32, gg_33, gg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_8 * fg_29[k]
                  + f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pb_x[k] * gg_34[k];

        t_70[k] = pb_z[k] * gg_32[k];

        t_71[k] = f_8 * fg_12[k]
                  + pb_y[k] * gg_33[k];

        t_72[k] = f_5 * gf0_7[k]
                  - f_6 * gf1_7[k]
                  + pb_z[k] * gg_33[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, fg_30, fg_31, fg_32, fg_33, \
                         gg_34, gg_35, gg_37, gg_38, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_8 * fg_30[k]
                  + pb_x[k] * gg_35[k];

        t_74[k] = pb_z[k] * gg_34[k];

        t_75[k] = f_8 * fg_31[k]
                  + pb_x[k] * gg_37[k];

        t_76[k] = f_8 * fg_32[k]
                  + pb_x[k] * gg_38[k];

        t_77[k] = f_8 * fg_33[k]
                  + pb_x[k] * gg_39[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, dh0_1, dh1_1, fh_30, gf0_9, \
                         gf0_10, gf1_9, gf1_10, gg_35, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_1[k]
                  + pa_x[k] * fh_30[k];

        t_79[k] = pb_z[k] * gg_35[k];

        t_80[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pb_z[k] * gg_36[k];

        t_81[k] = f_5 * gf0_10[k]
                  - f_6 * gf1_10[k]
                  + pb_z[k] * gg_37[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, fg_16, fh_12, \
                         fh_17, fh_18, gf0_11, gf1_11, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_8 * fg_16[k]
                  + pb_y[k] * gg_39[k];

        t_83[k] = f_1 * gf0_11[k]
                  - f_2 * gf1_11[k]
                  + pb_z[k] * gg_39[k];

        t_84[k] = pa_y[k] * fh_17[k];

        t_85[k] = pa_z[k] * fh_12[k];

        t_86[k] = pa_y[k] * fh_18[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, fg_11, fg_18, \
                         fh_13, fh_14, fh_19, gg_40, gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * fh_13[k];

        t_88[k] = f_7 * fg_18[k]
                  + pb_y[k] * gg_40[k];

        t_89[k] = pa_y[k] * fh_19[k];

        t_90[k] = pa_z[k] * fh_14[k];

        t_91[k] = f_7 * fg_11[k]
                  + pb_z[k] * gg_41[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, fg_20, fg_37, fh_15, \
                         fh_20, gg_42, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * fg_20[k]
                  + pb_y[k] * gg_42[k];

        t_93[k] = pa_y[k] * fh_20[k];

        t_94[k] = pa_z[k] * fh_15[k];

        t_95[k] = f_8 * fg_37[k]
                  + pb_x[k] * gg_44[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pa_z, pb_x, fg_38, fg_39, fh_16, fh_21, \
                         gg_45, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_8 * fg_38[k]
                  + pb_x[k] * gg_45[k];

        t_97[k] = f_8 * fg_39[k]
                  + pb_x[k] * gg_46[k];

        t_98[k] = pa_y[k] * fh_21[k];

        t_99[k] = pa_z[k] * fh_16[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_y, pb_z, fg_13, fg_23, fg_24, \
                         fg_25, fh_22, fh_23, gg_43, gg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_7 * fg_13[k]
                   + pb_z[k] * gg_43[k];

        t_101[k] = f_9 * fg_23[k]
                   + pa_y[k] * fh_22[k];

        t_102[k] = f_8 * fg_24[k]
                   + pa_y[k] * fh_23[k];

        t_103[k] = f_7 * fg_25[k]
                   + pb_y[k] * gg_47[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pa_z, pb_y, pb_z, dh0_0, dh1_0, \
                         fg_17, fh_17, fh_24, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * fh_24[k];

        t_105[k] = f_11 * dh0_0[k]
                   - f_12 * dh1_0[k]
                   + pa_z[k] * fh_17[k];

        t_106[k] = pb_y[k] * gg_48[k];

        t_107[k] = f_8 * fg_17[k]
                   + pb_z[k] * gg_48[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_x, pb_y, fg_43, gf0_12, gf0_14, gf1_12, \
                         gf1_14, gg_49, gg_50, gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_3 * gf0_12[k]
                   - f_4 * gf1_12[k]
                   + pb_y[k] * gg_49[k];

        t_109[k] = pb_y[k] * gg_50[k];

        t_110[k] = f_8 * fg_43[k]
                   + f_5 * gf0_14[k]
                   - f_6 * gf1_14[k]
                   + pb_x[k] * gg_52[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_y, pb_z, fg_19, fg_44, gf0_13, \
                         gf0_17, gf1_13, gf1_17, gg_51, gg_52, gg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * gf0_13[k]
                   - f_6 * gf1_13[k]
                   + pb_y[k] * gg_51[k];

        t_112[k] = f_8 * fg_19[k]
                   + pb_z[k] * gg_51[k];

        t_113[k] = pb_y[k] * gg_52[k];

        t_114[k] = f_8 * fg_44[k]
                   + f_3 * gf0_17[k]
                   - f_4 * gf1_17[k]
                   + pb_x[k] * gg_53[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, fg_45, fg_46, fg_47, \
                         fg_48, gg_53, gg_54, gg_55, gg_56, gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_8 * fg_45[k]
                   + pb_x[k] * gg_54[k];

        t_116[k] = f_8 * fg_46[k]
                   + pb_x[k] * gg_55[k];

        t_117[k] = f_8 * fg_47[k]
                   + pb_x[k] * gg_56[k];

        t_118[k] = pb_y[k] * gg_53[k];

        t_119[k] = f_8 * fg_48[k]
                   + pb_x[k] * gg_58[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, fg_21, gf0_15, gf0_16, \
                         gf0_17, gf1_15, gf1_16, gf1_17, gg_54, gg_56, \
                         gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * gf0_15[k]
                   - f_2 * gf1_15[k]
                   + pb_y[k] * gg_54[k];

        t_121[k] = f_8 * fg_21[k]
                   + pb_z[k] * gg_54[k];

        t_122[k] = f_5 * gf0_16[k]
                   - f_6 * gf1_16[k]
                   + pb_y[k] * gg_56[k];

        t_123[k] = f_3 * gf0_17[k]
                   - f_4 * gf1_17[k]
                   + pb_y[k] * gg_57[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_x, pb_y, pb_z, dh0_2, dh1_2, \
                         fg_26, fg_49, fh_36, fh_37, gg_58, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * gg_58[k];

        t_125[k] = f_11 * dh0_2[k]
                   - f_12 * dh1_2[k]
                   + pa_x[k] * fh_36[k];

        t_126[k] = f_10 * fg_49[k]
                   + pa_x[k] * fh_37[k];

        t_127[k] = f_9 * fg_26[k]
                   + pb_y[k] * gg_59[k];

        t_128[k] = pb_z[k] * gg_59[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pa_x, pb_z, fg_51, fg_52, fg_53, \
                         fh_39, fh_40, fh_41, gg_60, gg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_9 * fg_51[k]
                   + pa_x[k] * fh_39[k];

        t_130[k] = pb_z[k] * gg_60[k];

        t_131[k] = f_9 * fg_52[k]
                   + pa_x[k] * fh_40[k];

        t_132[k] = f_8 * fg_53[k]
                   + pa_x[k] * fh_41[k];

        t_133[k] = pb_z[k] * gg_61[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_y, pb_z, fg_28, fg_54, \
                         fg_55, fh_42, gg_62, gg_63, gg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_9 * fg_28[k]
                   + pb_y[k] * gg_62[k];

        t_135[k] = f_8 * fg_54[k]
                   + pa_x[k] * fh_42[k];

        t_136[k] = f_7 * fg_55[k]
                   + pb_x[k] * gg_64[k];

        t_137[k] = pb_z[k] * gg_63[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pa_x, pb_x, pb_z, fg_57, fg_58, \
                         fg_59, fh_43, gg_64, gg_65, gg_66, gg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_7 * fg_57[k]
                   + pb_x[k] * gg_65[k];

        t_139[k] = f_7 * fg_58[k]
                   + pb_x[k] * gg_66[k];

        t_140[k] = f_7 * fg_59[k]
                   + pb_x[k] * gg_67[k];

        t_141[k] = pa_x[k] * fh_43[k];

        t_142[k] = pb_z[k] * gg_64[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, pa_x, pa_z, fh_25, fh_26, \
                         fh_44, fh_45, fh_46, fh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = pa_x[k] * fh_44[k];

        t_144[k] = pa_x[k] * fh_45[k];

        t_145[k] = pa_x[k] * fh_46[k];

        t_146[k] = pa_x[k] * fh_47[k];

        t_147[k] = pa_z[k] * fh_25[k];

        t_148[k] = pa_z[k] * fh_26[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_x, pa_z, pb_y, pb_z, fg_26, fg_34, \
                         fg_63, fh_27, fh_48, gg_68, gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_7 * fg_26[k]
                   + pb_z[k] * gg_68[k];

        t_150[k] = pa_z[k] * fh_27[k];

        t_151[k] = f_8 * fg_34[k]
                   + pb_y[k] * gg_69[k];

        t_152[k] = f_9 * fg_63[k]
                   + pa_x[k] * fh_48[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_x, pa_z, pb_y, pb_z, fg_27, fg_36, \
                         fg_64, fh_28, fh_49, gg_70, gg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * fh_28[k];

        t_154[k] = f_7 * fg_27[k]
                   + pb_z[k] * gg_70[k];

        t_155[k] = f_8 * fg_36[k]
                   + pb_y[k] * gg_71[k];

        t_156[k] = f_8 * fg_64[k]
                   + pa_x[k] * fh_49[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_z, pb_x, fg_66, fg_67, fg_68, \
                         fg_69, fh_29, gg_72, gg_73, gg_74, gg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pa_z[k] * fh_29[k];

        t_158[k] = f_7 * fg_66[k]
                   + pb_x[k] * gg_72[k];

        t_159[k] = f_7 * fg_67[k]
                   + pb_x[k] * gg_73[k];

        t_160[k] = f_7 * fg_68[k]
                   + pb_x[k] * gg_74[k];

        t_161[k] = f_7 * fg_69[k]
                   + pb_x[k] * gg_75[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pa_y, fh_31, \
                         fh_50, fh_51, fh_52, fh_53, fh_54, fh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_x[k] * fh_50[k];

        t_163[k] = pa_x[k] * fh_51[k];

        t_164[k] = pa_x[k] * fh_52[k];

        t_165[k] = pa_x[k] * fh_53[k];

        t_166[k] = pa_x[k] * fh_54[k];

        t_167[k] = pa_x[k] * fh_55[k];

        t_168[k] = pa_y[k] * fh_31[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_x, pa_y, pb_y, fg_40, fg_41, \
                         fg_72, fh_32, fh_33, fh_56, gg_76, gg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_7 * fg_40[k]
                   + pb_y[k] * gg_76[k];

        t_170[k] = pa_y[k] * fh_32[k];

        t_171[k] = f_9 * fg_72[k]
                   + pa_x[k] * fh_56[k];

        t_172[k] = f_7 * fg_41[k]
                   + pb_y[k] * gg_77[k];

        t_173[k] = pa_y[k] * fh_33[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_x, pa_y, pb_y, pb_z, fg_35, fg_43, \
                         fg_74, fh_34, fh_57, gg_78, gg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_8 * fg_74[k]
                   + pa_x[k] * fh_57[k];

        t_175[k] = f_8 * fg_35[k]
                   + pb_z[k] * gg_78[k];

        t_176[k] = f_7 * fg_43[k]
                   + pb_y[k] * gg_79[k];

        t_177[k] = pa_y[k] * fh_34[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_y, pb_x, fg_75, fg_76, fg_77, \
                         fg_78, fh_35, gg_80, gg_81, gg_82, gg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_7 * fg_75[k]
                   + pb_x[k] * gg_80[k];

        t_179[k] = f_7 * fg_76[k]
                   + pb_x[k] * gg_81[k];

        t_180[k] = f_7 * fg_77[k]
                   + pb_x[k] * gg_82[k];

        t_181[k] = f_7 * fg_78[k]
                   + pb_x[k] * gg_83[k];

        t_182[k] = pa_y[k] * fh_35[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, pa_x, fg_80, fh_58, \
                         fh_59, fh_60, fh_61, fh_62, fh_63, fh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_x[k] * fh_58[k];

        t_184[k] = pa_x[k] * fh_59[k];

        t_185[k] = pa_x[k] * fh_60[k];

        t_186[k] = pa_x[k] * fh_61[k];

        t_187[k] = pa_x[k] * fh_62[k];

        t_188[k] = pa_x[k] * fh_63[k];

        t_189[k] = f_10 * fg_80[k]
                   + pa_x[k] * fh_64[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pa_x, pb_y, pb_z, fg_40, fg_83, \
                         fg_84, fh_66, fh_67, gg_84, gg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pb_y[k] * gg_84[k];

        t_191[k] = f_9 * fg_40[k]
                   + pb_z[k] * gg_84[k];

        t_192[k] = f_9 * fg_83[k]
                   + pa_x[k] * fh_66[k];

        t_193[k] = pb_y[k] * gg_85[k];

        t_194[k] = f_9 * fg_84[k]
                   + pa_x[k] * fh_67[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_x, pb_y, pb_z, fg_42, fg_85, fg_86, \
                         fh_68, fh_69, gg_86, gg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_8 * fg_85[k]
                   + pa_x[k] * fh_68[k];

        t_196[k] = f_9 * fg_42[k]
                   + pb_z[k] * gg_86[k];

        t_197[k] = pb_y[k] * gg_87[k];

        t_198[k] = f_8 * fg_86[k]
                   + pa_x[k] * fh_69[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pb_x, pb_y, fg_87, fg_88, fg_89, \
                         fg_91, gg_88, gg_89, gg_90, gg_91, gg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_7 * fg_87[k]
                   + pb_x[k] * gg_89[k];

        t_200[k] = f_7 * fg_88[k]
                   + pb_x[k] * gg_90[k];

        t_201[k] = f_7 * fg_89[k]
                   + pb_x[k] * gg_91[k];

        t_202[k] = pb_y[k] * gg_88[k];

        t_203[k] = f_7 * fg_91[k]
                   + pb_x[k] * gg_92[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, t_209, pa_x, pb_y, fh_70, fh_71, \
                         fh_72, fh_73, fh_74, gg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_x[k] * fh_70[k];

        t_205[k] = pa_x[k] * fh_71[k];

        t_206[k] = pa_x[k] * fh_72[k];

        t_207[k] = pa_x[k] * fh_73[k];

        t_208[k] = pb_y[k] * gg_92[k];

        t_209[k] = pa_x[k] * fh_74[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pb_x, pb_y, pb_z, fg_49, gf0_18, \
                         gf0_19, gf1_18, gf1_19, gg_93, gg_94, gg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_1 * gf0_18[k]
                   - f_2 * gf1_18[k]
                   + pb_x[k] * gg_93[k];

        t_211[k] = f_0 * fg_49[k]
                   + pb_y[k] * gg_93[k];

        t_212[k] = pb_z[k] * gg_93[k];

        t_213[k] = f_5 * gf0_19[k]
                   - f_6 * gf1_19[k]
                   + pb_x[k] * gg_95[k];

        t_214[k] = pb_z[k] * gg_94[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pb_x, pb_y, pb_z, fg_52, gf0_20, gf0_21, \
                         gf1_20, gf1_21, gg_95, gg_96, gg_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_5 * gf0_20[k]
                   - f_6 * gf1_20[k]
                   + pb_x[k] * gg_96[k];

        t_216[k] = f_3 * gf0_21[k]
                   - f_4 * gf1_21[k]
                   + pb_x[k] * gg_97[k];

        t_217[k] = pb_z[k] * gg_95[k];

        t_218[k] = f_0 * fg_52[k]
                   + pb_y[k] * gg_96[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, pb_x, gf0_23, gf1_23, \
                         gg_98, gg_99, gg_100, gg_101, gg_102, gg_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_3 * gf0_23[k]
                   - f_4 * gf1_23[k]
                   + pb_x[k] * gg_98[k];

        t_220[k] = pb_x[k] * gg_99[k];

        t_221[k] = pb_x[k] * gg_100[k];

        t_222[k] = pb_x[k] * gg_101[k];

        t_223[k] = pb_x[k] * gg_102[k];

        t_224[k] = pb_x[k] * gg_103[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pb_z, fg_55, gf0_21, gf0_22, \
                         gf1_21, gf1_22, gg_99, gg_100, gg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * fg_55[k]
                   + f_1 * gf0_21[k]
                   - f_2 * gf1_21[k]
                   + pb_y[k] * gg_99[k];

        t_226[k] = pb_z[k] * gg_99[k];

        t_227[k] = f_3 * gf0_21[k]
                   - f_4 * gf1_21[k]
                   + pb_z[k] * gg_100[k];

        t_228[k] = f_5 * gf0_22[k]
                   - f_6 * gf1_22[k]
                   + pb_z[k] * gg_101[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_z, pb_y, pb_z, fg_49, fg_59, \
                         fh_37, fh_38, gf0_23, gf1_23, gg_103, gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * fg_59[k]
                   + pb_y[k] * gg_103[k];

        t_230[k] = f_1 * gf0_23[k]
                   - f_2 * gf1_23[k]
                   + pb_z[k] * gg_103[k];

        t_231[k] = pa_z[k] * fh_37[k];

        t_232[k] = pa_z[k] * fh_38[k];

        t_233[k] = f_7 * fg_49[k]
                   + pb_z[k] * gg_104[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pa_z, pb_y, pb_z, fg_50, fg_51, \
                         fg_61, fh_39, fh_40, fh_41, gg_105, gg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * fh_39[k];

        t_235[k] = f_9 * fg_61[k]
                   + pb_y[k] * gg_105[k];

        t_236[k] = f_8 * fg_50[k]
                   + pa_z[k] * fh_40[k];

        t_237[k] = pa_z[k] * fh_41[k];

        t_238[k] = f_7 * fg_51[k]
                   + pb_z[k] * gg_106[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_z, pb_x, pb_y, fg_52, fg_63, \
                         fh_42, gg_107, gg_108, gg_109, gg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_9 * fg_63[k]
                   + pb_y[k] * gg_107[k];

        t_240[k] = f_9 * fg_52[k]
                   + pa_z[k] * fh_42[k];

        t_241[k] = pb_x[k] * gg_108[k];

        t_242[k] = pb_x[k] * gg_109[k];

        t_243[k] = pb_x[k] * gg_110[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, pa_z, pb_x, pb_z, fg_55, fg_56, \
                         fh_43, fh_44, gg_108, gg_111, gg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_x[k] * gg_111[k];

        t_245[k] = pb_x[k] * gg_112[k];

        t_246[k] = pa_z[k] * fh_43[k];

        t_247[k] = f_7 * fg_55[k]
                   + pb_z[k] * gg_108[k];

        t_248[k] = f_8 * fg_56[k]
                   + pa_z[k] * fh_44[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_z, pb_x, pb_y, fg_57, fg_59, fg_69, \
                         fh_45, fh_47, gf0_24, gf1_24, gg_112, gg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_9 * fg_57[k]
                   + pa_z[k] * fh_45[k];

        t_250[k] = f_9 * fg_69[k]
                   + pb_y[k] * gg_112[k];

        t_251[k] = f_10 * fg_59[k]
                   + pa_z[k] * fh_47[k];

        t_252[k] = f_1 * gf0_24[k]
                   - f_2 * gf1_24[k]
                   + pb_x[k] * gg_113[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pb_x, pb_y, pb_z, fg_60, fg_70, fg_71, \
                         gf0_25, gf1_25, gg_113, gg_114, gg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_8 * fg_70[k]
                   + pb_y[k] * gg_113[k];

        t_254[k] = f_8 * fg_60[k]
                   + pb_z[k] * gg_113[k];

        t_255[k] = f_5 * gf0_25[k]
                   - f_6 * gf1_25[k]
                   + pb_x[k] * gg_115[k];

        t_256[k] = f_8 * fg_71[k]
                   + pb_y[k] * gg_114[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pb_x, pb_y, pb_z, fg_62, fg_73, gf0_26, \
                         gf0_27, gf1_26, gf1_27, gg_115, gg_116, \
                         gg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_5 * gf0_26[k]
                   - f_6 * gf1_26[k]
                   + pb_x[k] * gg_116[k];

        t_258[k] = f_3 * gf0_27[k]
                   - f_4 * gf1_27[k]
                   + pb_x[k] * gg_117[k];

        t_259[k] = f_8 * fg_62[k]
                   + pb_z[k] * gg_115[k];

        t_260[k] = f_8 * fg_73[k]
                   + pb_y[k] * gg_116[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, t_266, pb_x, gf0_29, gf1_29, \
                         gg_118, gg_119, gg_120, gg_121, gg_122, \
                         gg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * gf0_29[k]
                   - f_4 * gf1_29[k]
                   + pb_x[k] * gg_118[k];

        t_262[k] = pb_x[k] * gg_119[k];

        t_263[k] = pb_x[k] * gg_120[k];

        t_264[k] = pb_x[k] * gg_121[k];

        t_265[k] = pb_x[k] * gg_122[k];

        t_266[k] = pb_x[k] * gg_123[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_z, pb_y, pb_z, dh0_1, dh1_1, fg_65, fg_77, \
                         fh_50, gf0_28, gf1_28, gg_119, gg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_11 * dh0_1[k]
                   - f_12 * dh1_1[k]
                   + pa_z[k] * fh_50[k];

        t_268[k] = f_8 * fg_65[k]
                   + pb_z[k] * gg_119[k];

        t_269[k] = f_8 * fg_77[k]
                   + f_5 * gf0_28[k]
                   - f_6 * gf1_28[k]
                   + pb_y[k] * gg_121[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pa_y, pb_y, dh0_2, dh1_2, fg_78, fg_79, \
                         fh_63, fh_64, gf0_29, gf1_29, gg_122, gg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_8 * fg_78[k]
                   + f_3 * gf0_29[k]
                   - f_4 * gf1_29[k]
                   + pb_y[k] * gg_122[k];

        t_271[k] = f_8 * fg_79[k]
                   + pb_y[k] * gg_123[k];

        t_272[k] = f_11 * dh0_2[k]
                   - f_12 * dh1_2[k]
                   + pa_y[k] * fh_63[k];

        t_273[k] = pa_y[k] * fh_64[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pa_y, pb_y, fg_80, fg_81, fg_82, \
                         fh_65, fh_66, fh_67, gg_124, gg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_7 * fg_80[k]
                   + pb_y[k] * gg_124[k];

        t_275[k] = pa_y[k] * fh_65[k];

        t_276[k] = f_8 * fg_81[k]
                   + pa_y[k] * fh_66[k];

        t_277[k] = f_7 * fg_82[k]
                   + pb_y[k] * gg_125[k];

        t_278[k] = pa_y[k] * fh_67[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pb_y, pb_z, fg_72, fg_83, fg_84, \
                         fh_68, fh_69, gg_126, gg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_9 * fg_83[k]
                   + pa_y[k] * fh_68[k];

        t_280[k] = f_9 * fg_72[k]
                   + pb_z[k] * gg_126[k];

        t_281[k] = f_7 * fg_84[k]
                   + pb_y[k] * gg_127[k];

        t_282[k] = pa_y[k] * fh_69[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, t_288, pa_y, pb_x, fg_87, fh_70, \
                         gg_128, gg_129, gg_130, gg_131, gg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = pb_x[k] * gg_128[k];

        t_284[k] = pb_x[k] * gg_129[k];

        t_285[k] = pb_x[k] * gg_130[k];

        t_286[k] = pb_x[k] * gg_131[k];

        t_287[k] = pb_x[k] * gg_132[k];

        t_288[k] = f_10 * fg_87[k]
                   + pa_y[k] * fh_70[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_y, pb_y, pb_z, fg_75, fg_89, fg_90, \
                         fg_91, fh_72, fh_73, gg_128, gg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_9 * fg_75[k]
                   + pb_z[k] * gg_128[k];

        t_290[k] = f_9 * fg_89[k]
                   + pa_y[k] * fh_72[k];

        t_291[k] = f_8 * fg_90[k]
                   + pa_y[k] * fh_73[k];

        t_292[k] = f_7 * fg_91[k]
                   + pb_y[k] * gg_132[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_y, pb_x, pb_y, pb_z, fg_80, fh_74, \
                         gf0_30, gf1_30, gg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_y[k] * fh_74[k];

        t_294[k] = f_1 * gf0_30[k]
                   - f_2 * gf1_30[k]
                   + pb_x[k] * gg_133[k];

        t_295[k] = pb_y[k] * gg_133[k];

        t_296[k] = f_0 * fg_80[k]
                   + pb_z[k] * gg_133[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pb_x, pb_y, gf0_31, gf0_32, gf0_33, \
                         gf1_31, gf1_32, gf1_33, gg_134, gg_135, gg_136, \
                         gg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * gf0_31[k]
                   - f_6 * gf1_31[k]
                   + pb_x[k] * gg_135[k];

        t_298[k] = pb_y[k] * gg_134[k];

        t_299[k] = f_5 * gf0_32[k]
                   - f_6 * gf1_32[k]
                   + pb_x[k] * gg_136[k];

        t_300[k] = f_3 * gf0_33[k]
                   - f_4 * gf1_33[k]
                   + pb_x[k] * gg_137[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, t_305, pb_x, pb_y, pb_z, fg_83, gf0_35, \
                         gf1_35, gg_135, gg_136, gg_138, gg_139, \
                         gg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_0 * fg_83[k]
                   + pb_z[k] * gg_135[k];

        t_302[k] = pb_y[k] * gg_136[k];

        t_303[k] = f_3 * gf0_35[k]
                   - f_4 * gf1_35[k]
                   + pb_x[k] * gg_138[k];

        t_304[k] = pb_x[k] * gg_139[k];

        t_305[k] = pb_x[k] * gg_140[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, pb_x, pb_y, pb_z, fg_87, gf0_33, \
                         gf1_33, gg_139, gg_141, gg_142, gg_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = pb_x[k] * gg_141[k];

        t_307[k] = pb_x[k] * gg_142[k];

        t_308[k] = pb_x[k] * gg_143[k];

        t_309[k] = f_1 * gf0_33[k]
                   - f_2 * gf1_33[k]
                   + pb_y[k] * gg_139[k];

        t_310[k] = f_0 * fg_87[k]
                   + pb_z[k] * gg_139[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_y, pb_z, fg_91, gf0_34, gf0_35, \
                         gf1_34, gf1_35, gg_141, gg_142, gg_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_5 * gf0_34[k]
                   - f_6 * gf1_34[k]
                   + pb_y[k] * gg_141[k];

        t_312[k] = f_3 * gf0_35[k]
                   - f_4 * gf1_35[k]
                   + pb_y[k] * gg_142[k];

        t_313[k] = pb_y[k] * gg_143[k];

        t_314[k] = f_0 * fg_91[k]
                   + f_1 * gf0_35[k]
                   - f_2 * gf1_35[k]
                   + pb_z[k] * gg_143[k];
    }
}

auto
compute_prim_gh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_56 = buffer.data(dh1 + 56);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_77 = buffer.data(fh + 77);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_35 = buffer.data(gf0 + 35);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_31 = buffer.data(gf1 + 31);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_35 = buffer.data(gf1 + 35);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fg_0, gf0_0, gf1_0, gg_0, \
                         gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pb_y[k] * gg_0[k];

        t_2[k] = pb_z[k] * gg_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_4[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, fg_5, gf0_1, gf0_2, gf1_1, \
                         gf1_2, gg_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = pb_y[k] * gg_4[k];

        t_7[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_8[k] = f_0 * fg_5[k]
                 + pb_x[k] * gg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, fg_9, gf0_3, gf0_4, gf1_3, gf1_4, gg_5, \
                         gg_6, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * fg_9[k]
                 + pb_x[k] * gg_8[k];

        t_10[k] = f_1 * gf0_3[k]
                  - f_2 * gf1_3[k]
                  + pb_y[k] * gg_5[k];

        t_11[k] = f_5 * gf0_4[k]
                  - f_6 * gf1_4[k]
                  + pb_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_y, pb_y, pb_z, fg_0, fh_0, gf0_5, \
                         gf1_5, gg_7, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pb_y[k] * gg_7[k];

        t_13[k] = pb_y[k] * gg_8[k];

        t_14[k] = f_1 * gf0_5[k]
                  - f_2 * gf1_5[k]
                  + pb_z[k] * gg_8[k];

        t_15[k] = pa_y[k] * fh_0[k];

        t_16[k] = f_7 * fg_0[k]
                  + pb_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_x, fg_1, fg_3, fg_11, fh_3, \
                         fh_4, fh_5, fh_7, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * fg_1[k]
                  + pa_y[k] * fh_3[k];

        t_18[k] = pa_y[k] * fh_4[k];

        t_19[k] = f_9 * fg_3[k]
                  + pa_y[k] * fh_5[k];

        t_20[k] = pa_y[k] * fh_7[k];

        t_21[k] = f_9 * fg_11[k]
                  + pb_x[k] * gg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, fg_5, fg_7, fg_8, fg_9, \
                         fh_8, fh_9, fh_10, fh_12, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * fg_5[k]
                  + pa_y[k] * fh_8[k];

        t_23[k] = f_9 * fg_7[k]
                  + pa_y[k] * fh_9[k];

        t_24[k] = f_8 * fg_8[k]
                  + pa_y[k] * fh_10[k];

        t_25[k] = f_7 * fg_9[k]
                  + pb_y[k] * gg_11[k];

        t_26[k] = pa_y[k] * fh_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_z, pb_z, fg_0, fg_2, fh_0, fh_3, \
                         fh_4, fh_5, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * fh_0[k];

        t_28[k] = f_7 * fg_0[k]
                  + pb_z[k] * gg_12[k];

        t_29[k] = pa_z[k] * fh_3[k];

        t_30[k] = f_8 * fg_2[k]
                  + pa_z[k] * fh_4[k];

        t_31[k] = pa_z[k] * fh_5[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_z, fg_4, fg_5, fg_18, fh_7, \
                         fh_8, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * fg_4[k]
                  + pa_z[k] * fh_7[k];

        t_33[k] = f_9 * fg_18[k]
                  + pb_x[k] * gg_14[k];

        t_34[k] = pa_z[k] * fh_8[k];

        t_35[k] = f_7 * fg_5[k]
                  + pb_z[k] * gg_13[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, dh0_0, dh1_0, fg_6, fg_7, fg_9, \
                         fh_9, fh_10, fh_12, fh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_8 * fg_6[k]
                  + pa_z[k] * fh_9[k];

        t_37[k] = f_9 * fg_7[k]
                  + pa_z[k] * fh_10[k];

        t_38[k] = f_10 * fg_9[k]
                  + pa_z[k] * fh_12[k];

        t_39[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_y[k] * fh_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pb_y, pb_z, fg_10, fg_20, gf0_6, gf0_8, \
                         gf1_6, gf1_8, gg_15, gg_16, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_8 * fg_10[k]
                  + pb_y[k] * gg_15[k];

        t_41[k] = pb_z[k] * gg_15[k];

        t_42[k] = f_8 * fg_20[k]
                  + f_5 * gf0_8[k]
                  - f_6 * gf1_8[k]
                  + pb_x[k] * gg_17[k];

        t_43[k] = f_3 * gf0_6[k]
                  - f_4 * gf1_6[k]
                  + pb_z[k] * gg_16[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_z, fg_21, fg_22, gf0_7, gf0_9, \
                         gf1_7, gf1_9, gg_17, gg_18, gg_19, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_8 * fg_21[k]
                  + f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pb_x[k] * gg_19[k];

        t_45[k] = pb_z[k] * gg_17[k];

        t_46[k] = f_5 * gf0_7[k]
                  - f_6 * gf1_7[k]
                  + pb_z[k] * gg_18[k];

        t_47[k] = f_8 * fg_22[k]
                  + pb_x[k] * gg_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_z, dh0_1, dh1_30, fh_27, gf0_9, \
                         gf0_10, gf1_9, gf1_10, gg_20, gg_21, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_30[k]
                  + pa_x[k] * fh_27[k];

        t_49[k] = pb_z[k] * gg_20[k];

        t_50[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pb_z[k] * gg_21[k];

        t_51[k] = f_5 * gf0_10[k]
                  - f_6 * gf1_10[k]
                  + pb_z[k] * gg_22[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, fg_12, fh_14, \
                         fh_18, fh_19, gf0_11, gf1_11, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * fg_12[k]
                  + pb_y[k] * gg_23[k];

        t_53[k] = f_1 * gf0_11[k]
                  - f_2 * gf1_11[k]
                  + pb_z[k] * gg_23[k];

        t_54[k] = pa_y[k] * fh_18[k];

        t_55[k] = pa_z[k] * fh_14[k];

        t_56[k] = pa_y[k] * fh_19[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_y, pa_z, pb_z, fg_11, fg_16, fh_15, \
                         fh_16, fh_20, fh_21, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * fh_15[k];

        t_58[k] = pa_y[k] * fh_20[k];

        t_59[k] = pa_z[k] * fh_16[k];

        t_60[k] = f_7 * fg_11[k]
                  + pb_z[k] * gg_24[k];

        t_61[k] = f_9 * fg_16[k]
                  + pa_y[k] * fh_21[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, dh0_0, dh1_0, fg_17, fg_18, \
                         fh_17, fh_22, fh_23, gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * fg_17[k]
                  + pa_y[k] * fh_22[k];

        t_63[k] = f_7 * fg_18[k]
                  + pb_y[k] * gg_25[k];

        t_64[k] = pa_y[k] * fh_23[k];

        t_65[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_z[k] * fh_17[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pb_y, pb_z, fg_13, fg_24, gf0_12, \
                         gf0_14, gf1_12, gf1_14, gg_26, gg_27, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_y[k] * gg_26[k];

        t_67[k] = f_8 * fg_13[k]
                  + pb_z[k] * gg_26[k];

        t_68[k] = f_3 * gf0_12[k]
                  - f_4 * gf1_12[k]
                  + pb_y[k] * gg_27[k];

        t_69[k] = f_8 * fg_24[k]
                  + f_5 * gf0_14[k]
                  - f_6 * gf1_14[k]
                  + pb_x[k] * gg_29[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pb_y, fg_25, fg_26, gf0_13, gf0_17, \
                         gf1_13, gf1_17, gg_28, gg_29, gg_30, gg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_13[k]
                  + pb_y[k] * gg_28[k];

        t_71[k] = pb_y[k] * gg_29[k];

        t_72[k] = f_8 * fg_25[k]
                  + f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pb_x[k] * gg_30[k];

        t_73[k] = f_8 * fg_26[k]
                  + pb_x[k] * gg_34[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_y, pb_z, fg_15, gf0_15, gf0_16, gf0_17, \
                         gf1_15, gf1_16, gf1_17, gg_31, gg_32, gg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * gf0_15[k]
                  - f_2 * gf1_15[k]
                  + pb_y[k] * gg_31[k];

        t_75[k] = f_8 * fg_15[k]
                  + pb_z[k] * gg_31[k];

        t_76[k] = f_5 * gf0_16[k]
                  - f_6 * gf1_16[k]
                  + pb_y[k] * gg_32[k];

        t_77[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pb_y[k] * gg_33[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_y, dh0_2, dh1_56, fg_19, fg_27, \
                         fh_32, fh_33, gg_34, gg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_y[k] * gg_34[k];

        t_79[k] = f_11 * dh0_2[k]
                  - f_12 * dh1_56[k]
                  + pa_x[k] * fh_32[k];

        t_80[k] = f_10 * fg_27[k]
                  + pa_x[k] * fh_33[k];

        t_81[k] = f_9 * fg_19[k]
                  + pb_y[k] * gg_35[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, fg_29, fg_30, fg_31, fg_32, fh_34, \
                         fh_35, fh_36, fh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_9 * fg_29[k]
                  + pa_x[k] * fh_34[k];

        t_83[k] = f_9 * fg_30[k]
                  + pa_x[k] * fh_35[k];

        t_84[k] = f_8 * fg_31[k]
                  + pa_x[k] * fh_36[k];

        t_85[k] = f_8 * fg_32[k]
                  + pa_x[k] * fh_37[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pa_x, pb_x, fg_33, fh_41, fh_43, \
                         fh_44, fh_45, fh_46, gg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_7 * fg_33[k]
                  + pb_x[k] * gg_36[k];

        t_87[k] = pa_x[k] * fh_41[k];

        t_88[k] = pa_x[k] * fh_43[k];

        t_89[k] = pa_x[k] * fh_44[k];

        t_90[k] = pa_x[k] * fh_45[k];

        t_91[k] = pa_x[k] * fh_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pa_x, pa_z, pb_z, fg_19, fg_38, fh_24, \
                         fh_25, fh_26, fh_47, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * fh_24[k];

        t_93[k] = f_7 * fg_19[k]
                  + pb_z[k] * gg_37[k];

        t_94[k] = pa_z[k] * fh_25[k];

        t_95[k] = f_9 * fg_38[k]
                  + pa_x[k] * fh_47[k];

        t_96[k] = pa_z[k] * fh_26[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, t_102, pa_x, fg_39, fh_48, fh_50, \
                         fh_51, fh_52, fh_53, fh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_8 * fg_39[k]
                  + pa_x[k] * fh_48[k];

        t_98[k] = pa_x[k] * fh_50[k];

        t_99[k] = pa_x[k] * fh_51[k];

        t_100[k] = pa_x[k] * fh_52[k];

        t_101[k] = pa_x[k] * fh_53[k];

        t_102[k] = pa_x[k] * fh_54[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, pa_x, pa_y, fg_43, fg_44, \
                         fh_28, fh_29, fh_30, fh_31, fh_55, fh_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pa_y[k] * fh_28[k];

        t_104[k] = pa_y[k] * fh_29[k];

        t_105[k] = f_9 * fg_43[k]
                   + pa_x[k] * fh_55[k];

        t_106[k] = pa_y[k] * fh_30[k];

        t_107[k] = f_8 * fg_44[k]
                   + pa_x[k] * fh_56[k];

        t_108[k] = pa_y[k] * fh_31[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, pa_x, fg_49, fh_57, fh_58, \
                         fh_59, fh_60, fh_61, fh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * fh_57[k];

        t_110[k] = pa_x[k] * fh_58[k];

        t_111[k] = pa_x[k] * fh_59[k];

        t_112[k] = pa_x[k] * fh_60[k];

        t_113[k] = pa_x[k] * fh_61[k];

        t_114[k] = f_10 * fg_49[k]
                   + pa_x[k] * fh_63[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_x, pb_z, fg_23, fg_51, fg_52, fg_53, \
                         fh_65, fh_66, fh_67, gg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_9 * fg_23[k]
                   + pb_z[k] * gg_38[k];

        t_116[k] = f_9 * fg_51[k]
                   + pa_x[k] * fh_65[k];

        t_117[k] = f_9 * fg_52[k]
                   + pa_x[k] * fh_66[k];

        t_118[k] = f_8 * fg_53[k]
                   + pa_x[k] * fh_67[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, pa_x, pb_x, fg_54, fg_58, \
                         fh_68, fh_72, fh_73, fh_74, fh_75, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_8 * fg_54[k]
                   + pa_x[k] * fh_68[k];

        t_120[k] = f_7 * fg_58[k]
                   + pb_x[k] * gg_39[k];

        t_121[k] = pa_x[k] * fh_72[k];

        t_122[k] = pa_x[k] * fh_73[k];

        t_123[k] = pa_x[k] * fh_74[k];

        t_124[k] = pa_x[k] * fh_75[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_x, pb_x, pb_y, fg_27, fh_77, gf0_18, \
                         gf0_19, gf1_18, gf1_19, gg_40, gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_x[k] * fh_77[k];

        t_126[k] = f_1 * gf0_18[k]
                   - f_2 * gf1_18[k]
                   + pb_x[k] * gg_40[k];

        t_127[k] = f_0 * fg_27[k]
                   + pb_y[k] * gg_40[k];

        t_128[k] = f_5 * gf0_19[k]
                   - f_6 * gf1_19[k]
                   + pb_x[k] * gg_41[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pb_x, gf0_20, gf0_21, gf0_23, gf1_20, \
                         gf1_21, gf1_23, gg_42, gg_43, gg_44, gg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_5 * gf0_20[k]
                   - f_6 * gf1_20[k]
                   + pb_x[k] * gg_42[k];

        t_130[k] = f_3 * gf0_21[k]
                   - f_4 * gf1_21[k]
                   + pb_x[k] * gg_43[k];

        t_131[k] = f_3 * gf0_23[k]
                   - f_4 * gf1_23[k]
                   + pb_x[k] * gg_44[k];

        t_132[k] = pb_x[k] * gg_45[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pb_x, pb_y, pb_z, fg_33, gf0_21, \
                         gf1_21, gg_45, gg_46, gg_47, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_x[k] * gg_47[k];

        t_134[k] = pb_x[k] * gg_48[k];

        t_135[k] = f_0 * fg_33[k]
                   + f_1 * gf0_21[k]
                   - f_2 * gf1_21[k]
                   + pb_y[k] * gg_45[k];

        t_136[k] = pb_z[k] * gg_45[k];

        t_137[k] = f_3 * gf0_21[k]
                   - f_4 * gf1_21[k]
                   + pb_z[k] * gg_46[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_z, pb_y, pb_z, fg_36, fh_33, gf0_22, \
                         gf0_23, gf1_22, gf1_23, gg_47, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * gf0_22[k]
                   - f_6 * gf1_22[k]
                   + pb_z[k] * gg_47[k];

        t_139[k] = f_0 * fg_36[k]
                   + pb_y[k] * gg_48[k];

        t_140[k] = f_1 * gf0_23[k]
                   - f_2 * gf1_23[k]
                   + pb_z[k] * gg_48[k];

        t_141[k] = pa_z[k] * fh_33[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, pa_z, pb_z, fg_27, fg_28, fg_30, \
                         fh_34, fh_35, fh_36, fh_37, gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_7 * fg_27[k]
                   + pb_z[k] * gg_49[k];

        t_143[k] = pa_z[k] * fh_34[k];

        t_144[k] = f_8 * fg_28[k]
                   + pa_z[k] * fh_35[k];

        t_145[k] = pa_z[k] * fh_36[k];

        t_146[k] = f_9 * fg_30[k]
                   + pa_z[k] * fh_37[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pb_z, fg_33, fg_34, fg_35, fh_41, \
                         fh_43, fh_44, gg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pa_z[k] * fh_41[k];

        t_148[k] = f_7 * fg_33[k]
                   + pb_z[k] * gg_50[k];

        t_149[k] = f_8 * fg_34[k]
                   + pa_z[k] * fh_43[k];

        t_150[k] = f_9 * fg_35[k]
                   + pa_z[k] * fh_44[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_z, pb_x, pb_y, pb_z, fg_36, fg_37, \
                         fg_42, fh_46, gf0_24, gf1_24, gg_51, gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_9 * fg_42[k]
                   + pb_y[k] * gg_51[k];

        t_152[k] = f_10 * fg_36[k]
                   + pa_z[k] * fh_46[k];

        t_153[k] = f_1 * gf0_24[k]
                   - f_2 * gf1_24[k]
                   + pb_x[k] * gg_52[k];

        t_154[k] = f_8 * fg_37[k]
                   + pb_z[k] * gg_52[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pb_x, gf0_25, gf0_26, gf0_27, gf1_25, gf1_26, \
                         gf1_27, gg_53, gg_54, gg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_5 * gf0_25[k]
                   - f_6 * gf1_25[k]
                   + pb_x[k] * gg_53[k];

        t_156[k] = f_5 * gf0_26[k]
                   - f_6 * gf1_26[k]
                   + pb_x[k] * gg_54[k];

        t_157[k] = f_3 * gf0_27[k]
                   - f_4 * gf1_27[k]
                   + pb_x[k] * gg_55[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_x, dh0_1, dh1_30, fh_49, \
                         gf0_29, gf1_29, gg_56, gg_57, gg_58, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * gf0_29[k]
                   - f_4 * gf1_29[k]
                   + pb_x[k] * gg_56[k];

        t_159[k] = pb_x[k] * gg_57[k];

        t_160[k] = pb_x[k] * gg_58[k];

        t_161[k] = pb_x[k] * gg_60[k];

        t_162[k] = f_11 * dh0_1[k]
                   - f_12 * dh1_30[k]
                   + pa_z[k] * fh_49[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pb_y, pb_z, fg_40, fg_46, fg_47, gf0_28, gf0_29, \
                         gf1_28, gf1_29, gg_57, gg_58, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_8 * fg_40[k]
                   + pb_z[k] * gg_57[k];

        t_164[k] = f_8 * fg_46[k]
                   + f_5 * gf0_28[k]
                   - f_6 * gf1_28[k]
                   + pb_y[k] * gg_58[k];

        t_165[k] = f_8 * fg_47[k]
                   + f_3 * gf0_29[k]
                   - f_4 * gf1_29[k]
                   + pb_y[k] * gg_59[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pa_y, pb_y, dh0_2, dh1_56, fg_48, \
                         fg_50, fh_62, fh_63, fh_64, fh_65, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_8 * fg_48[k]
                   + pb_y[k] * gg_60[k];

        t_167[k] = f_11 * dh0_2[k]
                   - f_12 * dh1_56[k]
                   + pa_y[k] * fh_62[k];

        t_168[k] = pa_y[k] * fh_63[k];

        t_169[k] = pa_y[k] * fh_64[k];

        t_170[k] = f_8 * fg_50[k]
                   + pa_y[k] * fh_65[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, t_175, pa_y, pb_z, fg_45, fg_51, fg_55, \
                         fh_66, fh_67, fh_68, fh_72, gg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * fh_66[k];

        t_172[k] = f_9 * fg_51[k]
                   + pa_y[k] * fh_67[k];

        t_173[k] = pa_y[k] * fh_68[k];

        t_174[k] = f_10 * fg_55[k]
                   + pa_y[k] * fh_72[k];

        t_175[k] = f_9 * fg_45[k]
                   + pb_z[k] * gg_61[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, fg_56, fg_57, fg_58, fh_74, \
                         fh_75, fh_77, gg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * fg_56[k]
                   + pa_y[k] * fh_74[k];

        t_177[k] = f_8 * fg_57[k]
                   + pa_y[k] * fh_75[k];

        t_178[k] = f_7 * fg_58[k]
                   + pb_y[k] * gg_62[k];

        t_179[k] = pa_y[k] * fh_77[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pb_x, pb_z, fg_49, gf0_30, gf0_31, \
                         gf0_32, gf1_30, gf1_31, gf1_32, gg_63, gg_64, \
                         gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * gf0_30[k]
                   - f_2 * gf1_30[k]
                   + pb_x[k] * gg_63[k];

        t_181[k] = f_0 * fg_49[k]
                   + pb_z[k] * gg_63[k];

        t_182[k] = f_5 * gf0_31[k]
                   - f_6 * gf1_31[k]
                   + pb_x[k] * gg_64[k];

        t_183[k] = f_5 * gf0_32[k]
                   - f_6 * gf1_32[k]
                   + pb_x[k] * gg_65[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pb_x, gf0_33, gf0_35, gf1_33, \
                         gf1_35, gg_66, gg_67, gg_68, gg_69, gg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * gf0_33[k]
                   - f_4 * gf1_33[k]
                   + pb_x[k] * gg_66[k];

        t_185[k] = f_3 * gf0_35[k]
                   - f_4 * gf1_35[k]
                   + pb_x[k] * gg_67[k];

        t_186[k] = pb_x[k] * gg_68[k];

        t_187[k] = pb_x[k] * gg_69[k];

        t_188[k] = pb_x[k] * gg_71[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_y, pb_z, fg_55, gf0_33, gf0_34, \
                         gf0_35, gf1_33, gf1_34, gf1_35, gg_68, gg_69, \
                         gg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_1 * gf0_33[k]
                   - f_2 * gf1_33[k]
                   + pb_y[k] * gg_68[k];

        t_190[k] = f_0 * fg_55[k]
                   + pb_z[k] * gg_68[k];

        t_191[k] = f_5 * gf0_34[k]
                   - f_6 * gf1_34[k]
                   + pb_y[k] * gg_69[k];

        t_192[k] = f_3 * gf0_35[k]
                   - f_4 * gf1_35[k]
                   + pb_y[k] * gg_70[k];
    }

#pragma omp simd aligned(t_193, t_194, pb_y, pb_z, fg_58, gf0_35, gf1_35, \
                         gg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_y[k] * gg_71[k];

        t_194[k] = f_0 * fg_58[k]
                   + f_1 * gf0_35[k]
                   - f_2 * gf1_35[k]
                   + pb_z[k] * gg_71[k];
    }
}

auto
compute_prim_gh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_37 = buffer.data(gf0 + 37);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_41 = buffer.data(gf0 + 41);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_37 = buffer.data(gf1 + 37);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_40 = buffer.data(gf1 + 40);
    const auto *gf1_41 = buffer.data(gf1 + 41);
    const auto *gf1_46 = buffer.data(gf1 + 46);
    const auto *gf1_47 = buffer.data(gf1 + 47);
    const auto *gf1_48 = buffer.data(gf1 + 48);
    const auto *gf1_49 = buffer.data(gf1 + 49);
    const auto *gf1_50 = buffer.data(gf1 + 50);
    const auto *gf1_51 = buffer.data(gf1 + 51);
    const auto *gf1_56 = buffer.data(gf1 + 56);
    const auto *gf1_58 = buffer.data(gf1 + 58);
    const auto *gf1_59 = buffer.data(gf1 + 59);
    const auto *gf1_60 = buffer.data(gf1 + 60);
    const auto *gf1_61 = buffer.data(gf1 + 61);
    const auto *gf1_62 = buffer.data(gf1 + 62);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf0_0, gf0_1, gf1_0, \
                         gf1_1, gg_0, gg_1, gg_2, gg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_2[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];

        t_3[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, fg_5, fg_8, gf0_2, gf0_3, \
                         gf1_2, gf1_3, gg_4, gg_5, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_5[k] = f_0 * fg_5[k]
                 + pb_x[k] * gg_5[k];

        t_6[k] = f_0 * fg_8[k]
                 + pb_x[k] * gg_8[k];

        t_7[k] = f_1 * gf0_3[k]
                 - f_2 * gf1_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pb_y, pb_z, fh_0, gf0_4, gf0_5, gf1_5, \
                         gf1_6, gg_6, gg_7, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gf0_4[k]
                 - f_6 * gf1_5[k]
                 + pb_y[k] * gg_6[k];

        t_9[k] = f_3 * gf0_5[k]
                 - f_4 * gf1_6[k]
                 + pb_y[k] * gg_7[k];

        t_10[k] = f_1 * gf0_5[k]
                  - f_2 * gf1_6[k]
                  + pb_z[k] * gg_8[k];

        t_11[k] = pa_y[k] * fh_0[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, fg_0, fg_1, fg_3, fg_11, \
                         fh_1, fh_3, gg_9, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * fg_0[k]
                  + pb_y[k] * gg_9[k];

        t_13[k] = f_8 * fg_1[k]
                  + pa_y[k] * fh_1[k];

        t_14[k] = f_9 * fg_3[k]
                  + pa_y[k] * fh_3[k];

        t_15[k] = f_9 * fg_11[k]
                  + pb_x[k] * gg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_z, fg_0, fg_2, fg_5, fh_0, \
                         fh_2, fh_5, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_10 * fg_5[k]
                  + pa_y[k] * fh_5[k];

        t_17[k] = pa_z[k] * fh_0[k];

        t_18[k] = f_7 * fg_0[k]
                  + pb_z[k] * gg_12[k];

        t_19[k] = f_8 * fg_2[k]
                  + pa_z[k] * fh_2[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_x, fg_4, fg_6, fg_7, fg_16, fh_4, \
                         fh_6, fh_7, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_9 * fg_4[k]
                  + pa_z[k] * fh_4[k];

        t_21[k] = f_9 * fg_16[k]
                  + pb_x[k] * gg_16[k];

        t_22[k] = f_8 * fg_6[k]
                  + pa_z[k] * fh_6[k];

        t_23[k] = f_9 * fg_7[k]
                  + pa_z[k] * fh_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pa_z, pb_y, dh0_0, dh1_0, fg_8, fg_9, fh_8, \
                         fh_9, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_10 * fg_8[k]
                  + pa_z[k] * fh_8[k];

        t_25[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_y[k] * fh_9[k];

        t_26[k] = f_8 * fg_9[k]
                  + pb_y[k] * gg_17[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pb_z, fg_18, fg_19, gf0_8, gf0_10, gf0_11, \
                         gf1_14, gf1_16, gf1_17, gg_18, gg_19, gg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_8 * fg_18[k]
                  + f_5 * gf0_10[k]
                  - f_6 * gf1_16[k]
                  + pb_x[k] * gg_19[k];

        t_28[k] = f_3 * gf0_8[k]
                  - f_4 * gf1_14[k]
                  + pb_z[k] * gg_18[k];

        t_29[k] = f_8 * fg_19[k]
                  + f_3 * gf0_11[k]
                  - f_4 * gf1_17[k]
                  + pb_x[k] * gg_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pb_x, pb_z, dh0_1, dh1_1, fg_20, fh_11, \
                         gf0_9, gf1_15, gg_20, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * gf0_9[k]
                  - f_6 * gf1_15[k]
                  + pb_z[k] * gg_20[k];

        t_31[k] = f_8 * fg_20[k]
                  + pb_x[k] * gg_22[k];

        t_32[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_1[k]
                  + pa_x[k] * fh_11[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_z, gf0_11, gf0_12, gf0_13, gf1_17, gf1_18, \
                         gf1_19, gg_23, gg_24, gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_17[k]
                  + pb_z[k] * gg_23[k];

        t_34[k] = f_5 * gf0_12[k]
                  - f_6 * gf1_18[k]
                  + pb_z[k] * gg_24[k];

        t_35[k] = f_1 * gf0_13[k]
                  - f_2 * gf1_19[k]
                  + pb_z[k] * gg_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_z, pb_y, pb_z, dh0_0, dh1_0, fg_12, fh_10, \
                         gf0_14, gf1_22, gg_26, gg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_z[k] * fh_10[k];

        t_37[k] = f_8 * fg_12[k]
                  + pb_z[k] * gg_26[k];

        t_38[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_22[k]
                  + pb_y[k] * gg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, fg_23, fg_24, gf0_15, gf0_16, gf0_19, \
                         gf1_23, gf1_24, gf1_27, gg_29, gg_30, gg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * fg_23[k]
                  + f_5 * gf0_16[k]
                  - f_6 * gf1_24[k]
                  + pb_x[k] * gg_30[k];

        t_40[k] = f_5 * gf0_15[k]
                  - f_6 * gf1_23[k]
                  + pb_y[k] * gg_29[k];

        t_41[k] = f_8 * fg_24[k]
                  + f_3 * gf0_19[k]
                  - f_4 * gf1_27[k]
                  + pb_x[k] * gg_31[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, fg_25, gf0_17, gf0_18, gf1_25, gf1_26, \
                         gg_32, gg_33, gg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_8 * fg_25[k]
                  + pb_x[k] * gg_35[k];

        t_43[k] = f_1 * gf0_17[k]
                  - f_2 * gf1_25[k]
                  + pb_y[k] * gg_32[k];

        t_44[k] = f_5 * gf0_18[k]
                  - f_6 * gf1_26[k]
                  + pb_y[k] * gg_33[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_y, dh0_2, dh1_2, fg_17, fg_26, \
                         fh_12, fh_13, gf0_19, gf1_27, gg_34, gg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * gf0_19[k]
                  - f_4 * gf1_27[k]
                  + pb_y[k] * gg_34[k];

        t_46[k] = f_11 * dh0_2[k]
                  - f_12 * dh1_2[k]
                  + pa_x[k] * fh_12[k];

        t_47[k] = f_10 * fg_26[k]
                  + pa_x[k] * fh_13[k];

        t_48[k] = f_9 * fg_17[k]
                  + pb_y[k] * gg_36[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_x, pb_x, fg_28, fg_30, fg_31, fg_48, \
                         fh_14, fh_16, fh_18, fh_24, gg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * fg_28[k]
                  + pa_x[k] * fh_14[k];

        t_50[k] = f_8 * fg_30[k]
                  + pa_x[k] * fh_16[k];

        t_51[k] = f_7 * fg_31[k]
                  + pb_x[k] * gg_38[k];

        t_52[k] = pa_x[k] * fh_18[k];

        t_53[k] = f_10 * fg_48[k]
                  + pa_x[k] * fh_24[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_x, pb_z, fg_21, fg_52, fg_53, fg_58, \
                         fh_26, fh_28, gg_39, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * fg_21[k]
                  + pb_z[k] * gg_39[k];

        t_55[k] = f_9 * fg_52[k]
                  + pa_x[k] * fh_26[k];

        t_56[k] = f_8 * fg_53[k]
                  + pa_x[k] * fh_28[k];

        t_57[k] = f_7 * fg_58[k]
                  + pb_x[k] * gg_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pb_x, pb_y, fg_26, fh_32, gf0_22, \
                         gf0_23, gf1_35, gf1_37, gg_43, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_x[k] * fh_32[k];

        t_59[k] = f_1 * gf0_22[k]
                  - f_2 * gf1_35[k]
                  + pb_x[k] * gg_43[k];

        t_60[k] = f_0 * fg_26[k]
                  + pb_y[k] * gg_43[k];

        t_61[k] = f_5 * gf0_23[k]
                  - f_6 * gf1_37[k]
                  + pb_x[k] * gg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pb_x, gf0_24, gf0_25, gf0_27, gf1_38, gf1_39, \
                         gf1_41, gg_45, gg_46, gg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_5 * gf0_24[k]
                  - f_6 * gf1_38[k]
                  + pb_x[k] * gg_45[k];

        t_63[k] = f_3 * gf0_25[k]
                  - f_4 * gf1_39[k]
                  + pb_x[k] * gg_46[k];

        t_64[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_41[k]
                  + pb_x[k] * gg_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pb_z, fg_31, fg_35, gf0_25, gf0_26, \
                         gf1_39, gf1_40, gg_48, gg_49, gg_50, gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * fg_31[k]
                  + f_1 * gf0_25[k]
                  - f_2 * gf1_39[k]
                  + pb_y[k] * gg_48[k];

        t_66[k] = f_3 * gf0_25[k]
                  - f_4 * gf1_39[k]
                  + pb_z[k] * gg_49[k];

        t_67[k] = f_5 * gf0_26[k]
                  - f_6 * gf1_40[k]
                  + pb_z[k] * gg_50[k];

        t_68[k] = f_0 * fg_35[k]
                  + pb_y[k] * gg_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pb_z, fg_27, fg_29, fh_15, fh_17, \
                         fh_18, gf0_27, gf1_41, gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * gf0_27[k]
                  - f_2 * gf1_41[k]
                  + pb_z[k] * gg_52[k];

        t_70[k] = f_8 * fg_27[k]
                  + pa_z[k] * fh_15[k];

        t_71[k] = f_9 * fg_29[k]
                  + pa_z[k] * fh_17[k];

        t_72[k] = pa_z[k] * fh_18[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pb_y, pb_z, fg_31, fg_32, fg_33, fg_41, \
                         fh_19, fh_20, gg_54, gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * fg_31[k]
                  + pb_z[k] * gg_54[k];

        t_74[k] = f_8 * fg_32[k]
                  + pa_z[k] * fh_19[k];

        t_75[k] = f_9 * fg_33[k]
                  + pa_z[k] * fh_20[k];

        t_76[k] = f_9 * fg_41[k]
                  + pb_y[k] * gg_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_z, pb_x, fg_35, fh_21, gf0_29, gf0_30, gf1_46, \
                         gf1_47, gg_59, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_10 * fg_35[k]
                  + pa_z[k] * fh_21[k];

        t_78[k] = f_1 * gf0_29[k]
                  - f_2 * gf1_46[k]
                  + pb_x[k] * gg_59[k];

        t_79[k] = f_5 * gf0_30[k]
                  - f_6 * gf1_47[k]
                  + pb_x[k] * gg_60[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, gf0_31, gf0_32, gf0_34, gf1_48, gf1_49, \
                         gf1_51, gg_61, gg_62, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_5 * gf0_31[k]
                  - f_6 * gf1_48[k]
                  + pb_x[k] * gg_61[k];

        t_81[k] = f_3 * gf0_32[k]
                  - f_4 * gf1_49[k]
                  + pb_x[k] * gg_62[k];

        t_82[k] = f_3 * gf0_34[k]
                  - f_4 * gf1_51[k]
                  + pb_x[k] * gg_63[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_z, pb_y, pb_z, dh0_1, dh1_1, fg_37, fg_45, \
                         fh_22, gf0_33, gf1_50, gg_64, gg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_1[k]
                  + pa_z[k] * fh_22[k];

        t_84[k] = f_8 * fg_37[k]
                  + pb_z[k] * gg_64[k];

        t_85[k] = f_8 * fg_45[k]
                  + f_5 * gf0_33[k]
                  - f_6 * gf1_50[k]
                  + pb_y[k] * gg_66[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_y, pb_y, dh0_2, dh1_2, fg_46, fg_47, fh_23, \
                         gf0_34, gf1_51, gg_67, gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_8 * fg_46[k]
                  + f_3 * gf0_34[k]
                  - f_4 * gf1_51[k]
                  + pb_y[k] * gg_67[k];

        t_87[k] = f_8 * fg_47[k]
                  + pb_y[k] * gg_68[k];

        t_88[k] = f_11 * dh0_2[k]
                  - f_12 * dh1_2[k]
                  + pa_y[k] * fh_23[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pb_z, fg_43, fg_49, fg_51, fg_54, \
                         fh_25, fh_27, fh_29, gg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_8 * fg_49[k]
                  + pa_y[k] * fh_25[k];

        t_90[k] = f_9 * fg_51[k]
                  + pa_y[k] * fh_27[k];

        t_91[k] = f_10 * fg_54[k]
                  + pa_y[k] * fh_29[k];

        t_92[k] = f_9 * fg_43[k]
                  + pb_z[k] * gg_70[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_y, fg_56, fg_57, fg_58, fh_30, \
                         fh_31, fh_32, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_9 * fg_56[k]
                  + pa_y[k] * fh_30[k];

        t_94[k] = f_8 * fg_57[k]
                  + pa_y[k] * fh_31[k];

        t_95[k] = f_7 * fg_58[k]
                  + pb_y[k] * gg_74[k];

        t_96[k] = pa_y[k] * fh_32[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_z, fg_48, gf0_36, gf0_37, gf0_38, \
                         gf1_56, gf1_58, gf1_59, gg_75, gg_77, gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * gf0_36[k]
                  - f_2 * gf1_56[k]
                  + pb_x[k] * gg_75[k];

        t_98[k] = f_0 * fg_48[k]
                  + pb_z[k] * gg_75[k];

        t_99[k] = f_5 * gf0_37[k]
                  - f_6 * gf1_58[k]
                  + pb_x[k] * gg_77[k];

        t_100[k] = f_5 * gf0_38[k]
                   - f_6 * gf1_59[k]
                   + pb_x[k] * gg_78[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, pb_y, pb_z, fg_54, gf0_39, gf0_41, \
                         gf1_60, gf1_62, gg_79, gg_80, gg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_3 * gf0_39[k]
                   - f_4 * gf1_60[k]
                   + pb_x[k] * gg_79[k];

        t_102[k] = f_3 * gf0_41[k]
                   - f_4 * gf1_62[k]
                   + pb_x[k] * gg_80[k];

        t_103[k] = f_1 * gf0_39[k]
                   - f_2 * gf1_60[k]
                   + pb_y[k] * gg_81[k];

        t_104[k] = f_0 * fg_54[k]
                   + pb_z[k] * gg_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pb_y, pb_z, fg_58, gf0_40, gf0_41, gf1_61, \
                         gf1_62, gg_83, gg_84, gg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_5 * gf0_40[k]
                   - f_6 * gf1_61[k]
                   + pb_y[k] * gg_83[k];

        t_106[k] = f_3 * gf0_41[k]
                   - f_4 * gf1_62[k]
                   + pb_y[k] * gg_84[k];

        t_107[k] = f_0 * fg_58[k]
                   + f_1 * gf0_41[k]
                   - f_2 * gf1_62[k]
                   + pb_z[k] * gg_85[k];
    }
}

auto
compute_prim_gh_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_79 = buffer.data(fh + 79);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_37 = buffer.data(gf0 + 37);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_41 = buffer.data(gf0 + 41);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_31 = buffer.data(gf1 + 31);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_40 = buffer.data(gf1 + 40);
    const auto *gf1_41 = buffer.data(gf1 + 41);
    const auto *gf1_42 = buffer.data(gf1 + 42);
    const auto *gf1_43 = buffer.data(gf1 + 43);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fg_0, gf0_0, gf1_0, gg_0, \
                         gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pb_y[k] * gg_0[k];

        t_2[k] = pb_z[k] * gg_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_4[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gf0_1, gf0_2, gf0_3, gf1_1, \
                         gf1_2, gf1_3, gg_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = pb_z[k] * gg_3[k];

        t_7[k] = pb_y[k] * gg_4[k];

        t_8[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_9[k] = f_1 * gf0_3[k]
                 - f_2 * gf1_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pb_z, gf0_4, gf0_5, gf1_4, gf1_5, \
                         gg_5, gg_7, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * gg_5[k];

        t_11[k] = f_5 * gf0_4[k]
                  - f_6 * gf1_4[k]
                  + pb_y[k] * gg_7[k];

        t_12[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pb_y[k] * gg_8[k];

        t_13[k] = pb_y[k] * gg_9[k];

        t_14[k] = f_1 * gf0_5[k]
                  - f_2 * gf1_5[k]
                  + pb_z[k] * gg_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_y, fg_1, fg_3, fg_5, fh_0, \
                         fh_3, fh_4, fh_5, fh_7, fh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * fh_0[k];

        t_16[k] = f_7 * fg_1[k]
                  + pa_y[k] * fh_3[k];

        t_17[k] = pa_y[k] * fh_4[k];

        t_18[k] = f_8 * fg_3[k]
                  + pa_y[k] * fh_5[k];

        t_19[k] = pa_y[k] * fh_7[k];

        t_20[k] = f_9 * fg_5[k]
                  + pa_y[k] * fh_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pa_z, pb_y, fg_7, fg_8, fg_9, \
                         fh_0, fh_10, fh_11, fh_12, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * fg_7[k]
                  + pa_y[k] * fh_10[k];

        t_22[k] = f_7 * fg_8[k]
                  + pa_y[k] * fh_11[k];

        t_23[k] = f_10 * fg_9[k]
                  + pb_y[k] * gg_12[k];

        t_24[k] = pa_y[k] * fh_12[k];

        t_25[k] = pa_z[k] * fh_0[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_z, pb_y, pb_z, fg_0, fg_2, fh_3, \
                         fh_4, fh_5, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_10 * fg_0[k]
                  + pb_z[k] * gg_13[k];

        t_27[k] = pa_z[k] * fh_3[k];

        t_28[k] = f_7 * fg_2[k]
                  + pa_z[k] * fh_4[k];

        t_29[k] = pa_z[k] * fh_5[k];

        t_30[k] = pb_y[k] * gg_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_z, pb_z, fg_4, fg_5, fg_6, fg_7, \
                         fh_7, fh_8, fh_10, fh_11, gg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_8 * fg_4[k]
                  + pa_z[k] * fh_7[k];

        t_32[k] = pa_z[k] * fh_8[k];

        t_33[k] = f_10 * fg_5[k]
                  + pb_z[k] * gg_15[k];

        t_34[k] = f_7 * fg_6[k]
                  + pa_z[k] * fh_10[k];

        t_35[k] = f_8 * fg_7[k]
                  + pa_z[k] * fh_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, pb_y, pb_z, dh0_0, dh1_0, fg_9, \
                         fh_12, fh_13, gg_18, gg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * gg_18[k];

        t_37[k] = f_9 * fg_9[k]
                  + pa_z[k] * fh_12[k];

        t_38[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_y[k] * fh_13[k];

        t_39[k] = pb_z[k] * gg_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, pb_z, fg_19, fg_20, gf0_6, gf0_8, gf0_9, \
                         gf1_8, gf1_10, gf1_11, gg_20, gg_21, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * fg_19[k]
                  + f_5 * gf0_8[k]
                  - f_6 * gf1_10[k]
                  + pb_x[k] * gg_21[k];

        t_41[k] = f_3 * gf0_6[k]
                  - f_4 * gf1_8[k]
                  + pb_z[k] * gg_20[k];

        t_42[k] = f_7 * fg_20[k]
                  + f_3 * gf0_9[k]
                  - f_4 * gf1_11[k]
                  + pb_x[k] * gg_23[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_z, dh0_1, dh1_1, fg_21, fh_28, \
                         gf0_7, gf1_9, gg_21, gg_22, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_z[k] * gg_21[k];

        t_44[k] = f_5 * gf0_7[k]
                  - f_6 * gf1_9[k]
                  + pb_z[k] * gg_22[k];

        t_45[k] = f_7 * fg_21[k]
                  + pb_x[k] * gg_24[k];

        t_46[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_1[k]
                  + pa_x[k] * fh_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pb_z, fg_12, gf0_9, gf0_10, gf1_11, \
                         gf1_12, gg_24, gg_25, gg_26, gg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * gg_24[k];

        t_48[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_11[k]
                  + pb_z[k] * gg_25[k];

        t_49[k] = f_5 * gf0_10[k]
                  - f_6 * gf1_12[k]
                  + pb_z[k] * gg_26[k];

        t_50[k] = f_7 * fg_12[k]
                  + pb_y[k] * gg_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_y, pa_z, pb_z, fh_14, fh_15, fh_18, \
                         fh_19, gf0_11, gf1_13, gg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * gf0_11[k]
                  - f_2 * gf1_13[k]
                  + pb_z[k] * gg_27[k];

        t_52[k] = pa_y[k] * fh_18[k];

        t_53[k] = pa_z[k] * fh_14[k];

        t_54[k] = pa_y[k] * fh_19[k];

        t_55[k] = pa_z[k] * fh_15[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pa_z, pb_z, fg_11, fg_15, fg_16, \
                         fh_16, fh_20, fh_21, fh_22, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * fh_20[k];

        t_57[k] = pa_z[k] * fh_16[k];

        t_58[k] = f_10 * fg_11[k]
                  + pb_z[k] * gg_28[k];

        t_59[k] = f_8 * fg_15[k]
                  + pa_y[k] * fh_21[k];

        t_60[k] = f_7 * fg_16[k]
                  + pa_y[k] * fh_22[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pa_z, pb_y, dh0_0, dh1_0, fg_17, fh_17, \
                         fh_23, gg_29, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_10 * fg_17[k]
                  + pb_y[k] * gg_29[k];

        t_62[k] = pa_y[k] * fh_23[k];

        t_63[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_z[k] * fh_17[k];

        t_64[k] = pb_y[k] * gg_30[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_x, pb_y, pb_z, fg_13, fg_23, gf0_12, gf0_14, \
                         gf1_14, gf1_16, gg_30, gg_31, gg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_7 * fg_13[k]
                  + pb_z[k] * gg_30[k];

        t_66[k] = f_3 * gf0_12[k]
                  - f_4 * gf1_14[k]
                  + pb_y[k] * gg_31[k];

        t_67[k] = f_7 * fg_23[k]
                  + f_5 * gf0_14[k]
                  - f_6 * gf1_16[k]
                  + pb_x[k] * gg_33[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, fg_24, fg_25, gf0_13, gf0_17, \
                         gf1_15, gf1_19, gg_32, gg_33, gg_34, gg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_15[k]
                  + pb_y[k] * gg_32[k];

        t_69[k] = pb_y[k] * gg_33[k];

        t_70[k] = f_7 * fg_24[k]
                  + f_3 * gf0_17[k]
                  - f_4 * gf1_19[k]
                  + pb_x[k] * gg_34[k];

        t_71[k] = f_7 * fg_25[k]
                  + pb_x[k] * gg_38[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_y, pb_z, fg_14, gf0_15, gf0_16, gf0_17, \
                         gf1_17, gf1_18, gf1_19, gg_35, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * gf0_15[k]
                  - f_2 * gf1_17[k]
                  + pb_y[k] * gg_35[k];

        t_73[k] = f_7 * fg_14[k]
                  + pb_z[k] * gg_35[k];

        t_74[k] = f_5 * gf0_16[k]
                  - f_6 * gf1_18[k]
                  + pb_y[k] * gg_36[k];

        t_75[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_19[k]
                  + pb_y[k] * gg_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_x, pb_y, dh0_2, dh1_2, fg_26, fg_28, \
                         fh_34, fh_35, fh_37, gg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * gg_38[k];

        t_77[k] = f_11 * dh0_2[k]
                  - f_12 * dh1_2[k]
                  + pa_x[k] * fh_34[k];

        t_78[k] = f_9 * fg_26[k]
                  + pa_x[k] * fh_35[k];

        t_79[k] = f_8 * fg_28[k]
                  + pa_x[k] * fh_37[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_x, pb_x, fg_29, fg_30, fg_31, fg_32, \
                         fh_38, fh_39, fh_41, fh_43, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_8 * fg_29[k]
                  + pa_x[k] * fh_38[k];

        t_81[k] = f_7 * fg_30[k]
                  + pa_x[k] * fh_39[k];

        t_82[k] = f_7 * fg_31[k]
                  + pa_x[k] * fh_41[k];

        t_83[k] = f_10 * fg_32[k]
                  + pb_x[k] * gg_42[k];

        t_84[k] = pa_x[k] * fh_43[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, t_90, pa_x, pa_z, pb_z, fg_18, fh_24, \
                         fh_45, fh_46, fh_47, fh_48, gg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_x[k] * fh_45[k];

        t_86[k] = pa_x[k] * fh_46[k];

        t_87[k] = pa_x[k] * fh_47[k];

        t_88[k] = pa_x[k] * fh_48[k];

        t_89[k] = pa_z[k] * fh_24[k];

        t_90[k] = f_10 * fg_18[k]
                  + pb_z[k] * gg_43[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, t_96, pa_x, pa_z, fg_37, fg_38, fh_25, \
                         fh_26, fh_49, fh_50, fh_52, fh_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pa_z[k] * fh_25[k];

        t_92[k] = f_8 * fg_37[k]
                  + pa_x[k] * fh_49[k];

        t_93[k] = pa_z[k] * fh_26[k];

        t_94[k] = f_7 * fg_38[k]
                  + pa_x[k] * fh_50[k];

        t_95[k] = pa_x[k] * fh_52[k];

        t_96[k] = pa_x[k] * fh_53[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, t_102, pa_x, pa_y, fg_41, fh_29, \
                         fh_30, fh_54, fh_55, fh_56, fh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pa_x[k] * fh_54[k];

        t_98[k] = pa_x[k] * fh_55[k];

        t_99[k] = pa_x[k] * fh_56[k];

        t_100[k] = pa_y[k] * fh_29[k];

        t_101[k] = pa_y[k] * fh_30[k];

        t_102[k] = f_8 * fg_41[k]
                   + pa_x[k] * fh_57[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, pa_x, pa_y, fg_42, fh_31, \
                         fh_32, fh_58, fh_59, fh_60, fh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pa_y[k] * fh_31[k];

        t_104[k] = f_7 * fg_42[k]
                   + pa_x[k] * fh_58[k];

        t_105[k] = pa_y[k] * fh_32[k];

        t_106[k] = pa_x[k] * fh_59[k];

        t_107[k] = pa_x[k] * fh_60[k];

        t_108[k] = pa_x[k] * fh_61[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pa_x, pb_z, fg_22, fg_47, fg_49, \
                         fh_62, fh_63, fh_65, fh_68, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * fh_62[k];

        t_110[k] = pa_x[k] * fh_63[k];

        t_111[k] = f_9 * fg_47[k]
                   + pa_x[k] * fh_65[k];

        t_112[k] = f_8 * fg_22[k]
                   + pb_z[k] * gg_44[k];

        t_113[k] = f_8 * fg_49[k]
                   + pa_x[k] * fh_68[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pa_x, pb_x, fg_50, fg_51, fg_52, \
                         fg_56, fh_69, fh_70, fh_72, fh_74, gg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_8 * fg_50[k]
                   + pa_x[k] * fh_69[k];

        t_115[k] = f_7 * fg_51[k]
                   + pa_x[k] * fh_70[k];

        t_116[k] = f_7 * fg_52[k]
                   + pa_x[k] * fh_72[k];

        t_117[k] = f_10 * fg_56[k]
                   + pb_x[k] * gg_47[k];

        t_118[k] = pa_x[k] * fh_74[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, pa_x, pb_x, pb_z, fh_75, \
                         fh_76, fh_77, fh_79, gf0_22, gf1_24, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * fh_75[k];

        t_120[k] = pa_x[k] * fh_76[k];

        t_121[k] = pa_x[k] * fh_77[k];

        t_122[k] = pa_x[k] * fh_79[k];

        t_123[k] = f_1 * gf0_22[k]
                   - f_2 * gf1_24[k]
                   + pb_x[k] * gg_48[k];

        t_124[k] = pb_z[k] * gg_48[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_x, pb_z, gf0_23, gf0_24, gf0_25, \
                         gf1_25, gf1_26, gf1_27, gg_50, gg_51, gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_5 * gf0_23[k]
                   - f_6 * gf1_25[k]
                   + pb_x[k] * gg_50[k];

        t_126[k] = f_5 * gf0_24[k]
                   - f_6 * gf1_26[k]
                   + pb_x[k] * gg_51[k];

        t_127[k] = f_3 * gf0_25[k]
                   - f_4 * gf1_27[k]
                   + pb_x[k] * gg_52[k];

        t_128[k] = pb_z[k] * gg_50[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_y, fg_32, gf0_25, gf0_27, \
                         gf1_27, gf1_29, gg_53, gg_54, gg_56, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_3 * gf0_27[k]
                   - f_4 * gf1_29[k]
                   + pb_x[k] * gg_53[k];

        t_130[k] = pb_x[k] * gg_54[k];

        t_131[k] = pb_x[k] * gg_56[k];

        t_132[k] = pb_x[k] * gg_57[k];

        t_133[k] = f_0 * fg_32[k]
                   + f_1 * gf0_25[k]
                   - f_2 * gf1_27[k]
                   + pb_y[k] * gg_54[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pb_y, pb_z, fg_35, gf0_25, gf0_26, \
                         gf1_27, gf1_28, gg_54, gg_55, gg_56, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pb_z[k] * gg_54[k];

        t_135[k] = f_3 * gf0_25[k]
                   - f_4 * gf1_27[k]
                   + pb_z[k] * gg_55[k];

        t_136[k] = f_5 * gf0_26[k]
                   - f_6 * gf1_28[k]
                   + pb_z[k] * gg_56[k];

        t_137[k] = f_0 * fg_35[k]
                   + pb_y[k] * gg_57[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pa_z, pb_z, fg_26, fg_27, fh_35, \
                         fh_37, fh_38, gf0_27, gf1_29, gg_57, gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_1 * gf0_27[k]
                   - f_2 * gf1_29[k]
                   + pb_z[k] * gg_57[k];

        t_139[k] = pa_z[k] * fh_35[k];

        t_140[k] = f_10 * fg_26[k]
                   + pb_z[k] * gg_58[k];

        t_141[k] = pa_z[k] * fh_37[k];

        t_142[k] = f_7 * fg_27[k]
                   + pa_z[k] * fh_38[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, pa_z, pb_x, fg_29, fh_39, fh_41, \
                         fh_43, gg_62, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = pa_z[k] * fh_39[k];

        t_144[k] = f_8 * fg_29[k]
                   + pa_z[k] * fh_41[k];

        t_145[k] = pb_x[k] * gg_62[k];

        t_146[k] = pb_x[k] * gg_63[k];

        t_147[k] = pa_z[k] * fh_43[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_z, pb_y, pb_z, fg_32, fg_33, fg_34, \
                         fg_40, fh_45, fh_46, gg_61, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * fg_32[k]
                   + pb_z[k] * gg_61[k];

        t_149[k] = f_7 * fg_33[k]
                   + pa_z[k] * fh_45[k];

        t_150[k] = f_8 * fg_34[k]
                   + pa_z[k] * fh_46[k];

        t_151[k] = f_8 * fg_40[k]
                   + pb_y[k] * gg_63[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_z, pb_x, pb_z, fg_35, fg_36, fh_48, \
                         gf0_28, gf0_29, gf1_31, gf1_32, gg_64, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_9 * fg_35[k]
                   + pa_z[k] * fh_48[k];

        t_153[k] = f_1 * gf0_28[k]
                   - f_2 * gf1_31[k]
                   + pb_x[k] * gg_64[k];

        t_154[k] = f_7 * fg_36[k]
                   + pb_z[k] * gg_64[k];

        t_155[k] = f_5 * gf0_29[k]
                   - f_6 * gf1_32[k]
                   + pb_x[k] * gg_65[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, gf0_30, gf0_31, gf0_33, gf1_33, \
                         gf1_34, gf1_36, gg_66, gg_67, gg_68, gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_5 * gf0_30[k]
                   - f_6 * gf1_33[k]
                   + pb_x[k] * gg_66[k];

        t_157[k] = f_3 * gf0_31[k]
                   - f_4 * gf1_34[k]
                   + pb_x[k] * gg_67[k];

        t_158[k] = f_3 * gf0_33[k]
                   - f_4 * gf1_36[k]
                   + pb_x[k] * gg_68[k];

        t_159[k] = pb_x[k] * gg_69[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pb_x, pb_z, dh0_1, dh1_1, fg_39, \
                         fh_51, gg_69, gg_70, gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_x[k] * gg_70[k];

        t_161[k] = pb_x[k] * gg_72[k];

        t_162[k] = f_11 * dh0_1[k]
                   - f_12 * dh1_1[k]
                   + pa_z[k] * fh_51[k];

        t_163[k] = f_7 * fg_39[k]
                   + pb_z[k] * gg_69[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pb_y, fg_44, fg_45, fg_46, gf0_32, gf0_33, \
                         gf1_35, gf1_36, gg_70, gg_71, gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_7 * fg_44[k]
                   + f_5 * gf0_32[k]
                   - f_6 * gf1_35[k]
                   + pb_y[k] * gg_70[k];

        t_165[k] = f_7 * fg_45[k]
                   + f_3 * gf0_33[k]
                   - f_4 * gf1_36[k]
                   + pb_y[k] * gg_71[k];

        t_166[k] = f_7 * fg_46[k]
                   + pb_y[k] * gg_72[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, dh0_2, dh1_2, fg_48, fh_64, \
                         fh_65, fh_67, fh_68, fh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_11 * dh0_2[k]
                   - f_12 * dh1_2[k]
                   + pa_y[k] * fh_64[k];

        t_168[k] = pa_y[k] * fh_65[k];

        t_169[k] = pa_y[k] * fh_67[k];

        t_170[k] = f_7 * fg_48[k]
                   + pa_y[k] * fh_68[k];

        t_171[k] = pa_y[k] * fh_69[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, pa_y, pb_x, fg_49, fg_53, fh_70, \
                         fh_72, fh_74, gg_75, gg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_8 * fg_49[k]
                   + pa_y[k] * fh_70[k];

        t_173[k] = pa_y[k] * fh_72[k];

        t_174[k] = pb_x[k] * gg_75[k];

        t_175[k] = pb_x[k] * gg_76[k];

        t_176[k] = f_9 * fg_53[k]
                   + pa_y[k] * fh_74[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pb_y, pb_z, fg_43, fg_54, fg_55, \
                         fg_56, fh_76, fh_77, gg_75, gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_8 * fg_43[k]
                   + pb_z[k] * gg_75[k];

        t_178[k] = f_8 * fg_54[k]
                   + pa_y[k] * fh_76[k];

        t_179[k] = f_7 * fg_55[k]
                   + pa_y[k] * fh_77[k];

        t_180[k] = f_10 * fg_56[k]
                   + pb_y[k] * gg_78[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_y, pb_x, pb_y, pb_z, fg_47, fh_79, \
                         gf0_36, gf1_39, gg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = pa_y[k] * fh_79[k];

        t_182[k] = f_1 * gf0_36[k]
                   - f_2 * gf1_39[k]
                   + pb_x[k] * gg_79[k];

        t_183[k] = pb_y[k] * gg_79[k];

        t_184[k] = f_0 * fg_47[k]
                   + pb_z[k] * gg_79[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, pb_y, gf0_37, gf0_38, gf0_39, \
                         gf1_40, gf1_41, gf1_42, gg_81, gg_82, gg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_5 * gf0_37[k]
                   - f_6 * gf1_40[k]
                   + pb_x[k] * gg_81[k];

        t_186[k] = f_5 * gf0_38[k]
                   - f_6 * gf1_41[k]
                   + pb_x[k] * gg_82[k];

        t_187[k] = f_3 * gf0_39[k]
                   - f_4 * gf1_42[k]
                   + pb_x[k] * gg_83[k];

        t_188[k] = pb_y[k] * gg_82[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pb_x, pb_y, gf0_39, gf0_41, \
                         gf1_42, gf1_44, gg_84, gg_85, gg_86, gg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_3 * gf0_41[k]
                   - f_4 * gf1_44[k]
                   + pb_x[k] * gg_84[k];

        t_190[k] = pb_x[k] * gg_85[k];

        t_191[k] = pb_x[k] * gg_86[k];

        t_192[k] = pb_x[k] * gg_88[k];

        t_193[k] = f_1 * gf0_39[k]
                   - f_2 * gf1_42[k]
                   + pb_y[k] * gg_85[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pb_y, pb_z, fg_53, gf0_40, gf0_41, \
                         gf1_43, gf1_44, gg_85, gg_86, gg_87, gg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_0 * fg_53[k]
                   + pb_z[k] * gg_85[k];

        t_195[k] = f_5 * gf0_40[k]
                   - f_6 * gf1_43[k]
                   + pb_y[k] * gg_86[k];

        t_196[k] = f_3 * gf0_41[k]
                   - f_4 * gf1_44[k]
                   + pb_y[k] * gg_87[k];

        t_197[k] = pb_y[k] * gg_88[k];
    }

#pragma omp simd aligned(t_198, pb_z, fg_56, gf0_41, gf1_44, gg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_0 * fg_56[k]
                   + f_1 * gf0_41[k]
                   - f_2 * gf1_44[k]
                   + pb_z[k] * gg_88[k];
    }
}

auto
compute_prim_gh_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_26 = buffer.data(dh1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_44 = buffer.data(fh + 44);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_35 = buffer.data(gf0 + 35);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_41 = buffer.data(gf0 + 41);
    const auto *gf0_42 = buffer.data(gf0 + 42);
    const auto *gf0_43 = buffer.data(gf0 + 43);
    const auto *gf0_44 = buffer.data(gf0 + 44);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_37 = buffer.data(gf1 + 37);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_42 = buffer.data(gf1 + 42);
    const auto *gf1_43 = buffer.data(gf1 + 43);
    const auto *gf1_44 = buffer.data(gf1 + 44);
    const auto *gf1_45 = buffer.data(gf1 + 45);
    const auto *gf1_46 = buffer.data(gf1 + 46);
    const auto *gf1_47 = buffer.data(gf1 + 47);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fg_0, gf0_0, gf1_0, gg_0, \
                         gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pb_y[k] * gg_0[k];

        t_2[k] = pb_z[k] * gg_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_4[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gf0_1, gf0_2, gf0_3, gf1_1, gf1_2, \
                         gf1_3, gg_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = pb_y[k] * gg_4[k];

        t_7[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_8[k] = f_1 * gf0_3[k]
                 - f_2 * gf1_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, fh_0, gf0_4, gf0_5, \
                         gf1_4, gf1_5, gg_6, gg_7, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * gf0_4[k]
                 - f_6 * gf1_4[k]
                 + pb_y[k] * gg_6[k];

        t_10[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pb_y[k] * gg_7[k];

        t_11[k] = pb_y[k] * gg_8[k];

        t_12[k] = f_1 * gf0_5[k]
                  - f_2 * gf1_5[k]
                  + pb_z[k] * gg_8[k];

        t_13[k] = pa_y[k] * fh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, fg_1, fg_3, fg_5, fh_0, fh_3, \
                         fh_5, fh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * fg_1[k]
                  + pa_y[k] * fh_3[k];

        t_15[k] = f_8 * fg_3[k]
                  + pa_y[k] * fh_5[k];

        t_16[k] = f_9 * fg_5[k]
                  + pa_y[k] * fh_8[k];

        t_17[k] = pa_z[k] * fh_0[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_z, pb_z, fg_0, fg_2, fg_4, fg_6, fh_4, \
                         fh_7, fh_9, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_10 * fg_0[k]
                  + pb_z[k] * gg_11[k];

        t_19[k] = f_7 * fg_2[k]
                  + pa_z[k] * fh_4[k];

        t_20[k] = f_8 * fg_4[k]
                  + pa_z[k] * fh_7[k];

        t_21[k] = f_7 * fg_6[k]
                  + pa_z[k] * fh_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pa_z, pb_z, dh0_0, dh1_0, fg_7, fg_9, \
                         fh_10, fh_12, fh_13, gg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_8 * fg_7[k]
                  + pa_z[k] * fh_10[k];

        t_23[k] = f_9 * fg_9[k]
                  + pa_z[k] * fh_12[k];

        t_24[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_y[k] * fh_13[k];

        t_25[k] = pb_z[k] * gg_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_x, pb_z, fg_15, fg_16, gf0_8, gf0_10, gf0_11, \
                         gf1_9, gf1_11, gf1_12, gg_14, gg_15, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * fg_15[k]
                  + f_5 * gf0_10[k]
                  - f_6 * gf1_11[k]
                  + pb_x[k] * gg_15[k];

        t_27[k] = f_3 * gf0_8[k]
                  - f_4 * gf1_9[k]
                  + pb_z[k] * gg_14[k];

        t_28[k] = f_7 * fg_16[k]
                  + f_3 * gf0_11[k]
                  - f_4 * gf1_12[k]
                  + pb_x[k] * gg_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pb_x, pb_z, dh0_1, dh1_14, fg_17, \
                         fh_15, gf0_9, gf1_10, gg_15, gg_16, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * gg_15[k];

        t_30[k] = f_5 * gf0_9[k]
                  - f_6 * gf1_10[k]
                  + pb_z[k] * gg_16[k];

        t_31[k] = f_7 * fg_17[k]
                  + pb_x[k] * gg_18[k];

        t_32[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_14[k]
                  + pa_x[k] * fh_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, gf0_11, gf0_12, gf0_13, gf1_12, gf1_13, \
                         gf1_14, gg_18, gg_19, gg_20, gg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * gg_18[k];

        t_34[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_12[k]
                  + pb_z[k] * gg_19[k];

        t_35[k] = f_5 * gf0_12[k]
                  - f_6 * gf1_13[k]
                  + pb_z[k] * gg_20[k];

        t_36[k] = f_1 * gf0_13[k]
                  - f_2 * gf1_14[k]
                  + pb_z[k] * gg_21[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_y, pb_z, dh0_0, dh1_0, fg_12, fh_14, \
                         gf0_14, gf1_15, gg_22, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_z[k] * fh_14[k];

        t_38[k] = pb_y[k] * gg_22[k];

        t_39[k] = f_7 * fg_12[k]
                  + pb_z[k] * gg_22[k];

        t_40[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_15[k]
                  + pb_y[k] * gg_23[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_x, pb_y, fg_19, gf0_15, gf0_16, gf1_16, gf1_17, \
                         gg_24, gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * fg_19[k]
                  + f_5 * gf0_16[k]
                  - f_6 * gf1_17[k]
                  + pb_x[k] * gg_25[k];

        t_42[k] = f_5 * gf0_15[k]
                  - f_6 * gf1_16[k]
                  + pb_y[k] * gg_24[k];

        t_43[k] = pb_y[k] * gg_25[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_y, fg_20, fg_21, gf0_17, gf0_19, gf1_18, \
                         gf1_20, gg_26, gg_27, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_7 * fg_20[k]
                  + f_3 * gf0_19[k]
                  - f_4 * gf1_20[k]
                  + pb_x[k] * gg_26[k];

        t_45[k] = f_7 * fg_21[k]
                  + pb_x[k] * gg_30[k];

        t_46[k] = f_1 * gf0_17[k]
                  - f_2 * gf1_18[k]
                  + pb_y[k] * gg_27[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_y, dh0_2, dh1_26, fh_16, gf0_18, \
                         gf0_19, gf1_19, gf1_20, gg_28, gg_29, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * gf0_18[k]
                  - f_6 * gf1_19[k]
                  + pb_y[k] * gg_28[k];

        t_48[k] = f_3 * gf0_19[k]
                  - f_4 * gf1_20[k]
                  + pb_y[k] * gg_29[k];

        t_49[k] = pb_y[k] * gg_30[k];

        t_50[k] = f_11 * dh0_2[k]
                  - f_12 * dh1_26[k]
                  + pa_x[k] * fh_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_x, fg_22, fg_24, fg_26, fg_38, \
                         fh_17, fh_18, fh_20, fh_25, fh_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_9 * fg_22[k]
                  + pa_x[k] * fh_17[k];

        t_52[k] = f_8 * fg_24[k]
                  + pa_x[k] * fh_18[k];

        t_53[k] = f_7 * fg_26[k]
                  + pa_x[k] * fh_20[k];

        t_54[k] = pa_x[k] * fh_25[k];

        t_55[k] = f_9 * fg_38[k]
                  + pa_x[k] * fh_32[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_z, fg_18, fg_41, fg_43, fh_34, \
                         fh_36, fh_44, gg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_8 * fg_18[k]
                  + pb_z[k] * gg_33[k];

        t_57[k] = f_8 * fg_41[k]
                  + pa_x[k] * fh_34[k];

        t_58[k] = f_7 * fg_43[k]
                  + pa_x[k] * fh_36[k];

        t_59[k] = pa_x[k] * fh_44[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_x, gf0_24, gf0_25, gf0_26, gf1_25, gf1_26, \
                         gf1_27, gg_35, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * gf0_24[k]
                  - f_2 * gf1_25[k]
                  + pb_x[k] * gg_35[k];

        t_61[k] = f_5 * gf0_25[k]
                  - f_6 * gf1_26[k]
                  + pb_x[k] * gg_36[k];

        t_62[k] = f_5 * gf0_26[k]
                  - f_6 * gf1_27[k]
                  + pb_x[k] * gg_37[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pb_x, gf0_27, gf0_29, gf1_28, gf1_30, \
                         gg_38, gg_39, gg_40, gg_42, gg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_28[k]
                  + pb_x[k] * gg_38[k];

        t_64[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_30[k]
                  + pb_x[k] * gg_39[k];

        t_65[k] = pb_x[k] * gg_40[k];

        t_66[k] = pb_x[k] * gg_42[k];

        t_67[k] = pb_x[k] * gg_43[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_y, pb_z, fg_28, gf0_27, gf0_28, gf1_28, \
                         gf1_29, gg_40, gg_41, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * fg_28[k]
                  + f_1 * gf0_27[k]
                  - f_2 * gf1_28[k]
                  + pb_y[k] * gg_40[k];

        t_69[k] = pb_z[k] * gg_40[k];

        t_70[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_28[k]
                  + pb_z[k] * gg_41[k];

        t_71[k] = f_5 * gf0_28[k]
                  - f_6 * gf1_29[k]
                  + pb_z[k] * gg_42[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_y, pb_z, fg_23, fg_25, fg_31, fh_19, \
                         fh_21, gf0_29, gf1_30, gg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * fg_31[k]
                  + pb_y[k] * gg_43[k];

        t_73[k] = f_1 * gf0_29[k]
                  - f_2 * gf1_30[k]
                  + pb_z[k] * gg_43[k];

        t_74[k] = f_7 * fg_23[k]
                  + pa_z[k] * fh_19[k];

        t_75[k] = f_8 * fg_25[k]
                  + pa_z[k] * fh_21[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pb_z, fg_28, fg_29, fg_30, fh_25, \
                         fh_27, fh_28, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_z[k] * fh_25[k];

        t_77[k] = f_10 * fg_28[k]
                  + pb_z[k] * gg_44[k];

        t_78[k] = f_7 * fg_29[k]
                  + pa_z[k] * fh_27[k];

        t_79[k] = f_8 * fg_30[k]
                  + pa_z[k] * fh_28[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_z, pb_x, pb_y, fg_31, fg_33, fh_29, gf0_31, \
                         gf1_33, gg_45, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_8 * fg_33[k]
                  + pb_y[k] * gg_45[k];

        t_81[k] = f_9 * fg_31[k]
                  + pa_z[k] * fh_29[k];

        t_82[k] = f_1 * gf0_31[k]
                  - f_2 * gf1_33[k]
                  + pb_x[k] * gg_46[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, gf0_32, gf0_33, gf0_34, gf1_34, gf1_35, \
                         gf1_36, gg_47, gg_48, gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_5 * gf0_32[k]
                  - f_6 * gf1_34[k]
                  + pb_x[k] * gg_47[k];

        t_84[k] = f_5 * gf0_33[k]
                  - f_6 * gf1_35[k]
                  + pb_x[k] * gg_48[k];

        t_85[k] = f_3 * gf0_34[k]
                  - f_4 * gf1_36[k]
                  + pb_x[k] * gg_49[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_z, pb_x, dh0_1, dh1_14, fh_30, \
                         gf0_36, gf1_38, gg_50, gg_51, gg_52, gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * gf0_36[k]
                  - f_4 * gf1_38[k]
                  + pb_x[k] * gg_50[k];

        t_87[k] = pb_x[k] * gg_51[k];

        t_88[k] = pb_x[k] * gg_52[k];

        t_89[k] = pb_x[k] * gg_54[k];

        t_90[k] = f_11 * dh0_1[k]
                  - f_12 * dh1_14[k]
                  + pa_z[k] * fh_30[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_y, pb_z, fg_32, fg_35, fg_36, gf0_35, gf0_36, \
                         gf1_37, gf1_38, gg_51, gg_52, gg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_7 * fg_32[k]
                  + pb_z[k] * gg_51[k];

        t_92[k] = f_7 * fg_35[k]
                  + f_5 * gf0_35[k]
                  - f_6 * gf1_37[k]
                  + pb_y[k] * gg_52[k];

        t_93[k] = f_7 * fg_36[k]
                  + f_3 * gf0_36[k]
                  - f_4 * gf1_38[k]
                  + pb_y[k] * gg_53[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pa_y, pb_y, dh0_2, dh1_26, fg_37, fg_39, \
                         fg_40, fh_31, fh_33, fh_35, gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_7 * fg_37[k]
                  + pb_y[k] * gg_54[k];

        t_95[k] = f_11 * dh0_2[k]
                  - f_12 * dh1_26[k]
                  + pa_y[k] * fh_31[k];

        t_96[k] = f_7 * fg_39[k]
                  + pa_y[k] * fh_33[k];

        t_97[k] = f_8 * fg_40[k]
                  + pa_y[k] * fh_35[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_y, pb_z, fg_34, fg_44, fg_45, fg_46, \
                         fh_40, fh_41, fh_42, gg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_9 * fg_44[k]
                  + pa_y[k] * fh_40[k];

        t_99[k] = f_8 * fg_34[k]
                  + pb_z[k] * gg_55[k];

        t_100[k] = f_8 * fg_45[k]
                   + pa_y[k] * fh_41[k];

        t_101[k] = f_7 * fg_46[k]
                   + pa_y[k] * fh_42[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_y, pb_x, pb_y, pb_z, fg_38, fg_47, \
                         fh_44, gf0_39, gf1_42, gg_56, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_10 * fg_47[k]
                   + pb_y[k] * gg_56[k];

        t_103[k] = pa_y[k] * fh_44[k];

        t_104[k] = f_1 * gf0_39[k]
                   - f_2 * gf1_42[k]
                   + pb_x[k] * gg_57[k];

        t_105[k] = f_0 * fg_38[k]
                   + pb_z[k] * gg_57[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pb_x, gf0_40, gf0_41, gf0_42, gf1_43, gf1_44, \
                         gf1_45, gg_58, gg_59, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_5 * gf0_40[k]
                   - f_6 * gf1_43[k]
                   + pb_x[k] * gg_58[k];

        t_107[k] = f_5 * gf0_41[k]
                   - f_6 * gf1_44[k]
                   + pb_x[k] * gg_59[k];

        t_108[k] = f_3 * gf0_42[k]
                   - f_4 * gf1_45[k]
                   + pb_x[k] * gg_60[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, gf0_42, gf0_44, \
                         gf1_45, gf1_47, gg_61, gg_62, gg_63, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_3 * gf0_44[k]
                   - f_4 * gf1_47[k]
                   + pb_x[k] * gg_61[k];

        t_110[k] = pb_x[k] * gg_62[k];

        t_111[k] = pb_x[k] * gg_63[k];

        t_112[k] = pb_x[k] * gg_65[k];

        t_113[k] = f_1 * gf0_42[k]
                   - f_2 * gf1_45[k]
                   + pb_y[k] * gg_62[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_y, pb_z, fg_44, gf0_43, gf0_44, \
                         gf1_46, gf1_47, gg_62, gg_63, gg_64, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * fg_44[k]
                   + pb_z[k] * gg_62[k];

        t_115[k] = f_5 * gf0_43[k]
                   - f_6 * gf1_46[k]
                   + pb_y[k] * gg_63[k];

        t_116[k] = f_3 * gf0_44[k]
                   - f_4 * gf1_47[k]
                   + pb_y[k] * gg_64[k];

        t_117[k] = pb_y[k] * gg_65[k];
    }

#pragma omp simd aligned(t_118, pb_z, fg_47, gf0_44, gf1_47, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_0 * fg_47[k]
                   + f_1 * gf0_44[k]
                   - f_2 * gf1_47[k]
                   + pb_z[k] * gg_65[k];
    }
}

auto
compute_prim_gh_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_14 = buffer.data(fg + 14);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_20 = buffer.data(gg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dh0_0, dh1_0, fg_0, fh_0, fh_1, \
                         gf0_0, gf1_0, gg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pa_y[k] * fh_0[k];

        t_2[k] = pa_z[k] * fh_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pa_y[k] * fh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_x, dh0_1, dh1_1, fg_3, fg_4, fh_3, gf0_4, \
                         gf0_5, gf1_4, gf1_5, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fg_3[k]
                 + f_6 * gf0_4[k]
                 - f_7 * gf1_4[k]
                 + pb_x[k] * gg_4[k];

        t_5[k] = f_5 * fg_4[k]
                 + f_8 * gf0_5[k]
                 - f_9 * gf1_5[k]
                 + pb_x[k] * gg_5[k];

        t_6[k] = f_3 * dh0_1[k]
                 - f_4 * dh1_1[k]
                 + pa_x[k] * fh_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_x, dh0_0, dh1_0, fg_6, fg_7, fh_2, gf0_7, \
                         gf0_8, gf1_7, gf1_8, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pa_z[k] * fh_2[k];

        t_8[k] = f_5 * fg_6[k]
                 + f_6 * gf0_7[k]
                 - f_7 * gf1_7[k]
                 + pb_x[k] * gg_8[k];

        t_9[k] = f_5 * fg_7[k]
                 + f_8 * gf0_8[k]
                 - f_9 * gf1_8[k]
                 + pb_x[k] * gg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, dh0_2, dh1_2, fg_9, fh_4, fh_5, \
                         fh_8, gf0_11, gf1_11, gg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dh0_2[k]
                  - f_4 * dh1_2[k]
                  + pa_x[k] * fh_4[k];

        t_11[k] = pa_x[k] * fh_5[k];

        t_12[k] = pa_x[k] * fh_8[k];

        t_13[k] = f_0 * fg_9[k]
                  + f_1 * gf0_11[k]
                  - f_2 * gf1_11[k]
                  + pb_y[k] * gg_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_y, dh0_1, dh1_1, fg_11, fh_5, fh_6, \
                         gf0_14, gf1_14, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * fh_5[k];

        t_15[k] = f_3 * dh0_1[k]
                  - f_4 * dh1_1[k]
                  + pa_z[k] * fh_6[k];

        t_16[k] = f_5 * fg_11[k]
                  + f_6 * gf0_14[k]
                  - f_7 * gf1_14[k]
                  + pb_y[k] * gg_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pb_y, dh0_2, dh1_2, fg_12, fh_7, fh_8, \
                         gf0_15, gf1_15, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * fg_12[k]
                  + f_8 * gf0_15[k]
                  - f_9 * gf1_15[k]
                  + pb_y[k] * gg_17[k];

        t_18[k] = f_3 * dh0_2[k]
                  - f_4 * dh1_2[k]
                  + pa_y[k] * fh_7[k];

        t_19[k] = pa_y[k] * fh_8[k];
    }

#pragma omp simd aligned(t_20, pb_z, fg_14, gf0_17, gf1_17, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * fg_14[k]
                  + f_1 * gf0_17[k]
                  - f_2 * gf1_17[k]
                  + pb_z[k] * gg_20[k];
    }
}

auto
compute_prim_gh_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_14 = buffer.data(fg + 14);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_50 = buffer.data(gf1 + 50);
    const auto *gf1_51 = buffer.data(gf1 + 51);
    const auto *gf1_62 = buffer.data(gf1 + 62);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_81 = buffer.data(gg + 81);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dh0_0, dh1_0, fg_0, fh_0, fh_1, \
                         gf0_0, gf1_0, gg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pa_y[k] * fh_0[k];

        t_2[k] = pa_z[k] * fh_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pa_y[k] * fh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_x, dh0_1, dh1_1, fg_3, fg_4, fh_3, gf0_4, \
                         gf0_5, gf1_16, gf1_17, gg_18, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fg_3[k]
                 + f_6 * gf0_4[k]
                 - f_7 * gf1_16[k]
                 + pb_x[k] * gg_18[k];

        t_5[k] = f_5 * fg_4[k]
                 + f_8 * gf0_5[k]
                 - f_9 * gf1_17[k]
                 + pb_x[k] * gg_20[k];

        t_6[k] = f_3 * dh0_1[k]
                 - f_4 * dh1_1[k]
                 + pa_x[k] * fh_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_x, dh0_0, dh1_0, fg_6, fg_7, fh_2, gf0_7, \
                         gf0_8, gf1_24, gf1_27, gg_27, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pa_z[k] * fh_2[k];

        t_8[k] = f_5 * fg_6[k]
                 + f_6 * gf0_7[k]
                 - f_7 * gf1_24[k]
                 + pb_x[k] * gg_27[k];

        t_9[k] = f_5 * fg_7[k]
                 + f_8 * gf0_8[k]
                 - f_9 * gf1_27[k]
                 + pb_x[k] * gg_28[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, dh0_2, dh1_2, fg_9, fh_4, fh_5, \
                         fh_8, gf0_11, gf1_39, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dh0_2[k]
                  - f_4 * dh1_2[k]
                  + pa_x[k] * fh_4[k];

        t_11[k] = pa_x[k] * fh_5[k];

        t_12[k] = pa_x[k] * fh_8[k];

        t_13[k] = f_0 * fg_9[k]
                  + f_1 * gf0_11[k]
                  - f_2 * gf1_39[k]
                  + pb_y[k] * gg_46[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_y, dh0_1, dh1_1, fg_11, fh_5, fh_6, \
                         gf0_14, gf1_50, gg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * fh_5[k];

        t_15[k] = f_3 * dh0_1[k]
                  - f_4 * dh1_1[k]
                  + pa_z[k] * fh_6[k];

        t_16[k] = f_5 * fg_11[k]
                  + f_6 * gf0_14[k]
                  - f_7 * gf1_50[k]
                  + pb_y[k] * gg_62[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pb_y, dh0_2, dh1_2, fg_12, fh_7, fh_8, \
                         gf0_15, gf1_51, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * fg_12[k]
                  + f_8 * gf0_15[k]
                  - f_9 * gf1_51[k]
                  + pb_y[k] * gg_63[k];

        t_18[k] = f_3 * dh0_2[k]
                  - f_4 * dh1_2[k]
                  + pa_y[k] * fh_7[k];

        t_19[k] = pa_y[k] * fh_8[k];
    }

#pragma omp simd aligned(t_20, pb_z, fg_14, gf0_17, gf1_62, gg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * fg_14[k]
                  + f_1 * gf0_17[k]
                  - f_2 * gf1_62[k]
                  + pb_z[k] * gg_81[k];
    }
}

auto
compute_prim_gh_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_29 = buffer.data(fg + 29);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_35 = buffer.data(gf0 + 35);
    const auto *gf0_37 = buffer.data(gf0 + 37);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_41 = buffer.data(gf0 + 41);
    const auto *gf0_50 = buffer.data(gf0 + 50);
    const auto *gf0_51 = buffer.data(gf0 + 51);
    const auto *gf0_56 = buffer.data(gf0 + 56);
    const auto *gf0_58 = buffer.data(gf0 + 58);
    const auto *gf0_59 = buffer.data(gf0 + 59);
    const auto *gf0_60 = buffer.data(gf0 + 60);
    const auto *gf0_61 = buffer.data(gf0 + 61);
    const auto *gf0_62 = buffer.data(gf0 + 62);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_31 = buffer.data(gf1 + 31);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_42 = buffer.data(gf1 + 42);
    const auto *gf1_43 = buffer.data(gf1 + 43);
    const auto *gf1_47 = buffer.data(gf1 + 47);
    const auto *gf1_49 = buffer.data(gf1 + 49);
    const auto *gf1_50 = buffer.data(gf1 + 50);
    const auto *gf1_51 = buffer.data(gf1 + 51);
    const auto *gf1_52 = buffer.data(gf1 + 52);
    const auto *gf1_53 = buffer.data(gf1 + 53);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf0_0, gf0_1, gf1_0, \
                         gf1_1, gg_0, gg_1, gg_2, gg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_2[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];

        t_3[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, gf0_2, gf0_3, gf0_5, gf1_2, gf1_3, gf1_5, \
                         gg_4, gg_5, gg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_5[k] = f_1 * gf0_3[k]
                 - f_2 * gf1_3[k]
                 + pb_y[k] * gg_5[k];

        t_6[k] = f_5 * gf0_5[k]
                 - f_6 * gf1_5[k]
                 + pb_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_y, pa_z, pb_y, pb_z, fh_0, gf0_6, gf1_6, \
                         gg_7, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * gf0_6[k]
                 - f_4 * gf1_6[k]
                 + pb_y[k] * gg_7[k];

        t_8[k] = f_1 * gf0_6[k]
                 - f_2 * gf1_6[k]
                 + pb_z[k] * gg_8[k];

        t_9[k] = pa_y[k] * fh_0[k];

        t_10[k] = pa_z[k] * fh_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_y, pb_x, dh0_0, dh1_0, fg_8, fg_9, fh_1, gf0_16, \
                         gf0_17, gf1_13, gf1_14, gg_12, gg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * dh0_0[k]
                  - f_8 * dh1_0[k]
                  + pa_y[k] * fh_1[k];

        t_12[k] = f_9 * fg_8[k]
                  + f_5 * gf0_16[k]
                  - f_6 * gf1_13[k]
                  + pb_x[k] * gg_12[k];

        t_13[k] = f_9 * fg_9[k]
                  + f_3 * gf0_17[k]
                  - f_4 * gf1_14[k]
                  + pb_x[k] * gg_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_x, dh0_0, dh0_1, dh1_0, dh1_1, \
                         fg_11, fh_2, fh_3, gf0_24, gf1_19, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * dh0_1[k]
                  - f_8 * dh1_1[k]
                  + pa_x[k] * fh_3[k];

        t_15[k] = f_7 * dh0_0[k]
                  - f_8 * dh1_0[k]
                  + pa_z[k] * fh_2[k];

        t_16[k] = f_9 * fg_11[k]
                  + f_5 * gf0_24[k]
                  - f_6 * gf1_19[k]
                  + pb_x[k] * gg_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pb_x, dh0_2, dh1_2, fg_12, fh_4, fh_5, \
                         fh_8, gf0_27, gf1_22, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * fg_12[k]
                  + f_3 * gf0_27[k]
                  - f_4 * gf1_22[k]
                  + pb_x[k] * gg_17[k];

        t_18[k] = f_7 * dh0_2[k]
                  - f_8 * dh1_2[k]
                  + pa_x[k] * fh_4[k];

        t_19[k] = pa_x[k] * fh_5[k];

        t_20[k] = pa_x[k] * fh_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, gf0_35, gf0_37, gf0_38, gf1_29, gf1_31, \
                         gf1_32, gg_21, gg_22, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * gf0_35[k]
                  - f_2 * gf1_29[k]
                  + pb_x[k] * gg_21[k];

        t_22[k] = f_5 * gf0_37[k]
                  - f_6 * gf1_31[k]
                  + pb_x[k] * gg_22[k];

        t_23[k] = f_5 * gf0_38[k]
                  - f_6 * gf1_32[k]
                  + pb_x[k] * gg_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_y, pb_z, fg_17, gf0_39, gf0_41, \
                         gf1_33, gf1_35, gg_24, gg_25, gg_26, gg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * gf0_39[k]
                  - f_4 * gf1_33[k]
                  + pb_x[k] * gg_24[k];

        t_25[k] = f_3 * gf0_41[k]
                  - f_4 * gf1_35[k]
                  + pb_x[k] * gg_25[k];

        t_26[k] = f_0 * fg_17[k]
                  + f_1 * gf0_39[k]
                  - f_2 * gf1_33[k]
                  + pb_y[k] * gg_26[k];

        t_27[k] = f_3 * gf0_39[k]
                  - f_4 * gf1_33[k]
                  + pb_z[k] * gg_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_z, dh0_1, dh1_1, fh_5, fh_6, gf0_40, \
                         gf0_41, gf1_34, gf1_35, gg_28, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * gf0_40[k]
                  - f_6 * gf1_34[k]
                  + pb_z[k] * gg_28[k];

        t_29[k] = f_1 * gf0_41[k]
                  - f_2 * gf1_35[k]
                  + pb_z[k] * gg_29[k];

        t_30[k] = pa_z[k] * fh_5[k];

        t_31[k] = f_7 * dh0_1[k]
                  - f_8 * dh1_1[k]
                  + pa_z[k] * fh_6[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pb_y, dh0_2, dh1_2, fg_21, fg_22, fh_7, \
                         gf0_50, gf0_51, gf1_42, gf1_43, gg_32, gg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * fg_21[k]
                  + f_5 * gf0_50[k]
                  - f_6 * gf1_42[k]
                  + pb_y[k] * gg_32[k];

        t_33[k] = f_9 * fg_22[k]
                  + f_3 * gf0_51[k]
                  - f_4 * gf1_43[k]
                  + pb_y[k] * gg_33[k];

        t_34[k] = f_7 * dh0_2[k]
                  - f_8 * dh1_2[k]
                  + pa_y[k] * fh_7[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pb_x, fh_8, gf0_56, gf0_58, gf0_59, \
                         gf1_47, gf1_49, gf1_50, gg_36, gg_37, gg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_y[k] * fh_8[k];

        t_36[k] = f_1 * gf0_56[k]
                  - f_2 * gf1_47[k]
                  + pb_x[k] * gg_36[k];

        t_37[k] = f_5 * gf0_58[k]
                  - f_6 * gf1_49[k]
                  + pb_x[k] * gg_37[k];

        t_38[k] = f_5 * gf0_59[k]
                  - f_6 * gf1_50[k]
                  + pb_x[k] * gg_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pb_x, pb_y, gf0_60, gf0_61, gf0_62, gf1_51, \
                         gf1_52, gf1_53, gg_39, gg_40, gg_41, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * gf0_60[k]
                  - f_4 * gf1_51[k]
                  + pb_x[k] * gg_39[k];

        t_40[k] = f_3 * gf0_62[k]
                  - f_4 * gf1_53[k]
                  + pb_x[k] * gg_40[k];

        t_41[k] = f_1 * gf0_60[k]
                  - f_2 * gf1_51[k]
                  + pb_y[k] * gg_41[k];

        t_42[k] = f_5 * gf0_61[k]
                  - f_6 * gf1_52[k]
                  + pb_y[k] * gg_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, fg_29, gf0_62, gf1_53, gg_43, \
                         gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gf0_62[k]
                  - f_4 * gf1_53[k]
                  + pb_y[k] * gg_43[k];

        t_44[k] = f_0 * fg_29[k]
                  + f_1 * gf0_62[k]
                  - f_2 * gf1_53[k]
                  + pb_z[k] * gg_44[k];
    }
}

auto
compute_prim_gh_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_8 = buffer.data(dh1 + 8);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_86 = buffer.data(gg + 86);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dh0_0, dh1_0, fg_0, fh_0, fh_1, \
                         gf0_0, gf1_0, gg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pa_y[k] * fh_0[k];

        t_2[k] = pa_z[k] * fh_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pa_y[k] * fh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, fg_3, fg_4, fg_5, gf0_4, gf0_5, gf1_10, gf1_11, \
                         gg_20, gg_22, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * fg_3[k]
                 + f_6 * gf0_4[k]
                 - f_7 * gf1_10[k]
                 + pb_x[k] * gg_20[k];

        t_5[k] = f_5 * fg_4[k]
                 + f_8 * gf0_5[k]
                 - f_9 * gf1_11[k]
                 + pb_x[k] * gg_22[k];

        t_6[k] = f_5 * fg_5[k]
                 + pb_x[k] * gg_23[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, dh0_0, dh0_1, dh1_0, dh1_3, fg_6, \
                         fh_2, fh_6, gf0_7, gf1_16, gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dh0_1[k]
                 - f_4 * dh1_3[k]
                 + pa_x[k] * fh_6[k];

        t_8[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pa_z[k] * fh_2[k];

        t_9[k] = f_5 * fg_6[k]
                 + f_6 * gf0_7[k]
                 - f_7 * gf1_16[k]
                 + pb_x[k] * gg_32[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, dh0_2, dh1_8, fg_7, fg_8, fh_10, \
                         fh_11, gf0_8, gf1_19, gg_33, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fg_7[k]
                  + f_8 * gf0_8[k]
                  - f_9 * gf1_19[k]
                  + pb_x[k] * gg_33[k];

        t_11[k] = f_5 * fg_8[k]
                  + pb_x[k] * gg_37[k];

        t_12[k] = f_3 * dh0_2[k]
                  - f_4 * dh1_8[k]
                  + pa_x[k] * fh_10[k];

        t_13[k] = pa_x[k] * fh_11[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, pb_y, dh0_1, dh1_3, fg_9, fh_11, \
                         fh_12, fh_17, gf0_11, gf1_27, gg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * fh_17[k];

        t_15[k] = f_0 * fg_9[k]
                  + f_1 * gf0_11[k]
                  - f_2 * gf1_27[k]
                  + pb_y[k] * gg_53[k];

        t_16[k] = pa_z[k] * fh_11[k];

        t_17[k] = f_3 * dh0_1[k]
                  - f_4 * dh1_3[k]
                  + pa_z[k] * fh_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, fg_11, fg_12, fg_13, gf0_14, gf0_15, gf1_35, \
                         gf1_36, gg_68, gg_69, gg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * fg_11[k]
                  + f_6 * gf0_14[k]
                  - f_7 * gf1_35[k]
                  + pb_y[k] * gg_68[k];

        t_19[k] = f_5 * fg_12[k]
                  + f_8 * gf0_15[k]
                  - f_9 * gf1_36[k]
                  + pb_y[k] * gg_69[k];

        t_20[k] = f_5 * fg_13[k]
                  + pb_y[k] * gg_70[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_z, dh0_2, dh1_8, fg_14, fh_16, fh_17, \
                         gf0_17, gf1_44, gg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * dh0_2[k]
                  - f_4 * dh1_8[k]
                  + pa_y[k] * fh_16[k];

        t_22[k] = pa_y[k] * fh_17[k];

        t_23[k] = f_0 * fg_14[k]
                  + f_1 * gf0_17[k]
                  - f_2 * gf1_44[k]
                  + pb_z[k] * gg_86[k];
    }
}

auto
compute_prim_gh_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dh0, const size_t dh1,
                                     const size_t fg, const size_t fh, const size_t gf0,
                                     const size_t gf1, const size_t gg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_8 = buffer.data(dh0 + 8);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_8 = buffer.data(dh1 + 8);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_35 = buffer.data(gf0 + 35);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_41 = buffer.data(gf0 + 41);
    const auto *gf0_42 = buffer.data(gf0 + 42);
    const auto *gf0_43 = buffer.data(gf0 + 43);
    const auto *gf0_44 = buffer.data(gf0 + 44);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_31 = buffer.data(gf1 + 31);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_37 = buffer.data(gf1 + 37);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_40 = buffer.data(gf1 + 40);
    const auto *gf1_44 = buffer.data(gf1 + 44);
    const auto *gf1_46 = buffer.data(gf1 + 46);
    const auto *gf1_47 = buffer.data(gf1 + 47);
    const auto *gf1_48 = buffer.data(gf1 + 48);
    const auto *gf1_49 = buffer.data(gf1 + 49);
    const auto *gf1_50 = buffer.data(gf1 + 50);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fg_0, gf0_0, gf1_0, gg_0, \
                         gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pb_y[k] * gg_0[k];

        t_2[k] = pb_z[k] * gg_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_4[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gf0_1, gf0_2, gf0_3, gf1_1, \
                         gf1_2, gf1_3, gg_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = pb_z[k] * gg_3[k];

        t_7[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_8[k] = f_1 * gf0_3[k]
                 - f_2 * gf1_3[k]
                 + pb_y[k] * gg_5[k];

        t_9[k] = pb_z[k] * gg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, fh_0, gf0_4, gf0_5, gf1_5, \
                         gf1_6, gg_7, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gf0_4[k]
                  - f_6 * gf1_5[k]
                  + pb_y[k] * gg_7[k];

        t_11[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_6[k]
                  + pb_y[k] * gg_8[k];

        t_12[k] = f_1 * gf0_5[k]
                  - f_2 * gf1_6[k]
                  + pb_z[k] * gg_9[k];

        t_13[k] = pa_y[k] * fh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, fg_1, fg_3, fg_5, fh_0, fh_1, \
                         fh_3, fh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * fg_1[k]
                  + pa_y[k] * fh_1[k];

        t_15[k] = f_8 * fg_3[k]
                  + pa_y[k] * fh_3[k];

        t_16[k] = f_9 * fg_5[k]
                  + pa_y[k] * fh_5[k];

        t_17[k] = pa_z[k] * fh_0[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_z, pb_z, fg_0, fg_2, fg_4, fg_6, fh_2, \
                         fh_4, fh_6, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_10 * fg_0[k]
                  + pb_z[k] * gg_12[k];

        t_19[k] = f_7 * fg_2[k]
                  + pa_z[k] * fh_2[k];

        t_20[k] = f_8 * fg_4[k]
                  + pa_z[k] * fh_4[k];

        t_21[k] = f_7 * fg_6[k]
                  + pa_z[k] * fh_6[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, dh0_0, dh1_0, fg_7, fg_8, fh_7, fh_8, \
                         fh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_8 * fg_7[k]
                  + pa_z[k] * fh_7[k];

        t_23[k] = f_9 * fg_8[k]
                  + pa_z[k] * fh_8[k];

        t_24[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_y[k] * fh_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, pb_z, fg_14, fg_15, gf0_8, gf0_10, gf0_11, \
                         gf1_10, gf1_12, gf1_13, gg_15, gg_16, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * fg_14[k]
                  + f_5 * gf0_10[k]
                  - f_6 * gf1_12[k]
                  + pb_x[k] * gg_16[k];

        t_26[k] = f_3 * gf0_8[k]
                  - f_4 * gf1_10[k]
                  + pb_z[k] * gg_15[k];

        t_27[k] = f_7 * fg_15[k]
                  + f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pb_x[k] * gg_18[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pb_x, pb_z, dh0_3, dh1_3, fg_16, fh_14, \
                         gf0_9, gf1_11, gg_17, gg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * gf0_9[k]
                  - f_6 * gf1_11[k]
                  + pb_z[k] * gg_17[k];

        t_29[k] = f_7 * fg_16[k]
                  + pb_x[k] * gg_19[k];

        t_30[k] = f_11 * dh0_3[k]
                  - f_12 * dh1_3[k]
                  + pa_x[k] * fh_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_z, gf0_11, gf0_12, gf0_13, gf1_13, gf1_14, \
                         gf1_15, gg_20, gg_21, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pb_z[k] * gg_20[k];

        t_32[k] = f_5 * gf0_12[k]
                  - f_6 * gf1_14[k]
                  + pb_z[k] * gg_21[k];

        t_33[k] = f_1 * gf0_13[k]
                  - f_2 * gf1_15[k]
                  + pb_z[k] * gg_22[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_y, pb_z, dh0_0, dh1_0, fg_11, fh_10, \
                         gf0_14, gf1_16, gg_23, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * dh0_0[k]
                  - f_12 * dh1_0[k]
                  + pa_z[k] * fh_10[k];

        t_35[k] = f_7 * fg_11[k]
                  + pb_z[k] * gg_23[k];

        t_36[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_16[k]
                  + pb_y[k] * gg_24[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, fg_18, fg_19, gf0_15, gf0_16, gf0_19, \
                         gf1_17, gf1_18, gf1_21, gg_25, gg_26, gg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * fg_18[k]
                  + f_5 * gf0_16[k]
                  - f_6 * gf1_18[k]
                  + pb_x[k] * gg_26[k];

        t_38[k] = f_5 * gf0_15[k]
                  - f_6 * gf1_17[k]
                  + pb_y[k] * gg_25[k];

        t_39[k] = f_7 * fg_19[k]
                  + f_3 * gf0_19[k]
                  - f_4 * gf1_21[k]
                  + pb_x[k] * gg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, pb_y, fg_20, gf0_17, gf0_18, gf1_19, gf1_20, \
                         gg_28, gg_29, gg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * fg_20[k]
                  + pb_x[k] * gg_31[k];

        t_41[k] = f_1 * gf0_17[k]
                  - f_2 * gf1_19[k]
                  + pb_y[k] * gg_28[k];

        t_42[k] = f_5 * gf0_18[k]
                  - f_6 * gf1_20[k]
                  + pb_y[k] * gg_29[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_y, dh0_8, dh1_8, fg_21, fg_23, \
                         fh_18, fh_19, fh_20, gf0_19, gf1_21, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gf0_19[k]
                  - f_4 * gf1_21[k]
                  + pb_y[k] * gg_30[k];

        t_44[k] = f_11 * dh0_8[k]
                  - f_12 * dh1_8[k]
                  + pa_x[k] * fh_18[k];

        t_45[k] = f_9 * fg_21[k]
                  + pa_x[k] * fh_19[k];

        t_46[k] = f_8 * fg_23[k]
                  + pa_x[k] * fh_20[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, fg_25, fg_26, fg_36, fh_22, \
                         fh_24, fh_33, gg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_7 * fg_25[k]
                  + pa_x[k] * fh_22[k];

        t_48[k] = f_10 * fg_26[k]
                  + pb_x[k] * gg_35[k];

        t_49[k] = pa_x[k] * fh_24[k];

        t_50[k] = f_9 * fg_36[k]
                  + pa_x[k] * fh_33[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_x, pb_x, pb_z, fg_17, fg_39, fg_40, fg_44, \
                         fh_35, fh_37, gg_36, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * fg_17[k]
                  + pb_z[k] * gg_36[k];

        t_52[k] = f_8 * fg_39[k]
                  + pa_x[k] * fh_35[k];

        t_53[k] = f_7 * fg_40[k]
                  + pa_x[k] * fh_37[k];

        t_54[k] = f_10 * fg_44[k]
                  + pb_x[k] * gg_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pb_x, pb_z, fh_41, gf0_24, gf0_25, \
                         gf1_26, gf1_28, gg_40, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_x[k] * fh_41[k];

        t_56[k] = f_1 * gf0_24[k]
                  - f_2 * gf1_26[k]
                  + pb_x[k] * gg_40[k];

        t_57[k] = pb_z[k] * gg_40[k];

        t_58[k] = f_5 * gf0_25[k]
                  - f_6 * gf1_28[k]
                  + pb_x[k] * gg_42[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pb_x, pb_z, gf0_26, gf0_27, gf0_29, gf1_29, \
                         gf1_30, gf1_32, gg_42, gg_43, gg_44, gg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * gf0_26[k]
                  - f_6 * gf1_29[k]
                  + pb_x[k] * gg_43[k];

        t_60[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_30[k]
                  + pb_x[k] * gg_44[k];

        t_61[k] = pb_z[k] * gg_42[k];

        t_62[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_32[k]
                  + pb_x[k] * gg_45[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pb_x, pb_y, pb_z, fg_26, gf0_27, \
                         gf0_28, gf1_30, gf1_31, gg_46, gg_47, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * gg_46[k];

        t_64[k] = f_0 * fg_26[k]
                  + f_1 * gf0_27[k]
                  - f_2 * gf1_30[k]
                  + pb_y[k] * gg_46[k];

        t_65[k] = pb_z[k] * gg_46[k];

        t_66[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_30[k]
                  + pb_z[k] * gg_47[k];

        t_67[k] = f_5 * gf0_28[k]
                  - f_6 * gf1_31[k]
                  + pb_z[k] * gg_48[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, pb_z, fg_22, fg_24, fg_29, fh_21, \
                         fh_23, gf0_29, gf1_32, gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * fg_29[k]
                  + pb_y[k] * gg_49[k];

        t_69[k] = f_1 * gf0_29[k]
                  - f_2 * gf1_32[k]
                  + pb_z[k] * gg_49[k];

        t_70[k] = f_7 * fg_22[k]
                  + pa_z[k] * fh_21[k];

        t_71[k] = f_8 * fg_24[k]
                  + pa_z[k] * fh_23[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_z, fg_26, fg_27, fg_28, fh_24, \
                         fh_25, fh_26, gg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * fh_24[k];

        t_73[k] = f_10 * fg_26[k]
                  + pb_z[k] * gg_50[k];

        t_74[k] = f_7 * fg_27[k]
                  + pa_z[k] * fh_25[k];

        t_75[k] = f_8 * fg_28[k]
                  + pa_z[k] * fh_26[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pb_x, pb_y, fg_29, fg_31, fh_27, gf0_31, \
                         gf1_35, gg_51, gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * fg_31[k]
                  + pb_y[k] * gg_51[k];

        t_77[k] = f_9 * fg_29[k]
                  + pa_z[k] * fh_27[k];

        t_78[k] = f_1 * gf0_31[k]
                  - f_2 * gf1_35[k]
                  + pb_x[k] * gg_52[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_x, gf0_32, gf0_33, gf0_34, gf1_36, gf1_37, \
                         gf1_38, gg_53, gg_54, gg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * gf0_32[k]
                  - f_6 * gf1_36[k]
                  + pb_x[k] * gg_53[k];

        t_80[k] = f_5 * gf0_33[k]
                  - f_6 * gf1_37[k]
                  + pb_x[k] * gg_54[k];

        t_81[k] = f_3 * gf0_34[k]
                  - f_4 * gf1_38[k]
                  + pb_x[k] * gg_55[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_z, pb_x, pb_z, dh0_3, dh1_3, fg_30, fh_28, \
                         gf0_36, gf1_40, gg_56, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * gf0_36[k]
                  - f_4 * gf1_40[k]
                  + pb_x[k] * gg_56[k];

        t_83[k] = f_11 * dh0_3[k]
                  - f_12 * dh1_3[k]
                  + pa_z[k] * fh_28[k];

        t_84[k] = f_7 * fg_30[k]
                  + pb_z[k] * gg_57[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pb_y, fg_33, fg_34, fg_35, gf0_35, gf0_36, gf1_39, \
                         gf1_40, gg_58, gg_59, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_7 * fg_33[k]
                  + f_5 * gf0_35[k]
                  - f_6 * gf1_39[k]
                  + pb_y[k] * gg_58[k];

        t_86[k] = f_7 * fg_34[k]
                  + f_3 * gf0_36[k]
                  - f_4 * gf1_40[k]
                  + pb_y[k] * gg_59[k];

        t_87[k] = f_7 * fg_35[k]
                  + pb_y[k] * gg_60[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, dh0_8, dh1_8, fg_37, fg_38, fg_41, \
                         fh_32, fh_34, fh_36, fh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_11 * dh0_8[k]
                  - f_12 * dh1_8[k]
                  + pa_y[k] * fh_32[k];

        t_89[k] = f_7 * fg_37[k]
                  + pa_y[k] * fh_34[k];

        t_90[k] = f_8 * fg_38[k]
                  + pa_y[k] * fh_36[k];

        t_91[k] = f_9 * fg_41[k]
                  + pa_y[k] * fh_38[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pb_y, pb_z, fg_32, fg_42, fg_43, fg_44, \
                         fh_39, fh_40, gg_61, gg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_8 * fg_32[k]
                  + pb_z[k] * gg_61[k];

        t_93[k] = f_8 * fg_42[k]
                  + pa_y[k] * fh_39[k];

        t_94[k] = f_7 * fg_43[k]
                  + pa_y[k] * fh_40[k];

        t_95[k] = f_10 * fg_44[k]
                  + pb_y[k] * gg_64[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pb_x, pb_y, pb_z, fg_36, fh_41, gf0_39, \
                         gf1_44, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_y[k] * fh_41[k];

        t_97[k] = f_1 * gf0_39[k]
                  - f_2 * gf1_44[k]
                  + pb_x[k] * gg_65[k];

        t_98[k] = pb_y[k] * gg_65[k];

        t_99[k] = f_0 * fg_36[k]
                  + pb_z[k] * gg_65[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_x, pb_y, gf0_40, gf0_41, gf0_42, \
                         gf1_46, gf1_47, gf1_48, gg_67, gg_68, gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * gf0_40[k]
                   - f_6 * gf1_46[k]
                   + pb_x[k] * gg_67[k];

        t_101[k] = f_5 * gf0_41[k]
                   - f_6 * gf1_47[k]
                   + pb_x[k] * gg_68[k];

        t_102[k] = f_3 * gf0_42[k]
                   - f_4 * gf1_48[k]
                   + pb_x[k] * gg_69[k];

        t_103[k] = pb_y[k] * gg_68[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pb_y, pb_z, fg_41, gf0_42, gf0_44, \
                         gf1_48, gf1_50, gg_70, gg_71, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * gf0_44[k]
                   - f_4 * gf1_50[k]
                   + pb_x[k] * gg_70[k];

        t_105[k] = pb_x[k] * gg_74[k];

        t_106[k] = f_1 * gf0_42[k]
                   - f_2 * gf1_48[k]
                   + pb_y[k] * gg_71[k];

        t_107[k] = f_0 * fg_41[k]
                   + pb_z[k] * gg_71[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_y, pb_z, fg_44, gf0_43, gf0_44, \
                         gf1_49, gf1_50, gg_72, gg_73, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * gf0_43[k]
                   - f_6 * gf1_49[k]
                   + pb_y[k] * gg_72[k];

        t_109[k] = f_3 * gf0_44[k]
                   - f_4 * gf1_50[k]
                   + pb_y[k] * gg_73[k];

        t_110[k] = pb_y[k] * gg_74[k];

        t_111[k] = f_0 * fg_44[k]
                   + f_1 * gf0_44[k]
                   - f_2 * gf1_50[k]
                   + pb_z[k] * gg_74[k];
    }
}

auto
compute_prim_gh_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dh0, const size_t dh1,
                                      const size_t fg, const size_t fh, const size_t gf0,
                                      const size_t gf1, const size_t gg, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / p;

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

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_8 = buffer.data(dh0 + 8);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_8 = buffer.data(dh1 + 8);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_38 = buffer.data(fg + 38);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_29 = buffer.data(gf0 + 29);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_31 = buffer.data(gf0 + 31);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_44 = buffer.data(gf0 + 44);
    const auto *gf0_46 = buffer.data(gf0 + 46);
    const auto *gf0_47 = buffer.data(gf0 + 47);
    const auto *gf0_48 = buffer.data(gf0 + 48);
    const auto *gf0_49 = buffer.data(gf0 + 49);
    const auto *gf0_50 = buffer.data(gf0 + 50);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_40 = buffer.data(gf1 + 40);
    const auto *gf1_41 = buffer.data(gf1 + 41);
    const auto *gf1_42 = buffer.data(gf1 + 42);
    const auto *gf1_43 = buffer.data(gf1 + 43);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fg_0, gf0_0, gf1_0, gg_0, \
                         gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 + f_1 * gf0_0[k]
                 - f_2 * gf1_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = pb_y[k] * gg_0[k];

        t_2[k] = pb_z[k] * gg_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_y[k] * gg_1[k];

        t_4[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gf0_1, gf0_2, gf0_3, gf1_1, gf1_2, \
                         gf1_3, gg_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * gf0_1[k]
                 - f_6 * gf1_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = pb_y[k] * gg_4[k];

        t_7[k] = f_5 * gf0_2[k]
                 - f_6 * gf1_2[k]
                 + pb_z[k] * gg_4[k];

        t_8[k] = f_1 * gf0_3[k]
                 - f_2 * gf1_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, fh_0, gf0_5, gf0_6, \
                         gf1_4, gf1_5, gg_6, gg_7, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * gf0_5[k]
                 - f_6 * gf1_4[k]
                 + pb_y[k] * gg_6[k];

        t_10[k] = f_3 * gf0_6[k]
                  - f_4 * gf1_5[k]
                  + pb_y[k] * gg_7[k];

        t_11[k] = pb_y[k] * gg_8[k];

        t_12[k] = f_1 * gf0_6[k]
                  - f_2 * gf1_5[k]
                  + pb_z[k] * gg_8[k];

        t_13[k] = pa_y[k] * fh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_y, pa_z, pb_x, dh0_0, dh1_0, fg_11, fh_0, fh_1, \
                         gf0_12, gf1_10, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * fh_0[k];

        t_15[k] = f_7 * dh0_0[k]
                  - f_8 * dh1_0[k]
                  + pa_y[k] * fh_1[k];

        t_16[k] = f_9 * fg_11[k]
                  + f_5 * gf0_12[k]
                  - f_6 * gf1_10[k]
                  + pb_x[k] * gg_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pb_x, dh0_3, dh1_3, fg_12, fg_13, fh_3, \
                         gf0_13, gf1_11, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * fg_12[k]
                  + f_3 * gf0_13[k]
                  - f_4 * gf1_11[k]
                  + pb_x[k] * gg_13[k];

        t_18[k] = f_9 * fg_13[k]
                  + pb_x[k] * gg_14[k];

        t_19[k] = f_7 * dh0_3[k]
                  - f_8 * dh1_3[k]
                  + pa_x[k] * fh_3[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_x, dh0_0, dh1_0, fg_14, fg_15, fh_2, \
                         gf0_18, gf0_21, gf1_16, gf1_19, gg_16, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_7 * dh0_0[k]
                  - f_8 * dh1_0[k]
                  + pa_z[k] * fh_2[k];

        t_21[k] = f_9 * fg_14[k]
                  + f_5 * gf0_18[k]
                  - f_6 * gf1_16[k]
                  + pb_x[k] * gg_16[k];

        t_22[k] = f_9 * fg_15[k]
                  + f_3 * gf0_21[k]
                  - f_4 * gf1_19[k]
                  + pb_x[k] * gg_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, dh0_8, dh1_8, fg_16, fh_4, fh_5, \
                         fh_8, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_9 * fg_16[k]
                  + pb_x[k] * gg_18[k];

        t_24[k] = f_7 * dh0_8[k]
                  - f_8 * dh1_8[k]
                  + pa_x[k] * fh_4[k];

        t_25[k] = pa_x[k] * fh_5[k];

        t_26[k] = pa_x[k] * fh_8[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, gf0_26, gf0_28, gf0_29, gf1_24, gf1_25, \
                         gf1_26, gg_21, gg_22, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * gf0_26[k]
                  - f_2 * gf1_24[k]
                  + pb_x[k] * gg_21[k];

        t_28[k] = f_5 * gf0_28[k]
                  - f_6 * gf1_25[k]
                  + pb_x[k] * gg_22[k];

        t_29[k] = f_5 * gf0_29[k]
                  - f_6 * gf1_26[k]
                  + pb_x[k] * gg_23[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, gf0_30, gf0_32, gf1_27, gf1_29, \
                         gg_24, gg_25, gg_26, gg_28, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * gf0_30[k]
                  - f_4 * gf1_27[k]
                  + pb_x[k] * gg_24[k];

        t_31[k] = f_3 * gf0_32[k]
                  - f_4 * gf1_29[k]
                  + pb_x[k] * gg_25[k];

        t_32[k] = pb_x[k] * gg_26[k];

        t_33[k] = pb_x[k] * gg_28[k];

        t_34[k] = pb_x[k] * gg_29[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_y, pb_z, fg_22, gf0_30, gf0_31, gf1_27, \
                         gf1_28, gg_26, gg_27, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * fg_22[k]
                  + f_1 * gf0_30[k]
                  - f_2 * gf1_27[k]
                  + pb_y[k] * gg_26[k];

        t_36[k] = pb_z[k] * gg_26[k];

        t_37[k] = f_3 * gf0_30[k]
                  - f_4 * gf1_27[k]
                  + pb_z[k] * gg_27[k];

        t_38[k] = f_5 * gf0_31[k]
                  - f_6 * gf1_28[k]
                  + pb_z[k] * gg_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_z, pb_z, dh0_3, dh1_3, fh_5, fh_6, gf0_32, \
                         gf1_29, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * gf0_32[k]
                  - f_2 * gf1_29[k]
                  + pb_z[k] * gg_29[k];

        t_40[k] = pa_z[k] * fh_5[k];

        t_41[k] = f_7 * dh0_3[k]
                  - f_8 * dh1_3[k]
                  + pa_z[k] * fh_6[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_y, fg_27, fg_28, fg_29, gf0_39, gf0_40, gf1_35, \
                         gf1_36, gg_32, gg_33, gg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_9 * fg_27[k]
                  + f_5 * gf0_39[k]
                  - f_6 * gf1_35[k]
                  + pb_y[k] * gg_32[k];

        t_43[k] = f_9 * fg_28[k]
                  + f_3 * gf0_40[k]
                  - f_4 * gf1_36[k]
                  + pb_y[k] * gg_33[k];

        t_44[k] = f_9 * fg_29[k]
                  + pb_y[k] * gg_34[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pb_x, dh0_8, dh1_8, fh_7, fh_8, gf0_44, \
                         gf0_46, gf1_39, gf1_40, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_7 * dh0_8[k]
                  - f_8 * dh1_8[k]
                  + pa_y[k] * fh_7[k];

        t_46[k] = pa_y[k] * fh_8[k];

        t_47[k] = f_1 * gf0_44[k]
                  - f_2 * gf1_39[k]
                  + pb_x[k] * gg_36[k];

        t_48[k] = f_5 * gf0_46[k]
                  - f_6 * gf1_40[k]
                  + pb_x[k] * gg_37[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_x, gf0_47, gf0_48, gf0_50, gf1_41, gf1_42, \
                         gf1_44, gg_38, gg_39, gg_40, gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * gf0_47[k]
                  - f_6 * gf1_41[k]
                  + pb_x[k] * gg_38[k];

        t_50[k] = f_3 * gf0_48[k]
                  - f_4 * gf1_42[k]
                  + pb_x[k] * gg_39[k];

        t_51[k] = f_3 * gf0_50[k]
                  - f_4 * gf1_44[k]
                  + pb_x[k] * gg_40[k];

        t_52[k] = pb_x[k] * gg_41[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, gf0_48, gf0_49, gf1_42, gf1_43, \
                         gg_41, gg_42, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_x[k] * gg_42[k];

        t_54[k] = pb_x[k] * gg_44[k];

        t_55[k] = f_1 * gf0_48[k]
                  - f_2 * gf1_42[k]
                  + pb_y[k] * gg_41[k];

        t_56[k] = f_5 * gf0_49[k]
                  - f_6 * gf1_43[k]
                  + pb_y[k] * gg_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_y, pb_z, fg_38, gf0_50, gf1_44, gg_43, \
                         gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * gf0_50[k]
                  - f_4 * gf1_44[k]
                  + pb_y[k] * gg_43[k];

        t_58[k] = pb_y[k] * gg_44[k];

        t_59[k] = f_0 * fg_38[k]
                  + f_1 * gf0_50[k]
                  - f_2 * gf1_44[k]
                  + pb_z[k] * gg_44[k];
    }
}

}  // namespace simdt2ceri
