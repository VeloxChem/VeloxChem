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


#include "SimdKineticEnergyVrrRecGH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_gh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.5 / p;

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

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_24 = buffer.data(dh_s + 24);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_24 = buffer.data(dh + 24);

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
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_3 = buffer.data(gf_s + 3);
    const auto *gf_s_4 = buffer.data(gf_s + 4);
    const auto *gf_s_5 = buffer.data(gf_s + 5);
    const auto *gf_s_7 = buffer.data(gf_s + 7);
    const auto *gf_s_8 = buffer.data(gf_s + 8);
    const auto *gf_s_12 = buffer.data(gf_s + 12);
    const auto *gf_s_13 = buffer.data(gf_s + 13);
    const auto *gf_s_14 = buffer.data(gf_s + 14);
    const auto *gf_s_15 = buffer.data(gf_s + 15);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_17 = buffer.data(gf_s + 17);
    const auto *gf_s_18 = buffer.data(gf_s + 18);
    const auto *gf_s_19 = buffer.data(gf_s + 19);
    const auto *gf_s_20 = buffer.data(gf_s + 20);
    const auto *gf_s_23 = buffer.data(gf_s + 23);
    const auto *gf_s_24 = buffer.data(gf_s + 24);
    const auto *gf_s_25 = buffer.data(gf_s + 25);
    const auto *gf_s_26 = buffer.data(gf_s + 26);
    const auto *gf_s_27 = buffer.data(gf_s + 27);
    const auto *gf_s_28 = buffer.data(gf_s + 28);
    const auto *gf_s_29 = buffer.data(gf_s + 29);
    const auto *gf_s_30 = buffer.data(gf_s + 30);
    const auto *gf_s_39 = buffer.data(gf_s + 39);
    const auto *gf_s_40 = buffer.data(gf_s + 40);
    const auto *gf_s_41 = buffer.data(gf_s + 41);
    const auto *gf_s_42 = buffer.data(gf_s + 42);
    const auto *gf_s_43 = buffer.data(gf_s + 43);
    const auto *gf_s_44 = buffer.data(gf_s + 44);
    const auto *gf_s_45 = buffer.data(gf_s + 45);
    const auto *gf_s_46 = buffer.data(gf_s + 46);
    const auto *gf_s_47 = buffer.data(gf_s + 47);
    const auto *gf_s_49 = buffer.data(gf_s + 49);
    const auto *gf_s_53 = buffer.data(gf_s + 53);
    const auto *gf_s_54 = buffer.data(gf_s + 54);
    const auto *gf_s_55 = buffer.data(gf_s + 55);
    const auto *gf_s_56 = buffer.data(gf_s + 56);
    const auto *gf_s_57 = buffer.data(gf_s + 57);
    const auto *gf_s_58 = buffer.data(gf_s + 58);
    const auto *gf_s_59 = buffer.data(gf_s + 59);
    const auto *gf_s_60 = buffer.data(gf_s + 60);
    const auto *gf_s_61 = buffer.data(gf_s + 61);
    const auto *gf_s_62 = buffer.data(gf_s + 62);
    const auto *gf_s_63 = buffer.data(gf_s + 63);
    const auto *gf_s_71 = buffer.data(gf_s + 71);
    const auto *gf_s_72 = buffer.data(gf_s + 72);
    const auto *gf_s_73 = buffer.data(gf_s + 73);
    const auto *gf_s_74 = buffer.data(gf_s + 74);
    const auto *gf_s_75 = buffer.data(gf_s + 75);
    const auto *gf_s_76 = buffer.data(gf_s + 76);
    const auto *gf_s_77 = buffer.data(gf_s + 77);
    const auto *gf_s_78 = buffer.data(gf_s + 78);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_20 = buffer.data(gh_s + 20);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_31 = buffer.data(gh_s + 31);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_33 = buffer.data(gh_s + 33);
    const auto *gh_s_34 = buffer.data(gh_s + 34);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);
    const auto *gh_s_58 = buffer.data(gh_s + 58);
    const auto *gh_s_59 = buffer.data(gh_s + 59);
    const auto *gh_s_60 = buffer.data(gh_s + 60);
    const auto *gh_s_61 = buffer.data(gh_s + 61);
    const auto *gh_s_62 = buffer.data(gh_s + 62);
    const auto *gh_s_63 = buffer.data(gh_s + 63);
    const auto *gh_s_64 = buffer.data(gh_s + 64);
    const auto *gh_s_65 = buffer.data(gh_s + 65);
    const auto *gh_s_66 = buffer.data(gh_s + 66);
    const auto *gh_s_67 = buffer.data(gh_s + 67);
    const auto *gh_s_68 = buffer.data(gh_s + 68);
    const auto *gh_s_69 = buffer.data(gh_s + 69);
    const auto *gh_s_70 = buffer.data(gh_s + 70);
    const auto *gh_s_71 = buffer.data(gh_s + 71);
    const auto *gh_s_72 = buffer.data(gh_s + 72);
    const auto *gh_s_73 = buffer.data(gh_s + 73);
    const auto *gh_s_74 = buffer.data(gh_s + 74);
    const auto *gh_s_75 = buffer.data(gh_s + 75);
    const auto *gh_s_76 = buffer.data(gh_s + 76);
    const auto *gh_s_77 = buffer.data(gh_s + 77);
    const auto *gh_s_78 = buffer.data(gh_s + 78);
    const auto *gh_s_79 = buffer.data(gh_s + 79);
    const auto *gh_s_80 = buffer.data(gh_s + 80);
    const auto *gh_s_81 = buffer.data(gh_s + 81);
    const auto *gh_s_82 = buffer.data(gh_s + 82);
    const auto *gh_s_83 = buffer.data(gh_s + 83);
    const auto *gh_s_84 = buffer.data(gh_s + 84);
    const auto *gh_s_85 = buffer.data(gh_s + 85);
    const auto *gh_s_86 = buffer.data(gh_s + 86);
    const auto *gh_s_87 = buffer.data(gh_s + 87);
    const auto *gh_s_88 = buffer.data(gh_s + 88);
    const auto *gh_s_89 = buffer.data(gh_s + 89);
    const auto *gh_s_90 = buffer.data(gh_s + 90);
    const auto *gh_s_91 = buffer.data(gh_s + 91);
    const auto *gh_s_92 = buffer.data(gh_s + 92);
    const auto *gh_s_93 = buffer.data(gh_s + 93);
    const auto *gh_s_94 = buffer.data(gh_s + 94);
    const auto *gh_s_95 = buffer.data(gh_s + 95);
    const auto *gh_s_96 = buffer.data(gh_s + 96);
    const auto *gh_s_97 = buffer.data(gh_s + 97);
    const auto *gh_s_98 = buffer.data(gh_s + 98);
    const auto *gh_s_99 = buffer.data(gh_s + 99);
    const auto *gh_s_100 = buffer.data(gh_s + 100);
    const auto *gh_s_101 = buffer.data(gh_s + 101);
    const auto *gh_s_102 = buffer.data(gh_s + 102);
    const auto *gh_s_103 = buffer.data(gh_s + 103);
    const auto *gh_s_104 = buffer.data(gh_s + 104);
    const auto *gh_s_105 = buffer.data(gh_s + 105);
    const auto *gh_s_106 = buffer.data(gh_s + 106);
    const auto *gh_s_107 = buffer.data(gh_s + 107);
    const auto *gh_s_108 = buffer.data(gh_s + 108);
    const auto *gh_s_109 = buffer.data(gh_s + 109);
    const auto *gh_s_110 = buffer.data(gh_s + 110);
    const auto *gh_s_111 = buffer.data(gh_s + 111);
    const auto *gh_s_112 = buffer.data(gh_s + 112);
    const auto *gh_s_113 = buffer.data(gh_s + 113);
    const auto *gh_s_114 = buffer.data(gh_s + 114);
    const auto *gh_s_115 = buffer.data(gh_s + 115);
    const auto *gh_s_116 = buffer.data(gh_s + 116);
    const auto *gh_s_117 = buffer.data(gh_s + 117);
    const auto *gh_s_118 = buffer.data(gh_s + 118);
    const auto *gh_s_119 = buffer.data(gh_s + 119);
    const auto *gh_s_120 = buffer.data(gh_s + 120);
    const auto *gh_s_121 = buffer.data(gh_s + 121);
    const auto *gh_s_122 = buffer.data(gh_s + 122);
    const auto *gh_s_123 = buffer.data(gh_s + 123);
    const auto *gh_s_124 = buffer.data(gh_s + 124);
    const auto *gh_s_125 = buffer.data(gh_s + 125);
    const auto *gh_s_126 = buffer.data(gh_s + 126);
    const auto *gh_s_127 = buffer.data(gh_s + 127);
    const auto *gh_s_128 = buffer.data(gh_s + 128);
    const auto *gh_s_129 = buffer.data(gh_s + 129);
    const auto *gh_s_130 = buffer.data(gh_s + 130);
    const auto *gh_s_131 = buffer.data(gh_s + 131);
    const auto *gh_s_132 = buffer.data(gh_s + 132);
    const auto *gh_s_133 = buffer.data(gh_s + 133);
    const auto *gh_s_134 = buffer.data(gh_s + 134);
    const auto *gh_s_135 = buffer.data(gh_s + 135);
    const auto *gh_s_136 = buffer.data(gh_s + 136);
    const auto *gh_s_137 = buffer.data(gh_s + 137);
    const auto *gh_s_138 = buffer.data(gh_s + 138);
    const auto *gh_s_139 = buffer.data(gh_s + 139);
    const auto *gh_s_140 = buffer.data(gh_s + 140);
    const auto *gh_s_141 = buffer.data(gh_s + 141);
    const auto *gh_s_142 = buffer.data(gh_s + 142);
    const auto *gh_s_143 = buffer.data(gh_s + 143);
    const auto *gh_s_144 = buffer.data(gh_s + 144);
    const auto *gh_s_145 = buffer.data(gh_s + 145);
    const auto *gh_s_146 = buffer.data(gh_s + 146);
    const auto *gh_s_147 = buffer.data(gh_s + 147);
    const auto *gh_s_148 = buffer.data(gh_s + 148);
    const auto *gh_s_149 = buffer.data(gh_s + 149);
    const auto *gh_s_150 = buffer.data(gh_s + 150);
    const auto *gh_s_151 = buffer.data(gh_s + 151);
    const auto *gh_s_152 = buffer.data(gh_s + 152);
    const auto *gh_s_153 = buffer.data(gh_s + 153);
    const auto *gh_s_154 = buffer.data(gh_s + 154);
    const auto *gh_s_155 = buffer.data(gh_s + 155);
    const auto *gh_s_156 = buffer.data(gh_s + 156);
    const auto *gh_s_157 = buffer.data(gh_s + 157);
    const auto *gh_s_158 = buffer.data(gh_s + 158);
    const auto *gh_s_159 = buffer.data(gh_s + 159);
    const auto *gh_s_160 = buffer.data(gh_s + 160);
    const auto *gh_s_161 = buffer.data(gh_s + 161);
    const auto *gh_s_162 = buffer.data(gh_s + 162);
    const auto *gh_s_163 = buffer.data(gh_s + 163);
    const auto *gh_s_164 = buffer.data(gh_s + 164);
    const auto *gh_s_165 = buffer.data(gh_s + 165);
    const auto *gh_s_166 = buffer.data(gh_s + 166);
    const auto *gh_s_167 = buffer.data(gh_s + 167);
    const auto *gh_s_168 = buffer.data(gh_s + 168);
    const auto *gh_s_169 = buffer.data(gh_s + 169);
    const auto *gh_s_170 = buffer.data(gh_s + 170);
    const auto *gh_s_171 = buffer.data(gh_s + 171);
    const auto *gh_s_172 = buffer.data(gh_s + 172);
    const auto *gh_s_173 = buffer.data(gh_s + 173);
    const auto *gh_s_174 = buffer.data(gh_s + 174);
    const auto *gh_s_175 = buffer.data(gh_s + 175);
    const auto *gh_s_176 = buffer.data(gh_s + 176);
    const auto *gh_s_177 = buffer.data(gh_s + 177);
    const auto *gh_s_178 = buffer.data(gh_s + 178);
    const auto *gh_s_179 = buffer.data(gh_s + 179);
    const auto *gh_s_180 = buffer.data(gh_s + 180);
    const auto *gh_s_181 = buffer.data(gh_s + 181);
    const auto *gh_s_182 = buffer.data(gh_s + 182);
    const auto *gh_s_183 = buffer.data(gh_s + 183);
    const auto *gh_s_184 = buffer.data(gh_s + 184);
    const auto *gh_s_185 = buffer.data(gh_s + 185);
    const auto *gh_s_186 = buffer.data(gh_s + 186);
    const auto *gh_s_187 = buffer.data(gh_s + 187);
    const auto *gh_s_188 = buffer.data(gh_s + 188);
    const auto *gh_s_189 = buffer.data(gh_s + 189);
    const auto *gh_s_190 = buffer.data(gh_s + 190);
    const auto *gh_s_191 = buffer.data(gh_s + 191);
    const auto *gh_s_192 = buffer.data(gh_s + 192);
    const auto *gh_s_193 = buffer.data(gh_s + 193);
    const auto *gh_s_194 = buffer.data(gh_s + 194);
    const auto *gh_s_195 = buffer.data(gh_s + 195);
    const auto *gh_s_196 = buffer.data(gh_s + 196);
    const auto *gh_s_197 = buffer.data(gh_s + 197);
    const auto *gh_s_198 = buffer.data(gh_s + 198);
    const auto *gh_s_199 = buffer.data(gh_s + 199);
    const auto *gh_s_200 = buffer.data(gh_s + 200);
    const auto *gh_s_201 = buffer.data(gh_s + 201);
    const auto *gh_s_202 = buffer.data(gh_s + 202);
    const auto *gh_s_203 = buffer.data(gh_s + 203);
    const auto *gh_s_204 = buffer.data(gh_s + 204);
    const auto *gh_s_205 = buffer.data(gh_s + 205);
    const auto *gh_s_206 = buffer.data(gh_s + 206);
    const auto *gh_s_207 = buffer.data(gh_s + 207);
    const auto *gh_s_208 = buffer.data(gh_s + 208);
    const auto *gh_s_209 = buffer.data(gh_s + 209);
    const auto *gh_s_210 = buffer.data(gh_s + 210);
    const auto *gh_s_211 = buffer.data(gh_s + 211);
    const auto *gh_s_212 = buffer.data(gh_s + 212);
    const auto *gh_s_213 = buffer.data(gh_s + 213);
    const auto *gh_s_214 = buffer.data(gh_s + 214);
    const auto *gh_s_215 = buffer.data(gh_s + 215);
    const auto *gh_s_216 = buffer.data(gh_s + 216);
    const auto *gh_s_217 = buffer.data(gh_s + 217);
    const auto *gh_s_218 = buffer.data(gh_s + 218);
    const auto *gh_s_219 = buffer.data(gh_s + 219);
    const auto *gh_s_220 = buffer.data(gh_s + 220);
    const auto *gh_s_221 = buffer.data(gh_s + 221);
    const auto *gh_s_222 = buffer.data(gh_s + 222);
    const auto *gh_s_223 = buffer.data(gh_s + 223);
    const auto *gh_s_224 = buffer.data(gh_s + 224);
    const auto *gh_s_225 = buffer.data(gh_s + 225);
    const auto *gh_s_226 = buffer.data(gh_s + 226);
    const auto *gh_s_227 = buffer.data(gh_s + 227);
    const auto *gh_s_228 = buffer.data(gh_s + 228);
    const auto *gh_s_229 = buffer.data(gh_s + 229);
    const auto *gh_s_230 = buffer.data(gh_s + 230);
    const auto *gh_s_231 = buffer.data(gh_s + 231);
    const auto *gh_s_232 = buffer.data(gh_s + 232);
    const auto *gh_s_233 = buffer.data(gh_s + 233);
    const auto *gh_s_234 = buffer.data(gh_s + 234);
    const auto *gh_s_235 = buffer.data(gh_s + 235);
    const auto *gh_s_236 = buffer.data(gh_s + 236);
    const auto *gh_s_237 = buffer.data(gh_s + 237);
    const auto *gh_s_238 = buffer.data(gh_s + 238);
    const auto *gh_s_239 = buffer.data(gh_s + 239);
    const auto *gh_s_240 = buffer.data(gh_s + 240);
    const auto *gh_s_241 = buffer.data(gh_s + 241);
    const auto *gh_s_242 = buffer.data(gh_s + 242);
    const auto *gh_s_243 = buffer.data(gh_s + 243);
    const auto *gh_s_244 = buffer.data(gh_s + 244);
    const auto *gh_s_245 = buffer.data(gh_s + 245);
    const auto *gh_s_246 = buffer.data(gh_s + 246);
    const auto *gh_s_247 = buffer.data(gh_s + 247);
    const auto *gh_s_248 = buffer.data(gh_s + 248);
    const auto *gh_s_249 = buffer.data(gh_s + 249);
    const auto *gh_s_250 = buffer.data(gh_s + 250);
    const auto *gh_s_251 = buffer.data(gh_s + 251);
    const auto *gh_s_252 = buffer.data(gh_s + 252);
    const auto *gh_s_253 = buffer.data(gh_s + 253);
    const auto *gh_s_254 = buffer.data(gh_s + 254);
    const auto *gh_s_255 = buffer.data(gh_s + 255);
    const auto *gh_s_256 = buffer.data(gh_s + 256);
    const auto *gh_s_257 = buffer.data(gh_s + 257);
    const auto *gh_s_258 = buffer.data(gh_s + 258);
    const auto *gh_s_259 = buffer.data(gh_s + 259);
    const auto *gh_s_260 = buffer.data(gh_s + 260);
    const auto *gh_s_261 = buffer.data(gh_s + 261);
    const auto *gh_s_262 = buffer.data(gh_s + 262);
    const auto *gh_s_263 = buffer.data(gh_s + 263);
    const auto *gh_s_264 = buffer.data(gh_s + 264);
    const auto *gh_s_265 = buffer.data(gh_s + 265);
    const auto *gh_s_266 = buffer.data(gh_s + 266);
    const auto *gh_s_267 = buffer.data(gh_s + 267);
    const auto *gh_s_268 = buffer.data(gh_s + 268);
    const auto *gh_s_269 = buffer.data(gh_s + 269);
    const auto *gh_s_270 = buffer.data(gh_s + 270);
    const auto *gh_s_271 = buffer.data(gh_s + 271);
    const auto *gh_s_272 = buffer.data(gh_s + 272);
    const auto *gh_s_273 = buffer.data(gh_s + 273);
    const auto *gh_s_274 = buffer.data(gh_s + 274);
    const auto *gh_s_275 = buffer.data(gh_s + 275);
    const auto *gh_s_276 = buffer.data(gh_s + 276);
    const auto *gh_s_277 = buffer.data(gh_s + 277);
    const auto *gh_s_278 = buffer.data(gh_s + 278);
    const auto *gh_s_279 = buffer.data(gh_s + 279);
    const auto *gh_s_280 = buffer.data(gh_s + 280);
    const auto *gh_s_281 = buffer.data(gh_s + 281);
    const auto *gh_s_282 = buffer.data(gh_s + 282);
    const auto *gh_s_283 = buffer.data(gh_s + 283);
    const auto *gh_s_284 = buffer.data(gh_s + 284);
    const auto *gh_s_285 = buffer.data(gh_s + 285);
    const auto *gh_s_286 = buffer.data(gh_s + 286);
    const auto *gh_s_287 = buffer.data(gh_s + 287);
    const auto *gh_s_288 = buffer.data(gh_s + 288);
    const auto *gh_s_289 = buffer.data(gh_s + 289);
    const auto *gh_s_290 = buffer.data(gh_s + 290);
    const auto *gh_s_291 = buffer.data(gh_s + 291);
    const auto *gh_s_292 = buffer.data(gh_s + 292);
    const auto *gh_s_293 = buffer.data(gh_s + 293);
    const auto *gh_s_294 = buffer.data(gh_s + 294);
    const auto *gh_s_295 = buffer.data(gh_s + 295);
    const auto *gh_s_296 = buffer.data(gh_s + 296);
    const auto *gh_s_297 = buffer.data(gh_s + 297);
    const auto *gh_s_298 = buffer.data(gh_s + 298);
    const auto *gh_s_299 = buffer.data(gh_s + 299);
    const auto *gh_s_300 = buffer.data(gh_s + 300);
    const auto *gh_s_301 = buffer.data(gh_s + 301);
    const auto *gh_s_302 = buffer.data(gh_s + 302);
    const auto *gh_s_303 = buffer.data(gh_s + 303);
    const auto *gh_s_304 = buffer.data(gh_s + 304);
    const auto *gh_s_305 = buffer.data(gh_s + 305);
    const auto *gh_s_306 = buffer.data(gh_s + 306);
    const auto *gh_s_307 = buffer.data(gh_s + 307);
    const auto *gh_s_308 = buffer.data(gh_s + 308);
    const auto *gh_s_309 = buffer.data(gh_s + 309);
    const auto *gh_s_310 = buffer.data(gh_s + 310);
    const auto *gh_s_311 = buffer.data(gh_s + 311);
    const auto *gh_s_312 = buffer.data(gh_s + 312);
    const auto *gh_s_313 = buffer.data(gh_s + 313);
    const auto *gh_s_314 = buffer.data(gh_s + 314);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);

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
    const auto *gg_144 = buffer.data(gg + 144);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gh_s_3, gf_0, gg_0, gg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_2 * gh_s_1[k]
                 + pb_y[k] * gg_0[k];

        t_2[k] = f_2 * gh_s_2[k]
                 + pb_z[k] * gg_0[k];

        t_3[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_3[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, gf_s_0, gf_s_1, gh_s_4, gh_s_5, \
                         gh_s_6, gh_s_7, gf_0, gf_1, gg_2, gg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gh_s_4[k]
                 + pb_y[k] * gg_2[k];

        t_5[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_5[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];

        t_6[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_6[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_7[k] = f_2 * gh_s_7[k]
                 + pb_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, fg_5, gf_s_2, gh_s_8, gh_s_9, \
                         gh_s_10, gf_2, gg_4, gg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * gh_s_8[k]
                 + pb_y[k] * gg_4[k];

        t_9[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_9[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_4[k];

        t_10[k] = f_0 * fg_5[k]
                  + f_2 * gh_s_10[k]
                  + pb_x[k] * gg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, fg_6, gh_s_11, gh_s_12, gh_s_13, \
                         gg_5, gg_6, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * gh_s_11[k]
                  + pb_z[k] * gg_5[k];

        t_12[k] = f_0 * fg_6[k]
                  + f_2 * gh_s_12[k]
                  + pb_x[k] * gg_8[k];

        t_13[k] = f_2 * gh_s_13[k]
                  + pb_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, fg_7, gf_s_3, gh_s_14, gh_s_15, \
                         gh_s_16, gf_3, gg_7, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * fg_7[k]
                  + f_2 * gh_s_14[k]
                  + pb_x[k] * gg_10[k];

        t_15[k] = -f_1 * gf_s_3[k]
                  + f_2 * gh_s_15[k]
                  + f_0 * gf_3[k]
                  + pb_y[k] * gg_7[k];

        t_16[k] = f_2 * gh_s_16[k]
                  + pb_z[k] * gg_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, gf_s_4, gf_s_5, gh_s_17, gh_s_18, gh_s_19, \
                         gf_4, gf_5, gg_8, gg_9, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_5 * gf_s_4[k]
                  + f_2 * gh_s_17[k]
                  + f_6 * gf_4[k]
                  + pb_y[k] * gg_8[k];

        t_18[k] = -f_3 * gf_s_5[k]
                  + f_2 * gh_s_18[k]
                  + f_4 * gf_5[k]
                  + pb_y[k] * gg_9[k];

        t_19[k] = f_2 * gh_s_19[k]
                  + pb_y[k] * gg_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, fg_0, fh_0, gf_s_5, gh_s_20, \
                         gh_s_21, gh_s_22, gf_5, gg_10, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * gf_s_5[k]
                  + f_2 * gh_s_20[k]
                  + f_0 * gf_5[k]
                  + pb_z[k] * gg_10[k];

        t_21[k] = pa_y[k] * fh_0[k]
                  + f_2 * gh_s_21[k];

        t_22[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_22[k]
                  + pb_y[k] * gg_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_z, fg_1, fh_1, fh_2, gh_s_23, \
                         gh_s_24, gh_s_25, gh_s_26, gg_11, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * gh_s_23[k]
                  + pb_z[k] * gg_11[k];

        t_24[k] = f_6 * fg_1[k]
                  + pa_y[k] * fh_1[k]
                  + f_2 * gh_s_24[k];

        t_25[k] = f_2 * gh_s_25[k]
                  + pb_z[k] * gg_12[k];

        t_26[k] = pa_y[k] * fh_2[k]
                  + f_2 * gh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, fg_3, fg_4, fh_3, gh_s_27, \
                         gh_s_28, gh_s_29, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * fg_3[k]
                  + pa_y[k] * fh_3[k]
                  + f_2 * gh_s_27[k];

        t_28[k] = f_2 * gh_s_28[k]
                  + pb_z[k] * gg_13[k];

        t_29[k] = f_4 * fg_4[k]
                  + f_2 * gh_s_29[k]
                  + pb_y[k] * gg_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, fg_11, fh_5, gh_s_30, gh_s_31, \
                         gh_s_32, gg_15, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * fh_5[k]
                  + f_2 * gh_s_30[k];

        t_31[k] = f_7 * fg_11[k]
                  + f_2 * gh_s_31[k]
                  + pb_x[k] * gg_16[k];

        t_32[k] = f_2 * gh_s_32[k]
                  + pb_z[k] * gg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, fg_12, fg_13, fh_7, gh_s_33, gh_s_34, \
                         gh_s_35, gg_18, gg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * fg_12[k]
                  + f_2 * gh_s_33[k]
                  + pb_x[k] * gg_18[k];

        t_34[k] = f_7 * fg_13[k]
                  + f_2 * gh_s_34[k]
                  + pb_x[k] * gg_19[k];

        t_35[k] = pa_y[k] * fh_7[k]
                  + f_2 * gh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_z, dh_s_5, dh_5, fh_15, gf_s_7, gh_s_36, \
                         gh_s_37, gh_s_38, gf_7, gg_16, gg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * dh_s_5[k]
                  + f_6 * dh_5[k]
                  + pa_x[k] * fh_15[k]
                  + f_2 * gh_s_36[k];

        t_37[k] = f_2 * gh_s_37[k]
                  + pb_z[k] * gg_16[k];

        t_38[k] = -f_3 * gf_s_7[k]
                  + f_2 * gh_s_38[k]
                  + f_4 * gf_7[k]
                  + pb_z[k] * gg_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, pb_z, fg_7, fh_9, gf_s_8, gh_s_39, \
                         gh_s_40, gh_s_41, gf_8, gg_18, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_5 * gf_s_8[k]
                  + f_2 * gh_s_39[k]
                  + f_6 * gf_8[k]
                  + pb_z[k] * gg_18[k];

        t_40[k] = f_4 * fg_7[k]
                  + f_2 * gh_s_40[k]
                  + pb_y[k] * gg_20[k];

        t_41[k] = pa_y[k] * fh_9[k]
                  + f_2 * gh_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, fg_0, fh_0, fh_1, gh_s_42, \
                         gh_s_43, gh_s_44, gh_s_45, gg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * fh_0[k]
                  + f_2 * gh_s_42[k];

        t_43[k] = f_2 * gh_s_43[k]
                  + pb_y[k] * gg_21[k];

        t_44[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_44[k]
                  + pb_z[k] * gg_21[k];

        t_45[k] = pa_z[k] * fh_1[k]
                  + f_2 * gh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_y, fg_2, fg_3, fh_2, fh_3, fh_4, \
                         gh_s_46, gh_s_47, gh_s_48, gh_s_49, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * gh_s_46[k]
                  + pb_y[k] * gg_22[k];

        t_47[k] = f_6 * fg_2[k]
                  + pa_z[k] * fh_2[k]
                  + f_2 * gh_s_47[k];

        t_48[k] = pa_z[k] * fh_3[k]
                  + f_2 * gh_s_48[k];

        t_49[k] = f_4 * fg_3[k]
                  + pa_z[k] * fh_4[k]
                  + f_2 * gh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_z, pb_y, fg_4, fh_5, fh_6, gh_s_50, gh_s_51, \
                         gh_s_52, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * gh_s_50[k]
                  + pb_y[k] * gg_23[k];

        t_51[k] = f_7 * fg_4[k]
                  + pa_z[k] * fh_5[k]
                  + f_2 * gh_s_51[k];

        t_52[k] = pa_z[k] * fh_6[k]
                  + f_2 * gh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, fg_18, fg_19, gh_s_53, gh_s_54, \
                         gh_s_55, gg_24, gg_25, gg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * fg_18[k]
                  + f_2 * gh_s_53[k]
                  + pb_x[k] * gg_25[k];

        t_54[k] = f_7 * fg_19[k]
                  + f_2 * gh_s_54[k]
                  + pb_x[k] * gg_26[k];

        t_55[k] = f_2 * gh_s_55[k]
                  + pb_y[k] * gg_24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_z, pb_x, pb_y, fg_20, fh_8, gf_s_12, gh_s_56, \
                         gh_s_57, gh_s_58, gf_11, gg_25, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_7 * fg_20[k]
                  + f_2 * gh_s_56[k]
                  + pb_x[k] * gg_28[k];

        t_57[k] = pa_z[k] * fh_8[k]
                  + f_2 * gh_s_57[k];

        t_58[k] = -f_9 * gf_s_12[k]
                  + f_2 * gh_s_58[k]
                  + f_7 * gf_11[k]
                  + pb_y[k] * gg_25[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, gf_s_13, gf_s_14, gh_s_59, gh_s_60, gh_s_61, \
                         gf_12, gf_13, gg_26, gg_27, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_5 * gf_s_13[k]
                  + f_2 * gh_s_59[k]
                  + f_6 * gf_12[k]
                  + pb_y[k] * gg_26[k];

        t_60[k] = -f_3 * gf_s_14[k]
                  + f_2 * gh_s_60[k]
                  + f_4 * gf_13[k]
                  + pb_y[k] * gg_27[k];

        t_61[k] = f_2 * gh_s_61[k]
                  + pb_y[k] * gg_28[k];
    }

#pragma omp simd aligned(t_62, t_63, pa_x, pa_y, dh_s_0, dh_s_6, dh_0, dh_6, fh_10, fh_21, \
                         gh_s_62, gh_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -f_8 * dh_s_6[k]
                  + f_6 * dh_6[k]
                  + pa_x[k] * fh_21[k]
                  + f_2 * gh_s_62[k];

        t_63[k] = -f_10 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_10[k]
                  + f_2 * gh_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, pb_z, fg_8, fg_22, gf_s_17, gh_s_64, \
                         gh_s_65, gh_s_66, gf_16, gg_29, gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_6 * fg_8[k]
                  + f_2 * gh_s_64[k]
                  + pb_y[k] * gg_29[k];

        t_65[k] = f_2 * gh_s_65[k]
                  + pb_z[k] * gg_29[k];

        t_66[k] = f_6 * fg_22[k]
                  - f_5 * gf_s_17[k]
                  + f_2 * gh_s_66[k]
                  + f_6 * gf_16[k]
                  + pb_x[k] * gg_32[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pb_z, fg_24, gf_s_15, gf_s_18, gh_s_67, \
                         gh_s_68, gh_s_69, gf_14, gf_17, gg_30, gg_31, \
                         gg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * gh_s_67[k]
                  + pb_z[k] * gg_30[k];

        t_68[k] = -f_3 * gf_s_15[k]
                  + f_2 * gh_s_68[k]
                  + f_4 * gf_14[k]
                  + pb_z[k] * gg_31[k];

        t_69[k] = f_6 * fg_24[k]
                  - f_3 * gf_s_18[k]
                  + f_2 * gh_s_69[k]
                  + f_4 * gf_17[k]
                  + pb_x[k] * gg_34[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_y, pb_z, fg_10, gf_s_16, gh_s_70, gh_s_71, \
                         gh_s_72, gf_15, gg_32, gg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * gh_s_70[k]
                  + pb_z[k] * gg_32[k];

        t_71[k] = f_6 * fg_10[k]
                  + f_2 * gh_s_71[k]
                  + pb_y[k] * gg_33[k];

        t_72[k] = -f_5 * gf_s_16[k]
                  + f_2 * gh_s_72[k]
                  + f_6 * gf_15[k]
                  + pb_z[k] * gg_33[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, pb_z, fg_25, fg_26, gh_s_73, gh_s_74, \
                         gh_s_75, gg_34, gg_35, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * fg_25[k]
                  + f_2 * gh_s_73[k]
                  + pb_x[k] * gg_35[k];

        t_74[k] = f_2 * gh_s_74[k]
                  + pb_z[k] * gg_34[k];

        t_75[k] = f_6 * fg_26[k]
                  + f_2 * gh_s_75[k]
                  + pb_x[k] * gg_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_x, pb_x, dh_s_10, dh_10, fg_27, fg_28, fh_27, \
                         gh_s_76, gh_s_77, gh_s_78, gg_38, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_6 * fg_27[k]
                  + f_2 * gh_s_76[k]
                  + pb_x[k] * gg_38[k];

        t_77[k] = f_6 * fg_28[k]
                  + f_2 * gh_s_77[k]
                  + pb_x[k] * gg_39[k];

        t_78[k] = -f_10 * dh_s_10[k]
                  + f_4 * dh_10[k]
                  + pa_x[k] * fh_27[k]
                  + f_2 * gh_s_78[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_z, gf_s_18, gf_s_19, gh_s_79, gh_s_80, gh_s_81, \
                         gf_17, gf_18, gg_35, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_2 * gh_s_79[k]
                  + pb_z[k] * gg_35[k];

        t_80[k] = -f_3 * gf_s_18[k]
                  + f_2 * gh_s_80[k]
                  + f_4 * gf_17[k]
                  + pb_z[k] * gg_36[k];

        t_81[k] = -f_5 * gf_s_19[k]
                  + f_2 * gh_s_81[k]
                  + f_6 * gf_18[k]
                  + pb_z[k] * gg_37[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pb_y, pb_z, fg_14, fh_16, gf_s_20, gh_s_82, \
                         gh_s_83, gh_s_84, gf_19, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * fg_14[k]
                  + f_2 * gh_s_82[k]
                  + pb_y[k] * gg_39[k];

        t_83[k] = -f_1 * gf_s_20[k]
                  + f_2 * gh_s_83[k]
                  + f_0 * gf_19[k]
                  + pb_z[k] * gg_39[k];

        t_84[k] = pa_y[k] * fh_16[k]
                  + f_2 * gh_s_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pa_z, pb_y, fg_16, fh_11, fh_12, fh_17, \
                         gh_s_85, gh_s_86, gh_s_87, gh_s_88, gg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_z[k] * fh_11[k]
                  + f_2 * gh_s_85[k];

        t_86[k] = pa_y[k] * fh_17[k]
                  + f_2 * gh_s_86[k];

        t_87[k] = pa_z[k] * fh_12[k]
                  + f_2 * gh_s_87[k];

        t_88[k] = f_4 * fg_16[k]
                  + f_2 * gh_s_88[k]
                  + pb_y[k] * gg_40[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_y, pa_z, pb_z, fg_9, fh_13, fh_18, gh_s_89, \
                         gh_s_90, gh_s_91, gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * fh_18[k]
                  + f_2 * gh_s_89[k];

        t_90[k] = pa_z[k] * fh_13[k]
                  + f_2 * gh_s_90[k];

        t_91[k] = f_4 * fg_9[k]
                  + f_2 * gh_s_91[k]
                  + pb_z[k] * gg_41[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_y, pa_z, pb_y, fg_17, fh_14, fh_19, gh_s_92, \
                         gh_s_93, gh_s_94, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * fg_17[k]
                  + f_2 * gh_s_92[k]
                  + pb_y[k] * gg_42[k];

        t_93[k] = pa_y[k] * fh_19[k]
                  + f_2 * gh_s_93[k];

        t_94[k] = pa_z[k] * fh_14[k]
                  + f_2 * gh_s_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_x, fg_32, fg_33, fg_34, gh_s_95, gh_s_96, \
                         gh_s_97, gg_44, gg_45, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_6 * fg_32[k]
                  + f_2 * gh_s_95[k]
                  + pb_x[k] * gg_44[k];

        t_96[k] = f_6 * fg_33[k]
                  + f_2 * gh_s_96[k]
                  + pb_x[k] * gg_45[k];

        t_97[k] = f_6 * fg_34[k]
                  + f_2 * gh_s_97[k]
                  + pb_x[k] * gg_46[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_y, pa_z, pb_z, fg_11, fh_15, fh_20, gh_s_98, \
                         gh_s_99, gh_s_100, gg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pa_y[k] * fh_20[k]
                  + f_2 * gh_s_98[k];

        t_99[k] = pa_z[k] * fh_15[k]
                  + f_2 * gh_s_99[k];

        t_100[k] = f_4 * fg_11[k]
                   + f_2 * gh_s_100[k]
                   + pb_z[k] * gg_43[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_x, pb_y, dh_s_13, dh_s_14, dh_13, dh_14, \
                         fg_20, fh_28, fh_29, gh_s_101, gh_s_102, gh_s_103, \
                         gg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -f_10 * dh_s_13[k]
                   + f_4 * dh_13[k]
                   + pa_x[k] * fh_28[k]
                   + f_2 * gh_s_101[k];

        t_102[k] = -f_10 * dh_s_14[k]
                   + f_4 * dh_14[k]
                   + pa_x[k] * fh_29[k]
                   + f_2 * gh_s_102[k];

        t_103[k] = f_4 * fg_20[k]
                   + f_2 * gh_s_103[k]
                   + pb_y[k] * gg_47[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pa_z, pb_y, dh_s_0, dh_0, fh_16, fh_21, \
                         gh_s_104, gh_s_105, gh_s_106, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * fh_21[k]
                   + f_2 * gh_s_104[k];

        t_105[k] = -f_10 * dh_s_0[k]
                   + f_4 * dh_0[k]
                   + pa_z[k] * fh_16[k]
                   + f_2 * gh_s_105[k];

        t_106[k] = f_2 * gh_s_106[k]
                   + pb_y[k] * gg_48[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_y, pb_z, fg_15, gf_s_23, gh_s_107, gh_s_108, \
                         gh_s_109, gf_22, gg_48, gg_49, gg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * fg_15[k]
                   + f_2 * gh_s_107[k]
                   + pb_z[k] * gg_48[k];

        t_108[k] = -f_3 * gf_s_23[k]
                   + f_2 * gh_s_108[k]
                   + f_4 * gf_22[k]
                   + pb_y[k] * gg_49[k];

        t_109[k] = f_2 * gh_s_109[k]
                   + pb_y[k] * gg_50[k];
    }

#pragma omp simd aligned(t_110, t_111, pb_x, pb_y, fg_37, gf_s_24, gf_s_26, gh_s_110, \
                         gh_s_111, gf_23, gf_25, gg_51, gg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_6 * fg_37[k]
                   - f_5 * gf_s_26[k]
                   + f_2 * gh_s_110[k]
                   + f_6 * gf_25[k]
                   + pb_x[k] * gg_53[k];

        t_111[k] = -f_5 * gf_s_24[k]
                   + f_2 * gh_s_111[k]
                   + f_6 * gf_23[k]
                   + pb_y[k] * gg_51[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, pb_y, fg_38, gf_s_25, gf_s_30, gh_s_112, \
                         gh_s_113, gh_s_114, gf_24, gf_29, gg_52, gg_53, \
                         gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -f_3 * gf_s_25[k]
                   + f_2 * gh_s_112[k]
                   + f_4 * gf_24[k]
                   + pb_y[k] * gg_52[k];

        t_113[k] = f_2 * gh_s_113[k]
                   + pb_y[k] * gg_53[k];

        t_114[k] = f_6 * fg_38[k]
                   - f_3 * gf_s_30[k]
                   + f_2 * gh_s_114[k]
                   + f_4 * gf_29[k]
                   + pb_x[k] * gg_54[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_x, fg_39, fg_40, fg_41, gh_s_115, gh_s_116, \
                         gh_s_117, gg_55, gg_56, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_6 * fg_39[k]
                   + f_2 * gh_s_115[k]
                   + pb_x[k] * gg_55[k];

        t_116[k] = f_6 * fg_40[k]
                   + f_2 * gh_s_116[k]
                   + pb_x[k] * gg_56[k];

        t_117[k] = f_6 * fg_41[k]
                   + f_2 * gh_s_117[k]
                   + pb_x[k] * gg_57[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_y, fg_42, gf_s_27, gh_s_118, gh_s_119, \
                         gh_s_120, gf_26, gg_54, gg_55, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_2 * gh_s_118[k]
                   + pb_y[k] * gg_54[k];

        t_119[k] = f_6 * fg_42[k]
                   + f_2 * gh_s_119[k]
                   + pb_x[k] * gg_59[k];

        t_120[k] = -f_1 * gf_s_27[k]
                   + f_2 * gh_s_120[k]
                   + f_0 * gf_26[k]
                   + pb_y[k] * gg_55[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_y, gf_s_28, gf_s_29, gf_s_30, gh_s_121, \
                         gh_s_122, gh_s_123, gf_27, gf_28, gf_29, gg_56, gg_57, \
                         gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_9 * gf_s_28[k]
                   + f_2 * gh_s_121[k]
                   + f_7 * gf_27[k]
                   + pb_y[k] * gg_56[k];

        t_122[k] = -f_5 * gf_s_29[k]
                   + f_2 * gh_s_122[k]
                   + f_6 * gf_28[k]
                   + pb_y[k] * gg_57[k];

        t_123[k] = -f_3 * gf_s_30[k]
                   + f_2 * gh_s_123[k]
                   + f_4 * gf_29[k]
                   + pb_y[k] * gg_58[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pb_y, dh_s_24, dh_24, fg_43, fh_35, fh_36, \
                         gh_s_124, gh_s_125, gh_s_126, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * gh_s_124[k]
                   + pb_y[k] * gg_59[k];

        t_125[k] = -f_10 * dh_s_24[k]
                   + f_4 * dh_24[k]
                   + pa_x[k] * fh_35[k]
                   + f_2 * gh_s_125[k];

        t_126[k] = f_11 * fg_43[k]
                   + pa_x[k] * fh_36[k]
                   + f_2 * gh_s_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pb_y, pb_z, fg_21, fg_45, fh_38, \
                         gh_s_127, gh_s_128, gh_s_129, gh_s_130, gg_60, \
                         gg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_7 * fg_21[k]
                   + f_2 * gh_s_127[k]
                   + pb_y[k] * gg_60[k];

        t_128[k] = f_2 * gh_s_128[k]
                   + pb_z[k] * gg_60[k];

        t_129[k] = f_7 * fg_45[k]
                   + pa_x[k] * fh_38[k]
                   + f_2 * gh_s_129[k];

        t_130[k] = f_2 * gh_s_130[k]
                   + pb_z[k] * gg_61[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pa_x, pb_z, fg_47, fg_48, fh_40, fh_41, \
                         gh_s_131, gh_s_132, gh_s_133, gg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_7 * fg_47[k]
                   + pa_x[k] * fh_40[k]
                   + f_2 * gh_s_131[k];

        t_132[k] = f_6 * fg_48[k]
                   + pa_x[k] * fh_41[k]
                   + f_2 * gh_s_132[k];

        t_133[k] = f_2 * gh_s_133[k]
                   + pb_z[k] * gg_62[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_x, pb_x, pb_y, fg_23, fg_50, fg_51, fh_44, \
                         gh_s_134, gh_s_135, gh_s_136, gg_63, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_7 * fg_23[k]
                   + f_2 * gh_s_134[k]
                   + pb_y[k] * gg_63[k];

        t_135[k] = f_6 * fg_50[k]
                   + pa_x[k] * fh_44[k]
                   + f_2 * gh_s_135[k];

        t_136[k] = f_4 * fg_51[k]
                   + f_2 * gh_s_136[k]
                   + pb_x[k] * gg_65[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, pb_z, fg_53, fg_54, gh_s_137, gh_s_138, \
                         gh_s_139, gg_64, gg_66, gg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * gh_s_137[k]
                   + pb_z[k] * gg_64[k];

        t_138[k] = f_4 * fg_53[k]
                   + f_2 * gh_s_138[k]
                   + pb_x[k] * gg_66[k];

        t_139[k] = f_4 * fg_54[k]
                   + f_2 * gh_s_139[k]
                   + pb_x[k] * gg_67[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pb_x, pb_z, fg_55, fh_45, fh_46, \
                         gh_s_140, gh_s_141, gh_s_142, gh_s_143, gg_65, \
                         gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * fg_55[k]
                   + f_2 * gh_s_140[k]
                   + pb_x[k] * gg_68[k];

        t_141[k] = pa_x[k] * fh_45[k]
                   + f_2 * gh_s_141[k];

        t_142[k] = f_2 * gh_s_142[k]
                   + pb_z[k] * gg_65[k];

        t_143[k] = pa_x[k] * fh_46[k]
                   + f_2 * gh_s_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_x, pa_z, fh_22, fh_47, fh_48, fh_49, \
                         gh_s_144, gh_s_145, gh_s_146, gh_s_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pa_x[k] * fh_47[k]
                   + f_2 * gh_s_144[k];

        t_145[k] = pa_x[k] * fh_48[k]
                   + f_2 * gh_s_145[k];

        t_146[k] = pa_x[k] * fh_49[k]
                   + f_2 * gh_s_146[k];

        t_147[k] = pa_z[k] * fh_22[k]
                   + f_2 * gh_s_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_z, pb_z, fg_21, fh_23, fh_24, gh_s_148, \
                         gh_s_149, gh_s_150, gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pa_z[k] * fh_23[k]
                   + f_2 * gh_s_148[k];

        t_149[k] = f_4 * fg_21[k]
                   + f_2 * gh_s_149[k]
                   + pb_z[k] * gg_69[k];

        t_150[k] = pa_z[k] * fh_24[k]
                   + f_2 * gh_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_x, pa_z, pb_y, fg_29, fg_56, fh_25, fh_50, \
                         gh_s_151, gh_s_152, gh_s_153, gg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_6 * fg_29[k]
                   + f_2 * gh_s_151[k]
                   + pb_y[k] * gg_70[k];

        t_152[k] = f_7 * fg_56[k]
                   + pa_x[k] * fh_50[k]
                   + f_2 * gh_s_152[k];

        t_153[k] = pa_z[k] * fh_25[k]
                   + f_2 * gh_s_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_x, pb_y, pb_z, fg_22, fg_31, fg_57, fh_51, \
                         gh_s_154, gh_s_155, gh_s_156, gg_71, gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_4 * fg_22[k]
                   + f_2 * gh_s_154[k]
                   + pb_z[k] * gg_71[k];

        t_155[k] = f_6 * fg_31[k]
                   + f_2 * gh_s_155[k]
                   + pb_y[k] * gg_72[k];

        t_156[k] = f_6 * fg_57[k]
                   + pa_x[k] * fh_51[k]
                   + f_2 * gh_s_156[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pa_z, pb_x, fg_59, fg_60, fh_26, gh_s_157, \
                         gh_s_158, gh_s_159, gg_73, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pa_z[k] * fh_26[k]
                   + f_2 * gh_s_157[k];

        t_158[k] = f_4 * fg_59[k]
                   + f_2 * gh_s_158[k]
                   + pb_x[k] * gg_73[k];

        t_159[k] = f_4 * fg_60[k]
                   + f_2 * gh_s_159[k]
                   + pb_x[k] * gg_74[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_x, fg_61, fg_62, fh_52, fh_53, \
                         gh_s_160, gh_s_161, gh_s_162, gh_s_163, gg_75, \
                         gg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_4 * fg_61[k]
                   + f_2 * gh_s_160[k]
                   + pb_x[k] * gg_75[k];

        t_161[k] = f_4 * fg_62[k]
                   + f_2 * gh_s_161[k]
                   + pb_x[k] * gg_76[k];

        t_162[k] = pa_x[k] * fh_52[k]
                   + f_2 * gh_s_162[k];

        t_163[k] = pa_x[k] * fh_53[k]
                   + f_2 * gh_s_163[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_x, fh_54, fh_55, fh_56, fh_57, \
                         gh_s_164, gh_s_165, gh_s_166, gh_s_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_x[k] * fh_54[k]
                   + f_2 * gh_s_164[k];

        t_165[k] = pa_x[k] * fh_55[k]
                   + f_2 * gh_s_165[k];

        t_166[k] = pa_x[k] * fh_56[k]
                   + f_2 * gh_s_166[k];

        t_167[k] = pa_x[k] * fh_57[k]
                   + f_2 * gh_s_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pb_y, fg_35, fh_30, fh_31, gh_s_168, \
                         gh_s_169, gh_s_170, gg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * fh_30[k]
                   + f_2 * gh_s_168[k];

        t_169[k] = f_4 * fg_35[k]
                   + f_2 * gh_s_169[k]
                   + pb_y[k] * gg_77[k];

        t_170[k] = pa_y[k] * fh_31[k]
                   + f_2 * gh_s_170[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_x, pa_y, pb_y, fg_36, fg_63, fh_32, fh_58, \
                         gh_s_171, gh_s_172, gh_s_173, gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_7 * fg_63[k]
                   + pa_x[k] * fh_58[k]
                   + f_2 * gh_s_171[k];

        t_172[k] = f_4 * fg_36[k]
                   + f_2 * gh_s_172[k]
                   + pb_y[k] * gg_78[k];

        t_173[k] = pa_y[k] * fh_32[k]
                   + f_2 * gh_s_173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pa_x, pb_y, pb_z, fg_30, fg_37, fg_64, fh_59, \
                         gh_s_174, gh_s_175, gh_s_176, gg_79, gg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * fg_64[k]
                   + pa_x[k] * fh_59[k]
                   + f_2 * gh_s_174[k];

        t_175[k] = f_6 * fg_30[k]
                   + f_2 * gh_s_175[k]
                   + pb_z[k] * gg_79[k];

        t_176[k] = f_4 * fg_37[k]
                   + f_2 * gh_s_176[k]
                   + pb_y[k] * gg_80[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pb_x, fg_65, fg_66, fh_33, gh_s_177, \
                         gh_s_178, gh_s_179, gg_81, gg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_y[k] * fh_33[k]
                   + f_2 * gh_s_177[k];

        t_178[k] = f_4 * fg_65[k]
                   + f_2 * gh_s_178[k]
                   + pb_x[k] * gg_81[k];

        t_179[k] = f_4 * fg_66[k]
                   + f_2 * gh_s_179[k]
                   + pb_x[k] * gg_82[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_y, pb_x, fg_67, fg_68, fh_34, gh_s_180, \
                         gh_s_181, gh_s_182, gg_83, gg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * fg_67[k]
                   + f_2 * gh_s_180[k]
                   + pb_x[k] * gg_83[k];

        t_181[k] = f_4 * fg_68[k]
                   + f_2 * gh_s_181[k]
                   + pb_x[k] * gg_84[k];

        t_182[k] = pa_y[k] * fh_34[k]
                   + f_2 * gh_s_182[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pa_x, fh_60, fh_61, fh_62, fh_63, \
                         fh_64, gh_s_183, gh_s_184, gh_s_185, gh_s_186, \
                         gh_s_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_x[k] * fh_60[k]
                   + f_2 * gh_s_183[k];

        t_184[k] = pa_x[k] * fh_61[k]
                   + f_2 * gh_s_184[k];

        t_185[k] = pa_x[k] * fh_62[k]
                   + f_2 * gh_s_185[k];

        t_186[k] = pa_x[k] * fh_63[k]
                   + f_2 * gh_s_186[k];

        t_187[k] = pa_x[k] * fh_64[k]
                   + f_2 * gh_s_187[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_x, pb_y, pb_z, fg_35, fg_70, fh_65, \
                         fh_66, gh_s_188, gh_s_189, gh_s_190, gh_s_191, \
                         gg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_x[k] * fh_65[k]
                   + f_2 * gh_s_188[k];

        t_189[k] = f_11 * fg_70[k]
                   + pa_x[k] * fh_66[k]
                   + f_2 * gh_s_189[k];

        t_190[k] = f_2 * gh_s_190[k]
                   + pb_y[k] * gg_85[k];

        t_191[k] = f_7 * fg_35[k]
                   + f_2 * gh_s_191[k]
                   + pb_z[k] * gg_85[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_x, pb_y, fg_73, fg_75, fh_69, fh_71, \
                         gh_s_192, gh_s_193, gh_s_194, gg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_7 * fg_73[k]
                   + pa_x[k] * fh_69[k]
                   + f_2 * gh_s_192[k];

        t_193[k] = f_2 * gh_s_193[k]
                   + pb_y[k] * gg_86[k];

        t_194[k] = f_7 * fg_75[k]
                   + pa_x[k] * fh_71[k]
                   + f_2 * gh_s_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_x, pb_y, fg_76, fg_77, fh_72, fh_73, \
                         gh_s_195, gh_s_196, gh_s_197, gg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_6 * fg_76[k]
                   + pa_x[k] * fh_72[k]
                   + f_2 * gh_s_195[k];

        t_196[k] = f_6 * fg_77[k]
                   + pa_x[k] * fh_73[k]
                   + f_2 * gh_s_196[k];

        t_197[k] = f_2 * gh_s_197[k]
                   + pb_y[k] * gg_87[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_x, pb_x, fg_78, fg_79, fg_80, fh_75, \
                         gh_s_198, gh_s_199, gh_s_200, gg_89, gg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_6 * fg_78[k]
                   + pa_x[k] * fh_75[k]
                   + f_2 * gh_s_198[k];

        t_199[k] = f_4 * fg_79[k]
                   + f_2 * gh_s_199[k]
                   + pb_x[k] * gg_89[k];

        t_200[k] = f_4 * fg_80[k]
                   + f_2 * gh_s_200[k]
                   + pb_x[k] * gg_90[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_x, pb_y, fg_81, fg_83, gh_s_201, gh_s_202, \
                         gh_s_203, gg_88, gg_91, gg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_4 * fg_81[k]
                   + f_2 * gh_s_201[k]
                   + pb_x[k] * gg_91[k];

        t_202[k] = f_2 * gh_s_202[k]
                   + pb_y[k] * gg_88[k];

        t_203[k] = f_4 * fg_83[k]
                   + f_2 * gh_s_203[k]
                   + pb_x[k] * gg_92[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, fh_76, fh_77, fh_78, fh_79, \
                         gh_s_204, gh_s_205, gh_s_206, gh_s_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_x[k] * fh_76[k]
                   + f_2 * gh_s_204[k];

        t_205[k] = pa_x[k] * fh_77[k]
                   + f_2 * gh_s_205[k];

        t_206[k] = pa_x[k] * fh_78[k]
                   + f_2 * gh_s_206[k];

        t_207[k] = pa_x[k] * fh_79[k]
                   + f_2 * gh_s_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_x, pb_x, pb_y, fh_80, gf_s_39, gh_s_208, \
                         gh_s_209, gh_s_210, gf_35, gg_92, gg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_2 * gh_s_208[k]
                   + pb_y[k] * gg_92[k];

        t_209[k] = pa_x[k] * fh_80[k]
                   + f_2 * gh_s_209[k];

        t_210[k] = -f_1 * gf_s_39[k]
                   + f_2 * gh_s_210[k]
                   + f_0 * gf_35[k]
                   + pb_x[k] * gg_93[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pb_x, pb_z, gf_s_40, gf_s_41, gh_s_211, \
                         gh_s_212, gh_s_213, gf_36, gf_37, gg_93, gg_94, \
                         gg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -f_9 * gf_s_40[k]
                   + f_2 * gh_s_211[k]
                   + f_7 * gf_36[k]
                   + pb_x[k] * gg_94[k];

        t_212[k] = f_2 * gh_s_212[k]
                   + pb_z[k] * gg_93[k];

        t_213[k] = -f_5 * gf_s_41[k]
                   + f_2 * gh_s_213[k]
                   + f_6 * gf_37[k]
                   + pb_x[k] * gg_95[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, pb_z, gf_s_42, gf_s_43, gh_s_214, \
                         gh_s_215, gh_s_216, gf_38, gf_39, gg_94, gg_96, \
                         gg_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_2 * gh_s_214[k]
                   + pb_z[k] * gg_94[k];

        t_215[k] = -f_5 * gf_s_42[k]
                   + f_2 * gh_s_215[k]
                   + f_6 * gf_38[k]
                   + pb_x[k] * gg_96[k];

        t_216[k] = -f_3 * gf_s_43[k]
                   + f_2 * gh_s_216[k]
                   + f_4 * gf_39[k]
                   + pb_x[k] * gg_97[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pb_x, pb_z, gf_s_45, gf_s_46, gh_s_217, \
                         gh_s_218, gh_s_219, gf_41, gf_42, gg_95, gg_98, \
                         gg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_2 * gh_s_217[k]
                   + pb_z[k] * gg_95[k];

        t_218[k] = -f_3 * gf_s_45[k]
                   + f_2 * gh_s_218[k]
                   + f_4 * gf_41[k]
                   + pb_x[k] * gg_98[k];

        t_219[k] = -f_3 * gf_s_46[k]
                   + f_2 * gh_s_219[k]
                   + f_4 * gf_42[k]
                   + pb_x[k] * gg_99[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pb_x, gh_s_220, gh_s_221, \
                         gh_s_222, gh_s_223, gh_s_224, gg_100, gg_101, gg_102, gg_103, \
                         gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_2 * gh_s_220[k]
                   + pb_x[k] * gg_100[k];

        t_221[k] = f_2 * gh_s_221[k]
                   + pb_x[k] * gg_101[k];

        t_222[k] = f_2 * gh_s_222[k]
                   + pb_x[k] * gg_102[k];

        t_223[k] = f_2 * gh_s_223[k]
                   + pb_x[k] * gg_103[k];

        t_224[k] = f_2 * gh_s_224[k]
                   + pb_x[k] * gg_104[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pb_y, pb_z, fg_51, gf_s_43, gh_s_225, gh_s_226, \
                         gh_s_227, gf_39, gg_100, gg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * fg_51[k]
                   - f_1 * gf_s_43[k]
                   + f_2 * gh_s_225[k]
                   + f_0 * gf_39[k]
                   + pb_y[k] * gg_100[k];

        t_226[k] = f_2 * gh_s_226[k]
                   + pb_z[k] * gg_100[k];

        t_227[k] = -f_3 * gf_s_43[k]
                   + f_2 * gh_s_227[k]
                   + f_4 * gf_39[k]
                   + pb_z[k] * gg_101[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pb_y, pb_z, fg_55, gf_s_44, gf_s_46, gh_s_228, \
                         gh_s_229, gh_s_230, gf_40, gf_42, gg_102, \
                         gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = -f_5 * gf_s_44[k]
                   + f_2 * gh_s_228[k]
                   + f_6 * gf_40[k]
                   + pb_z[k] * gg_102[k];

        t_229[k] = f_0 * fg_55[k]
                   + f_2 * gh_s_229[k]
                   + pb_y[k] * gg_104[k];

        t_230[k] = -f_1 * gf_s_46[k]
                   + f_2 * gh_s_230[k]
                   + f_0 * gf_42[k]
                   + pb_z[k] * gg_104[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pa_z, pb_x, fh_36, fh_37, fh_38, gf_s_47, \
                         gh_s_231, gh_s_232, gh_s_233, gh_s_234, gf_43, \
                         gg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pa_z[k] * fh_36[k]
                   + f_2 * gh_s_231[k];

        t_232[k] = pa_z[k] * fh_37[k]
                   + f_2 * gh_s_232[k];

        t_233[k] = -f_9 * gf_s_47[k]
                   + f_2 * gh_s_233[k]
                   + f_7 * gf_43[k]
                   + pb_x[k] * gg_105[k];

        t_234[k] = pa_z[k] * fh_38[k]
                   + f_2 * gh_s_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_z, pb_x, fg_44, fh_39, fh_41, gf_s_49, \
                         gh_s_235, gh_s_236, gh_s_237, gf_44, gg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_4 * fg_44[k]
                   + pa_z[k] * fh_39[k]
                   + f_2 * gh_s_235[k];

        t_236[k] = -f_5 * gf_s_49[k]
                   + f_2 * gh_s_236[k]
                   + f_6 * gf_44[k]
                   + pb_x[k] * gg_106[k];

        t_237[k] = pa_z[k] * fh_41[k]
                   + f_2 * gh_s_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_z, pb_x, fg_45, fg_46, fh_42, fh_43, gf_s_53, \
                         gh_s_238, gh_s_239, gh_s_240, gf_46, gg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_4 * fg_45[k]
                   + pa_z[k] * fh_42[k]
                   + f_2 * gh_s_238[k];

        t_239[k] = f_6 * fg_46[k]
                   + pa_z[k] * fh_43[k]
                   + f_2 * gh_s_239[k];

        t_240[k] = -f_3 * gf_s_53[k]
                   + f_2 * gh_s_240[k]
                   + f_4 * gf_46[k]
                   + pb_x[k] * gg_107[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pb_x, gh_s_241, gh_s_242, \
                         gh_s_243, gh_s_244, gh_s_245, gg_108, gg_109, gg_110, gg_111, \
                         gg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_2 * gh_s_241[k]
                   + pb_x[k] * gg_108[k];

        t_242[k] = f_2 * gh_s_242[k]
                   + pb_x[k] * gg_109[k];

        t_243[k] = f_2 * gh_s_243[k]
                   + pb_x[k] * gg_110[k];

        t_244[k] = f_2 * gh_s_244[k]
                   + pb_x[k] * gg_111[k];

        t_245[k] = f_2 * gh_s_245[k]
                   + pb_x[k] * gg_112[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_z, pb_z, fg_51, fg_52, fh_45, fh_46, \
                         gh_s_246, gh_s_247, gh_s_248, gg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_z[k] * fh_45[k]
                   + f_2 * gh_s_246[k];

        t_247[k] = f_4 * fg_51[k]
                   + f_2 * gh_s_247[k]
                   + pb_z[k] * gg_108[k];

        t_248[k] = f_6 * fg_52[k]
                   + pa_z[k] * fh_46[k]
                   + f_2 * gh_s_248[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_y, pa_z, pb_y, dh_s_15, dh_15, fg_53, fg_62, \
                         fh_47, fh_57, gh_s_249, gh_s_250, gh_s_251, \
                         gg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_7 * fg_53[k]
                   + pa_z[k] * fh_47[k]
                   + f_2 * gh_s_249[k];

        t_250[k] = f_7 * fg_62[k]
                   + f_2 * gh_s_250[k]
                   + pb_y[k] * gg_112[k];

        t_251[k] = -f_8 * dh_s_15[k]
                   + f_6 * dh_15[k]
                   + pa_y[k] * fh_57[k]
                   + f_2 * gh_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pb_x, gf_s_54, gf_s_55, gf_s_56, gh_s_252, \
                         gh_s_253, gh_s_254, gf_47, gf_48, gf_49, gg_113, gg_114, \
                         gg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -f_1 * gf_s_54[k]
                   + f_2 * gh_s_252[k]
                   + f_0 * gf_47[k]
                   + pb_x[k] * gg_113[k];

        t_253[k] = -f_9 * gf_s_55[k]
                   + f_2 * gh_s_253[k]
                   + f_7 * gf_48[k]
                   + pb_x[k] * gg_114[k];

        t_254[k] = -f_9 * gf_s_56[k]
                   + f_2 * gh_s_254[k]
                   + f_7 * gf_49[k]
                   + pb_x[k] * gg_115[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pb_x, gf_s_57, gf_s_58, gf_s_59, gh_s_255, \
                         gh_s_256, gh_s_257, gf_50, gf_51, gf_52, gg_116, gg_117, \
                         gg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -f_5 * gf_s_57[k]
                   + f_2 * gh_s_255[k]
                   + f_6 * gf_50[k]
                   + pb_x[k] * gg_116[k];

        t_256[k] = -f_5 * gf_s_58[k]
                   + f_2 * gh_s_256[k]
                   + f_6 * gf_51[k]
                   + pb_x[k] * gg_117[k];

        t_257[k] = -f_5 * gf_s_59[k]
                   + f_2 * gh_s_257[k]
                   + f_6 * gf_52[k]
                   + pb_x[k] * gg_118[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pb_x, gf_s_60, gf_s_61, gf_s_62, gh_s_258, \
                         gh_s_259, gh_s_260, gf_53, gf_54, gf_55, gg_119, gg_120, \
                         gg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -f_3 * gf_s_60[k]
                   + f_2 * gh_s_258[k]
                   + f_4 * gf_53[k]
                   + pb_x[k] * gg_119[k];

        t_259[k] = -f_3 * gf_s_61[k]
                   + f_2 * gh_s_259[k]
                   + f_4 * gf_54[k]
                   + pb_x[k] * gg_120[k];

        t_260[k] = -f_3 * gf_s_62[k]
                   + f_2 * gh_s_260[k]
                   + f_4 * gf_55[k]
                   + pb_x[k] * gg_121[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_x, gf_s_63, gh_s_261, gh_s_262, \
                         gh_s_263, gh_s_264, gf_56, gg_122, gg_123, gg_124, \
                         gg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -f_3 * gf_s_63[k]
                   + f_2 * gh_s_261[k]
                   + f_4 * gf_56[k]
                   + pb_x[k] * gg_122[k];

        t_262[k] = f_2 * gh_s_262[k]
                   + pb_x[k] * gg_123[k];

        t_263[k] = f_2 * gh_s_263[k]
                   + pb_x[k] * gg_124[k];

        t_264[k] = f_2 * gh_s_264[k]
                   + pb_x[k] * gg_125[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pa_z, pb_x, dh_s_10, dh_10, fh_52, gh_s_265, \
                         gh_s_266, gh_s_267, gg_126, gg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_2 * gh_s_265[k]
                   + pb_x[k] * gg_126[k];

        t_266[k] = f_2 * gh_s_266[k]
                   + pb_x[k] * gg_127[k];

        t_267[k] = -f_10 * dh_s_10[k]
                   + f_4 * dh_10[k]
                   + pa_z[k] * fh_52[k]
                   + f_2 * gh_s_267[k];
    }

#pragma omp simd aligned(t_268, t_269, pb_y, pb_z, fg_58, fg_67, gf_s_62, gh_s_268, gh_s_269, \
                         gf_55, gg_123, gg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_6 * fg_58[k]
                   + f_2 * gh_s_268[k]
                   + pb_z[k] * gg_123[k];

        t_269[k] = f_6 * fg_67[k]
                   - f_5 * gf_s_62[k]
                   + f_2 * gh_s_269[k]
                   + f_6 * gf_55[k]
                   + pb_y[k] * gg_125[k];
    }

#pragma omp simd aligned(t_270, t_271, pb_y, fg_68, fg_69, gf_s_63, gh_s_270, gh_s_271, gf_56, \
                         gg_126, gg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_6 * fg_68[k]
                   - f_3 * gf_s_63[k]
                   + f_2 * gh_s_270[k]
                   + f_4 * gf_56[k]
                   + pb_y[k] * gg_126[k];

        t_271[k] = f_6 * fg_69[k]
                   + f_2 * gh_s_271[k]
                   + pb_y[k] * gg_127[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_y, dh_s_24, dh_24, fg_70, fh_65, \
                         fh_66, fh_67, fh_68, gh_s_272, gh_s_273, gh_s_274, \
                         gh_s_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -f_10 * dh_s_24[k]
                   + f_4 * dh_24[k]
                   + pa_y[k] * fh_65[k]
                   + f_2 * gh_s_272[k];

        t_273[k] = pa_y[k] * fh_66[k]
                   + f_2 * gh_s_273[k];

        t_274[k] = f_4 * fg_70[k]
                   + pa_y[k] * fh_67[k]
                   + f_2 * gh_s_274[k];

        t_275[k] = pa_y[k] * fh_68[k]
                   + f_2 * gh_s_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_y, fg_71, fg_72, fg_73, fh_69, fh_70, \
                         fh_71, fh_72, gh_s_276, gh_s_277, gh_s_278, \
                         gh_s_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * fg_71[k]
                   + pa_y[k] * fh_69[k]
                   + f_2 * gh_s_276[k];

        t_277[k] = f_4 * fg_72[k]
                   + pa_y[k] * fh_70[k]
                   + f_2 * gh_s_277[k];

        t_278[k] = pa_y[k] * fh_71[k]
                   + f_2 * gh_s_278[k];

        t_279[k] = f_7 * fg_73[k]
                   + pa_y[k] * fh_72[k]
                   + f_2 * gh_s_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, fg_74, fg_75, fh_73, fh_74, \
                         fh_75, gh_s_280, gh_s_281, gh_s_282, gh_s_283, \
                         gg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_6 * fg_74[k]
                   + pa_y[k] * fh_73[k]
                   + f_2 * gh_s_280[k];

        t_281[k] = f_4 * fg_75[k]
                   + pa_y[k] * fh_74[k]
                   + f_2 * gh_s_281[k];

        t_282[k] = pa_y[k] * fh_75[k]
                   + f_2 * gh_s_282[k];

        t_283[k] = f_2 * gh_s_283[k]
                   + pb_x[k] * gg_128[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pb_x, gh_s_284, gh_s_285, gh_s_286, \
                         gh_s_287, gg_129, gg_130, gg_131, gg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_2 * gh_s_284[k]
                   + pb_x[k] * gg_129[k];

        t_285[k] = f_2 * gh_s_285[k]
                   + pb_x[k] * gg_130[k];

        t_286[k] = f_2 * gh_s_286[k]
                   + pb_x[k] * gg_131[k];

        t_287[k] = f_2 * gh_s_287[k]
                   + pb_x[k] * gg_132[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_y, pb_z, fg_65, fg_79, fg_81, fh_76, fh_78, \
                         gh_s_288, gh_s_289, gh_s_290, gg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_11 * fg_79[k]
                   + pa_y[k] * fh_76[k]
                   + f_2 * gh_s_288[k];

        t_289[k] = f_7 * fg_65[k]
                   + f_2 * gh_s_289[k]
                   + pb_z[k] * gg_128[k];

        t_290[k] = f_7 * fg_81[k]
                   + pa_y[k] * fh_78[k]
                   + f_2 * gh_s_290[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pb_y, fg_82, fg_83, fh_79, fh_80, \
                         gh_s_291, gh_s_292, gh_s_293, gg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_6 * fg_82[k]
                   + pa_y[k] * fh_79[k]
                   + f_2 * gh_s_291[k];

        t_292[k] = f_4 * fg_83[k]
                   + f_2 * gh_s_292[k]
                   + pb_y[k] * gg_132[k];

        t_293[k] = pa_y[k] * fh_80[k]
                   + f_2 * gh_s_293[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pb_x, pb_y, gf_s_71, gf_s_72, gh_s_294, \
                         gh_s_295, gh_s_296, gf_59, gf_60, gg_133, \
                         gg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = -f_1 * gf_s_71[k]
                   + f_2 * gh_s_294[k]
                   + f_0 * gf_59[k]
                   + pb_x[k] * gg_133[k];

        t_295[k] = f_2 * gh_s_295[k]
                   + pb_y[k] * gg_133[k];

        t_296[k] = -f_9 * gf_s_72[k]
                   + f_2 * gh_s_296[k]
                   + f_7 * gf_60[k]
                   + pb_x[k] * gg_134[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pb_y, gf_s_73, gf_s_74, gh_s_297, \
                         gh_s_298, gh_s_299, gf_61, gf_62, gg_134, gg_135, \
                         gg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -f_5 * gf_s_73[k]
                   + f_2 * gh_s_297[k]
                   + f_6 * gf_61[k]
                   + pb_x[k] * gg_135[k];

        t_298[k] = f_2 * gh_s_298[k]
                   + pb_y[k] * gg_134[k];

        t_299[k] = -f_5 * gf_s_74[k]
                   + f_2 * gh_s_299[k]
                   + f_6 * gf_62[k]
                   + pb_x[k] * gg_136[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pb_x, pb_y, gf_s_75, gf_s_76, gh_s_300, \
                         gh_s_301, gh_s_302, gf_63, gf_64, gg_136, gg_137, \
                         gg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -f_3 * gf_s_75[k]
                   + f_2 * gh_s_300[k]
                   + f_4 * gf_63[k]
                   + pb_x[k] * gg_137[k];

        t_301[k] = -f_3 * gf_s_76[k]
                   + f_2 * gh_s_301[k]
                   + f_4 * gf_64[k]
                   + pb_x[k] * gg_138[k];

        t_302[k] = f_2 * gh_s_302[k]
                   + pb_y[k] * gg_136[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pb_x, gf_s_78, gh_s_303, gh_s_304, \
                         gh_s_305, gh_s_306, gf_66, gg_139, gg_140, gg_141, \
                         gg_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -f_3 * gf_s_78[k]
                   + f_2 * gh_s_303[k]
                   + f_4 * gf_66[k]
                   + pb_x[k] * gg_139[k];

        t_304[k] = f_2 * gh_s_304[k]
                   + pb_x[k] * gg_140[k];

        t_305[k] = f_2 * gh_s_305[k]
                   + pb_x[k] * gg_141[k];

        t_306[k] = f_2 * gh_s_306[k]
                   + pb_x[k] * gg_142[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, pb_x, pb_y, gf_s_75, gh_s_307, gh_s_308, \
                         gh_s_309, gf_63, gg_140, gg_143, gg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_2 * gh_s_307[k]
                   + pb_x[k] * gg_143[k];

        t_308[k] = f_2 * gh_s_308[k]
                   + pb_x[k] * gg_144[k];

        t_309[k] = -f_1 * gf_s_75[k]
                   + f_2 * gh_s_309[k]
                   + f_0 * gf_63[k]
                   + pb_y[k] * gg_140[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pb_y, gf_s_76, gf_s_77, gf_s_78, gh_s_310, \
                         gh_s_311, gh_s_312, gf_64, gf_65, gf_66, gg_141, gg_142, \
                         gg_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -f_9 * gf_s_76[k]
                   + f_2 * gh_s_310[k]
                   + f_7 * gf_64[k]
                   + pb_y[k] * gg_141[k];

        t_311[k] = -f_5 * gf_s_77[k]
                   + f_2 * gh_s_311[k]
                   + f_6 * gf_65[k]
                   + pb_y[k] * gg_142[k];

        t_312[k] = -f_3 * gf_s_78[k]
                   + f_2 * gh_s_312[k]
                   + f_4 * gf_66[k]
                   + pb_y[k] * gg_143[k];
    }

#pragma omp simd aligned(t_313, t_314, pb_y, pb_z, fg_83, gf_s_78, gh_s_313, gh_s_314, gf_66, \
                         gg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_2 * gh_s_313[k]
                   + pb_y[k] * gg_144[k];

        t_314[k] = f_0 * fg_83[k]
                   - f_1 * gf_s_78[k]
                   + f_2 * gh_s_314[k]
                   + f_0 * gf_66[k]
                   + pb_z[k] * gg_144[k];
    }
}

auto
compute_prim_gh_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_28 = buffer.data(dh_s + 28);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_37 = buffer.data(dh_s + 37);
    const auto *dh_s_39 = buffer.data(dh_s + 39);
    const auto *dh_s_55 = buffer.data(dh_s + 55);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_55 = buffer.data(dh + 55);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
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
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_120 = buffer.data(fh + 120);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_3 = buffer.data(gf_s + 3);
    const auto *gf_s_4 = buffer.data(gf_s + 4);
    const auto *gf_s_5 = buffer.data(gf_s + 5);
    const auto *gf_s_6 = buffer.data(gf_s + 6);
    const auto *gf_s_7 = buffer.data(gf_s + 7);
    const auto *gf_s_10 = buffer.data(gf_s + 10);
    const auto *gf_s_11 = buffer.data(gf_s + 11);
    const auto *gf_s_12 = buffer.data(gf_s + 12);
    const auto *gf_s_13 = buffer.data(gf_s + 13);
    const auto *gf_s_14 = buffer.data(gf_s + 14);
    const auto *gf_s_15 = buffer.data(gf_s + 15);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_17 = buffer.data(gf_s + 17);
    const auto *gf_s_18 = buffer.data(gf_s + 18);
    const auto *gf_s_19 = buffer.data(gf_s + 19);
    const auto *gf_s_20 = buffer.data(gf_s + 20);
    const auto *gf_s_21 = buffer.data(gf_s + 21);
    const auto *gf_s_22 = buffer.data(gf_s + 22);
    const auto *gf_s_23 = buffer.data(gf_s + 23);
    const auto *gf_s_24 = buffer.data(gf_s + 24);
    const auto *gf_s_25 = buffer.data(gf_s + 25);
    const auto *gf_s_26 = buffer.data(gf_s + 26);
    const auto *gf_s_32 = buffer.data(gf_s + 32);
    const auto *gf_s_33 = buffer.data(gf_s + 33);
    const auto *gf_s_34 = buffer.data(gf_s + 34);
    const auto *gf_s_35 = buffer.data(gf_s + 35);
    const auto *gf_s_36 = buffer.data(gf_s + 36);
    const auto *gf_s_37 = buffer.data(gf_s + 37);
    const auto *gf_s_38 = buffer.data(gf_s + 38);
    const auto *gf_s_39 = buffer.data(gf_s + 39);
    const auto *gf_s_40 = buffer.data(gf_s + 40);
    const auto *gf_s_42 = buffer.data(gf_s + 42);
    const auto *gf_s_46 = buffer.data(gf_s + 46);
    const auto *gf_s_47 = buffer.data(gf_s + 47);
    const auto *gf_s_48 = buffer.data(gf_s + 48);
    const auto *gf_s_49 = buffer.data(gf_s + 49);
    const auto *gf_s_50 = buffer.data(gf_s + 50);
    const auto *gf_s_51 = buffer.data(gf_s + 51);
    const auto *gf_s_52 = buffer.data(gf_s + 52);
    const auto *gf_s_53 = buffer.data(gf_s + 53);
    const auto *gf_s_54 = buffer.data(gf_s + 54);
    const auto *gf_s_55 = buffer.data(gf_s + 55);
    const auto *gf_s_56 = buffer.data(gf_s + 56);
    const auto *gf_s_64 = buffer.data(gf_s + 64);
    const auto *gf_s_65 = buffer.data(gf_s + 65);
    const auto *gf_s_66 = buffer.data(gf_s + 66);
    const auto *gf_s_67 = buffer.data(gf_s + 67);
    const auto *gf_s_68 = buffer.data(gf_s + 68);
    const auto *gf_s_69 = buffer.data(gf_s + 69);
    const auto *gf_s_70 = buffer.data(gf_s + 70);
    const auto *gf_s_71 = buffer.data(gf_s + 71);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);
    const auto *gh_s_58 = buffer.data(gh_s + 58);
    const auto *gh_s_59 = buffer.data(gh_s + 59);
    const auto *gh_s_60 = buffer.data(gh_s + 60);
    const auto *gh_s_61 = buffer.data(gh_s + 61);
    const auto *gh_s_62 = buffer.data(gh_s + 62);
    const auto *gh_s_63 = buffer.data(gh_s + 63);
    const auto *gh_s_64 = buffer.data(gh_s + 64);
    const auto *gh_s_65 = buffer.data(gh_s + 65);
    const auto *gh_s_66 = buffer.data(gh_s + 66);
    const auto *gh_s_67 = buffer.data(gh_s + 67);
    const auto *gh_s_68 = buffer.data(gh_s + 68);
    const auto *gh_s_69 = buffer.data(gh_s + 69);
    const auto *gh_s_70 = buffer.data(gh_s + 70);
    const auto *gh_s_71 = buffer.data(gh_s + 71);
    const auto *gh_s_72 = buffer.data(gh_s + 72);
    const auto *gh_s_73 = buffer.data(gh_s + 73);
    const auto *gh_s_74 = buffer.data(gh_s + 74);
    const auto *gh_s_75 = buffer.data(gh_s + 75);
    const auto *gh_s_76 = buffer.data(gh_s + 76);
    const auto *gh_s_77 = buffer.data(gh_s + 77);
    const auto *gh_s_78 = buffer.data(gh_s + 78);
    const auto *gh_s_79 = buffer.data(gh_s + 79);
    const auto *gh_s_80 = buffer.data(gh_s + 80);
    const auto *gh_s_81 = buffer.data(gh_s + 81);
    const auto *gh_s_82 = buffer.data(gh_s + 82);
    const auto *gh_s_83 = buffer.data(gh_s + 83);
    const auto *gh_s_84 = buffer.data(gh_s + 84);
    const auto *gh_s_85 = buffer.data(gh_s + 85);
    const auto *gh_s_86 = buffer.data(gh_s + 86);
    const auto *gh_s_88 = buffer.data(gh_s + 88);
    const auto *gh_s_89 = buffer.data(gh_s + 89);
    const auto *gh_s_90 = buffer.data(gh_s + 90);
    const auto *gh_s_92 = buffer.data(gh_s + 92);
    const auto *gh_s_93 = buffer.data(gh_s + 93);
    const auto *gh_s_94 = buffer.data(gh_s + 94);
    const auto *gh_s_95 = buffer.data(gh_s + 95);
    const auto *gh_s_96 = buffer.data(gh_s + 96);
    const auto *gh_s_97 = buffer.data(gh_s + 97);
    const auto *gh_s_98 = buffer.data(gh_s + 98);
    const auto *gh_s_99 = buffer.data(gh_s + 99);
    const auto *gh_s_100 = buffer.data(gh_s + 100);
    const auto *gh_s_101 = buffer.data(gh_s + 101);
    const auto *gh_s_102 = buffer.data(gh_s + 102);
    const auto *gh_s_103 = buffer.data(gh_s + 103);
    const auto *gh_s_104 = buffer.data(gh_s + 104);
    const auto *gh_s_105 = buffer.data(gh_s + 105);
    const auto *gh_s_106 = buffer.data(gh_s + 106);
    const auto *gh_s_107 = buffer.data(gh_s + 107);
    const auto *gh_s_108 = buffer.data(gh_s + 108);
    const auto *gh_s_109 = buffer.data(gh_s + 109);
    const auto *gh_s_110 = buffer.data(gh_s + 110);
    const auto *gh_s_111 = buffer.data(gh_s + 111);
    const auto *gh_s_112 = buffer.data(gh_s + 112);
    const auto *gh_s_113 = buffer.data(gh_s + 113);
    const auto *gh_s_114 = buffer.data(gh_s + 114);
    const auto *gh_s_115 = buffer.data(gh_s + 115);
    const auto *gh_s_116 = buffer.data(gh_s + 116);
    const auto *gh_s_117 = buffer.data(gh_s + 117);
    const auto *gh_s_118 = buffer.data(gh_s + 118);
    const auto *gh_s_119 = buffer.data(gh_s + 119);
    const auto *gh_s_120 = buffer.data(gh_s + 120);
    const auto *gh_s_121 = buffer.data(gh_s + 121);
    const auto *gh_s_123 = buffer.data(gh_s + 123);
    const auto *gh_s_126 = buffer.data(gh_s + 126);
    const auto *gh_s_130 = buffer.data(gh_s + 130);
    const auto *gh_s_131 = buffer.data(gh_s + 131);
    const auto *gh_s_132 = buffer.data(gh_s + 132);
    const auto *gh_s_133 = buffer.data(gh_s + 133);
    const auto *gh_s_134 = buffer.data(gh_s + 134);
    const auto *gh_s_135 = buffer.data(gh_s + 135);
    const auto *gh_s_136 = buffer.data(gh_s + 136);
    const auto *gh_s_137 = buffer.data(gh_s + 137);
    const auto *gh_s_138 = buffer.data(gh_s + 138);
    const auto *gh_s_139 = buffer.data(gh_s + 139);
    const auto *gh_s_140 = buffer.data(gh_s + 140);
    const auto *gh_s_141 = buffer.data(gh_s + 141);
    const auto *gh_s_142 = buffer.data(gh_s + 142);
    const auto *gh_s_143 = buffer.data(gh_s + 143);
    const auto *gh_s_144 = buffer.data(gh_s + 144);
    const auto *gh_s_145 = buffer.data(gh_s + 145);
    const auto *gh_s_146 = buffer.data(gh_s + 146);
    const auto *gh_s_147 = buffer.data(gh_s + 147);
    const auto *gh_s_148 = buffer.data(gh_s + 148);
    const auto *gh_s_149 = buffer.data(gh_s + 149);
    const auto *gh_s_150 = buffer.data(gh_s + 150);
    const auto *gh_s_151 = buffer.data(gh_s + 151);
    const auto *gh_s_152 = buffer.data(gh_s + 152);
    const auto *gh_s_153 = buffer.data(gh_s + 153);
    const auto *gh_s_154 = buffer.data(gh_s + 154);
    const auto *gh_s_156 = buffer.data(gh_s + 156);
    const auto *gh_s_159 = buffer.data(gh_s + 159);
    const auto *gh_s_163 = buffer.data(gh_s + 163);
    const auto *gh_s_164 = buffer.data(gh_s + 164);
    const auto *gh_s_165 = buffer.data(gh_s + 165);
    const auto *gh_s_166 = buffer.data(gh_s + 166);
    const auto *gh_s_167 = buffer.data(gh_s + 167);
    const auto *gh_s_168 = buffer.data(gh_s + 168);
    const auto *gh_s_169 = buffer.data(gh_s + 169);
    const auto *gh_s_170 = buffer.data(gh_s + 170);
    const auto *gh_s_171 = buffer.data(gh_s + 171);
    const auto *gh_s_172 = buffer.data(gh_s + 172);
    const auto *gh_s_173 = buffer.data(gh_s + 173);
    const auto *gh_s_174 = buffer.data(gh_s + 174);
    const auto *gh_s_175 = buffer.data(gh_s + 175);
    const auto *gh_s_176 = buffer.data(gh_s + 176);
    const auto *gh_s_177 = buffer.data(gh_s + 177);
    const auto *gh_s_178 = buffer.data(gh_s + 178);
    const auto *gh_s_179 = buffer.data(gh_s + 179);
    const auto *gh_s_180 = buffer.data(gh_s + 180);
    const auto *gh_s_181 = buffer.data(gh_s + 181);
    const auto *gh_s_182 = buffer.data(gh_s + 182);
    const auto *gh_s_183 = buffer.data(gh_s + 183);
    const auto *gh_s_184 = buffer.data(gh_s + 184);
    const auto *gh_s_185 = buffer.data(gh_s + 185);
    const auto *gh_s_186 = buffer.data(gh_s + 186);
    const auto *gh_s_187 = buffer.data(gh_s + 187);
    const auto *gh_s_188 = buffer.data(gh_s + 188);
    const auto *gh_s_189 = buffer.data(gh_s + 189);
    const auto *gh_s_190 = buffer.data(gh_s + 190);
    const auto *gh_s_201 = buffer.data(gh_s + 201);
    const auto *gh_s_202 = buffer.data(gh_s + 202);
    const auto *gh_s_203 = buffer.data(gh_s + 203);
    const auto *gh_s_204 = buffer.data(gh_s + 204);
    const auto *gh_s_205 = buffer.data(gh_s + 205);
    const auto *gh_s_206 = buffer.data(gh_s + 206);
    const auto *gh_s_207 = buffer.data(gh_s + 207);
    const auto *gh_s_208 = buffer.data(gh_s + 208);
    const auto *gh_s_209 = buffer.data(gh_s + 209);
    const auto *gh_s_210 = buffer.data(gh_s + 210);
    const auto *gh_s_211 = buffer.data(gh_s + 211);
    const auto *gh_s_212 = buffer.data(gh_s + 212);
    const auto *gh_s_213 = buffer.data(gh_s + 213);
    const auto *gh_s_214 = buffer.data(gh_s + 214);
    const auto *gh_s_215 = buffer.data(gh_s + 215);
    const auto *gh_s_216 = buffer.data(gh_s + 216);
    const auto *gh_s_217 = buffer.data(gh_s + 217);
    const auto *gh_s_218 = buffer.data(gh_s + 218);
    const auto *gh_s_219 = buffer.data(gh_s + 219);
    const auto *gh_s_220 = buffer.data(gh_s + 220);
    const auto *gh_s_221 = buffer.data(gh_s + 221);
    const auto *gh_s_222 = buffer.data(gh_s + 222);
    const auto *gh_s_223 = buffer.data(gh_s + 223);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gh_s_3, gf_0, gg_0, gg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_2 * gh_s_1[k]
                 + pb_y[k] * gg_0[k];

        t_2[k] = f_2 * gh_s_2[k]
                 + pb_z[k] * gg_0[k];

        t_3[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_3[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, gf_s_0, gf_s_1, gh_s_4, gh_s_5, gh_s_6, \
                         gf_0, gf_1, gg_2, gg_3, gg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_4[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];

        t_5[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_5[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = f_2 * gh_s_6[k]
                 + pb_y[k] * gg_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_z, fg_5, fg_8, gf_s_2, gh_s_7, gh_s_8, \
                         gh_s_9, gf_2, gg_4, gg_5, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_7[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_4[k];

        t_8[k] = f_0 * fg_5[k]
                 + f_2 * gh_s_8[k]
                 + pb_x[k] * gg_5[k];

        t_9[k] = f_0 * fg_8[k]
                 + f_2 * gh_s_9[k]
                 + pb_x[k] * gg_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, gf_s_3, gf_s_4, gf_s_5, gh_s_10, gh_s_11, \
                         gh_s_12, gf_3, gf_4, gf_5, gg_5, gg_6, gg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * gf_s_3[k]
                  + f_2 * gh_s_10[k]
                  + f_0 * gf_3[k]
                  + pb_y[k] * gg_5[k];

        t_11[k] = -f_5 * gf_s_4[k]
                  + f_2 * gh_s_11[k]
                  + f_6 * gf_4[k]
                  + pb_y[k] * gg_6[k];

        t_12[k] = -f_3 * gf_s_5[k]
                  + f_2 * gh_s_12[k]
                  + f_4 * gf_5[k]
                  + pb_y[k] * gg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, fh_0, gf_s_5, gh_s_13, gh_s_14, \
                         gh_s_15, gf_5, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * gh_s_13[k]
                  + pb_y[k] * gg_8[k];

        t_14[k] = -f_1 * gf_s_5[k]
                  + f_2 * gh_s_14[k]
                  + f_0 * gf_5[k]
                  + pb_z[k] * gg_8[k];

        t_15[k] = pa_y[k] * fh_0[k]
                  + f_2 * gh_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_y, fg_0, fg_1, fh_3, fh_4, gh_s_16, \
                         gh_s_17, gh_s_18, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_16[k]
                  + pb_y[k] * gg_9[k];

        t_17[k] = f_6 * fg_1[k]
                  + pa_y[k] * fh_3[k]
                  + f_2 * gh_s_17[k];

        t_18[k] = pa_y[k] * fh_4[k]
                  + f_2 * gh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pb_x, fg_3, fg_10, fh_5, fh_8, gh_s_19, \
                         gh_s_21, gh_s_22, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * fg_3[k]
                  + pa_y[k] * fh_5[k]
                  + f_2 * gh_s_19[k];

        t_20[k] = pa_y[k] * fh_8[k]
                  + f_2 * gh_s_21[k];

        t_21[k] = f_7 * fg_10[k]
                  + f_2 * gh_s_22[k]
                  + pb_x[k] * gg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pb_z, dh_s_14, dh_14, fh_18, gf_s_6, gh_s_23, \
                         gh_s_24, gh_s_25, gf_6, gg_10, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_8 * dh_s_14[k]
                  + f_6 * dh_14[k]
                  + pa_x[k] * fh_18[k]
                  + f_2 * gh_s_23[k];

        t_23[k] = f_2 * gh_s_24[k]
                  + pb_z[k] * gg_10[k];

        t_24[k] = -f_3 * gf_s_6[k]
                  + f_2 * gh_s_25[k]
                  + f_4 * gf_6[k]
                  + pb_z[k] * gg_11[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pb_y, pb_z, fg_8, fh_12, gf_s_7, gh_s_26, \
                         gh_s_27, gh_s_28, gf_7, gg_12, gg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_5 * gf_s_7[k]
                  + f_2 * gh_s_26[k]
                  + f_6 * gf_7[k]
                  + pb_z[k] * gg_12[k];

        t_26[k] = f_4 * fg_8[k]
                  + f_2 * gh_s_27[k]
                  + pb_y[k] * gg_13[k];

        t_27[k] = pa_y[k] * fh_12[k]
                  + f_2 * gh_s_28[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_z, pb_z, fg_0, fg_2, fh_0, fh_4, gh_s_29, \
                         gh_s_30, gh_s_32, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_z[k] * fh_0[k]
                  + f_2 * gh_s_29[k];

        t_29[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_30[k]
                  + pb_z[k] * gg_14[k];

        t_30[k] = f_6 * fg_2[k]
                  + pa_z[k] * fh_4[k]
                  + f_2 * gh_s_32[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_z, pb_x, pb_y, fg_4, fg_19, fh_8, gf_s_10, \
                         gh_s_35, gh_s_36, gh_s_37, gf_8, gg_15, \
                         gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * fg_4[k]
                  + pa_z[k] * fh_8[k]
                  + f_2 * gh_s_35[k];

        t_32[k] = f_7 * fg_19[k]
                  + f_2 * gh_s_36[k]
                  + pb_x[k] * gg_18[k];

        t_33[k] = -f_9 * gf_s_10[k]
                  + f_2 * gh_s_37[k]
                  + f_7 * gf_8[k]
                  + pb_y[k] * gg_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_y, gf_s_11, gf_s_12, gh_s_38, gh_s_39, gh_s_40, \
                         gf_9, gf_10, gg_16, gg_17, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -f_5 * gf_s_11[k]
                  + f_2 * gh_s_38[k]
                  + f_6 * gf_9[k]
                  + pb_y[k] * gg_16[k];

        t_35[k] = -f_3 * gf_s_12[k]
                  + f_2 * gh_s_39[k]
                  + f_4 * gf_10[k]
                  + pb_y[k] * gg_17[k];

        t_36[k] = f_2 * gh_s_40[k]
                  + pb_y[k] * gg_18[k];
    }

#pragma omp simd aligned(t_37, t_38, pa_x, pa_y, dh_s_0, dh_s_19, dh_0, dh_19, fh_13, fh_31, \
                         gh_s_41, gh_s_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_8 * dh_s_19[k]
                  + f_6 * dh_19[k]
                  + pa_x[k] * fh_31[k]
                  + f_2 * gh_s_41[k];

        t_38[k] = -f_10 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_13[k]
                  + f_2 * gh_s_42[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, pb_z, fg_9, fg_21, gf_s_15, gh_s_43, \
                         gh_s_44, gh_s_45, gf_13, gg_19, gg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_6 * fg_9[k]
                  + f_2 * gh_s_43[k]
                  + pb_y[k] * gg_19[k];

        t_40[k] = f_2 * gh_s_44[k]
                  + pb_z[k] * gg_19[k];

        t_41[k] = f_6 * fg_21[k]
                  - f_5 * gf_s_15[k]
                  + f_2 * gh_s_45[k]
                  + f_6 * gf_13[k]
                  + pb_x[k] * gg_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_z, fg_22, gf_s_13, gf_s_16, gh_s_46, \
                         gh_s_47, gh_s_48, gf_11, gf_14, gg_20, gg_21, \
                         gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -f_3 * gf_s_13[k]
                  + f_2 * gh_s_46[k]
                  + f_4 * gf_11[k]
                  + pb_z[k] * gg_20[k];

        t_43[k] = f_6 * fg_22[k]
                  - f_3 * gf_s_16[k]
                  + f_2 * gh_s_47[k]
                  + f_4 * gf_14[k]
                  + pb_x[k] * gg_23[k];

        t_44[k] = f_2 * gh_s_48[k]
                  + pb_z[k] * gg_21[k];
    }

#pragma omp simd aligned(t_45, t_46, pb_x, pb_z, fg_23, gf_s_14, gh_s_49, gh_s_50, gf_12, \
                         gg_22, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_5 * gf_s_14[k]
                  + f_2 * gh_s_49[k]
                  + f_6 * gf_12[k]
                  + pb_z[k] * gg_22[k];

        t_46[k] = f_6 * fg_23[k]
                  + f_2 * gh_s_50[k]
                  + pb_x[k] * gg_24[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pb_z, dh_s_28, dh_28, fh_38, gf_s_16, \
                         gh_s_51, gh_s_52, gh_s_53, gf_14, gg_24, \
                         gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_10 * dh_s_28[k]
                  + f_4 * dh_28[k]
                  + pa_x[k] * fh_38[k]
                  + f_2 * gh_s_51[k];

        t_48[k] = f_2 * gh_s_52[k]
                  + pb_z[k] * gg_24[k];

        t_49[k] = -f_3 * gf_s_16[k]
                  + f_2 * gh_s_53[k]
                  + f_4 * gf_14[k]
                  + pb_z[k] * gg_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_y, pb_z, fg_13, gf_s_17, gf_s_18, gh_s_54, \
                         gh_s_55, gh_s_56, gf_15, gf_16, gg_26, gg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_5 * gf_s_17[k]
                  + f_2 * gh_s_54[k]
                  + f_6 * gf_15[k]
                  + pb_z[k] * gg_26[k];

        t_51[k] = f_6 * fg_13[k]
                  + f_2 * gh_s_55[k]
                  + pb_y[k] * gg_27[k];

        t_52[k] = -f_1 * gf_s_18[k]
                  + f_2 * gh_s_56[k]
                  + f_0 * gf_16[k]
                  + pb_z[k] * gg_27[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, fh_14, fh_16, fh_24, fh_25, \
                         gh_s_57, gh_s_58, gh_s_59, gh_s_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_y[k] * fh_24[k]
                  + f_2 * gh_s_57[k];

        t_54[k] = pa_z[k] * fh_14[k]
                  + f_2 * gh_s_58[k];

        t_55[k] = pa_y[k] * fh_25[k]
                  + f_2 * gh_s_59[k];

        t_56[k] = pa_z[k] * fh_16[k]
                  + f_2 * gh_s_60[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pa_z, pb_z, fg_10, fh_18, fh_27, gh_s_61, \
                         gh_s_62, gh_s_63, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_y[k] * fh_27[k]
                  + f_2 * gh_s_61[k];

        t_58[k] = pa_z[k] * fh_18[k]
                  + f_2 * gh_s_62[k];

        t_59[k] = f_4 * fg_10[k]
                  + f_2 * gh_s_63[k]
                  + pb_z[k] * gg_28[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pb_y, dh_s_36, dh_s_37, dh_36, dh_37, fg_19, \
                         fh_49, fh_50, gh_s_64, gh_s_65, gh_s_66, \
                         gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -f_10 * dh_s_36[k]
                  + f_4 * dh_36[k]
                  + pa_x[k] * fh_49[k]
                  + f_2 * gh_s_64[k];

        t_61[k] = -f_10 * dh_s_37[k]
                  + f_4 * dh_37[k]
                  + pa_x[k] * fh_50[k]
                  + f_2 * gh_s_65[k];

        t_62[k] = f_4 * fg_19[k]
                  + f_2 * gh_s_66[k]
                  + pb_y[k] * gg_29[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pb_y, dh_s_0, dh_0, fh_23, fh_31, \
                         gh_s_67, gh_s_68, gh_s_69, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * fh_31[k]
                  + f_2 * gh_s_67[k];

        t_64[k] = -f_10 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_z[k] * fh_23[k]
                  + f_2 * gh_s_68[k];

        t_65[k] = f_2 * gh_s_69[k]
                  + pb_y[k] * gg_30[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_y, pb_z, fg_14, gf_s_19, gh_s_70, gh_s_71, \
                         gh_s_72, gf_17, gg_30, gg_31, gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_6 * fg_14[k]
                  + f_2 * gh_s_70[k]
                  + pb_z[k] * gg_30[k];

        t_67[k] = -f_3 * gf_s_19[k]
                  + f_2 * gh_s_71[k]
                  + f_4 * gf_17[k]
                  + pb_y[k] * gg_31[k];

        t_68[k] = f_2 * gh_s_72[k]
                  + pb_y[k] * gg_32[k];
    }

#pragma omp simd aligned(t_69, t_70, pb_x, pb_y, fg_25, gf_s_20, gf_s_22, gh_s_73, gh_s_74, \
                         gf_18, gf_20, gg_33, gg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * fg_25[k]
                  - f_5 * gf_s_22[k]
                  + f_2 * gh_s_73[k]
                  + f_6 * gf_20[k]
                  + pb_x[k] * gg_35[k];

        t_70[k] = -f_5 * gf_s_20[k]
                  + f_2 * gh_s_74[k]
                  + f_6 * gf_18[k]
                  + pb_y[k] * gg_33[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pb_x, pb_y, fg_26, gf_s_21, gf_s_26, gh_s_75, \
                         gh_s_76, gh_s_77, gf_19, gf_24, gg_34, gg_35, \
                         gg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -f_3 * gf_s_21[k]
                  + f_2 * gh_s_75[k]
                  + f_4 * gf_19[k]
                  + pb_y[k] * gg_34[k];

        t_72[k] = f_2 * gh_s_76[k]
                  + pb_y[k] * gg_35[k];

        t_73[k] = f_6 * fg_26[k]
                  - f_3 * gf_s_26[k]
                  + f_2 * gh_s_77[k]
                  + f_4 * gf_24[k]
                  + pb_x[k] * gg_36[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pb_x, pb_y, fg_27, gf_s_23, gf_s_24, gh_s_78, \
                         gh_s_79, gh_s_80, gf_21, gf_22, gg_37, gg_38, \
                         gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_6 * fg_27[k]
                  + f_2 * gh_s_78[k]
                  + pb_x[k] * gg_41[k];

        t_75[k] = -f_1 * gf_s_23[k]
                  + f_2 * gh_s_79[k]
                  + f_0 * gf_21[k]
                  + pb_y[k] * gg_37[k];

        t_76[k] = -f_9 * gf_s_24[k]
                  + f_2 * gh_s_80[k]
                  + f_7 * gf_22[k]
                  + pb_y[k] * gg_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, gf_s_25, gf_s_26, gh_s_81, gh_s_82, gh_s_83, \
                         gf_23, gf_24, gg_39, gg_40, gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -f_5 * gf_s_25[k]
                  + f_2 * gh_s_81[k]
                  + f_6 * gf_23[k]
                  + pb_y[k] * gg_39[k];

        t_78[k] = -f_3 * gf_s_26[k]
                  + f_2 * gh_s_82[k]
                  + f_4 * gf_24[k]
                  + pb_y[k] * gg_40[k];

        t_79[k] = f_2 * gh_s_83[k]
                  + pb_y[k] * gg_41[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_x, pb_y, dh_s_55, dh_55, fg_20, fg_28, fh_61, \
                         fh_62, gh_s_84, gh_s_85, gh_s_86, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -f_10 * dh_s_55[k]
                  + f_4 * dh_55[k]
                  + pa_x[k] * fh_61[k]
                  + f_2 * gh_s_84[k];

        t_81[k] = f_11 * fg_28[k]
                  + pa_x[k] * fh_62[k]
                  + f_2 * gh_s_85[k];

        t_82[k] = f_7 * fg_20[k]
                  + f_2 * gh_s_86[k]
                  + pb_y[k] * gg_42[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_x, fg_30, fg_31, fg_32, fh_64, fh_66, fh_67, \
                         gh_s_88, gh_s_89, gh_s_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_7 * fg_30[k]
                  + pa_x[k] * fh_64[k]
                  + f_2 * gh_s_88[k];

        t_84[k] = f_7 * fg_31[k]
                  + pa_x[k] * fh_66[k]
                  + f_2 * gh_s_89[k];

        t_85[k] = f_6 * fg_32[k]
                  + pa_x[k] * fh_67[k]
                  + f_2 * gh_s_90[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pb_x, fg_34, fg_35, fh_70, fh_75, \
                         fh_77, gh_s_92, gh_s_93, gh_s_94, gh_s_95, \
                         gg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_6 * fg_34[k]
                  + pa_x[k] * fh_70[k]
                  + f_2 * gh_s_92[k];

        t_87[k] = f_4 * fg_35[k]
                  + f_2 * gh_s_93[k]
                  + pb_x[k] * gg_43[k];

        t_88[k] = pa_x[k] * fh_75[k]
                  + f_2 * gh_s_94[k];

        t_89[k] = pa_x[k] * fh_77[k]
                  + f_2 * gh_s_95[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pa_z, fh_32, fh_78, fh_79, fh_80, \
                         gh_s_96, gh_s_97, gh_s_98, gh_s_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_x[k] * fh_78[k]
                  + f_2 * gh_s_96[k];

        t_91[k] = pa_x[k] * fh_79[k]
                  + f_2 * gh_s_97[k];

        t_92[k] = pa_x[k] * fh_80[k]
                  + f_2 * gh_s_98[k];

        t_93[k] = pa_z[k] * fh_32[k]
                  + f_2 * gh_s_99[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_z, pb_z, fg_20, fg_40, fh_33, fh_81, \
                         gh_s_100, gh_s_101, gh_s_102, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_4 * fg_20[k]
                  + f_2 * gh_s_100[k]
                  + pb_z[k] * gg_44[k];

        t_95[k] = pa_z[k] * fh_33[k]
                  + f_2 * gh_s_101[k];

        t_96[k] = f_7 * fg_40[k]
                  + pa_x[k] * fh_81[k]
                  + f_2 * gh_s_102[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pa_z, fg_41, fh_35, fh_82, fh_86, \
                         fh_87, gh_s_103, gh_s_104, gh_s_105, \
                         gh_s_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pa_z[k] * fh_35[k]
                  + f_2 * gh_s_103[k];

        t_98[k] = f_6 * fg_41[k]
                  + pa_x[k] * fh_82[k]
                  + f_2 * gh_s_104[k];

        t_99[k] = pa_x[k] * fh_86[k]
                  + f_2 * gh_s_105[k];

        t_100[k] = pa_x[k] * fh_87[k]
                   + f_2 * gh_s_106[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pa_y, fh_52, fh_88, fh_89, fh_90, \
                         gh_s_107, gh_s_108, gh_s_109, gh_s_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pa_x[k] * fh_88[k]
                   + f_2 * gh_s_107[k];

        t_102[k] = pa_x[k] * fh_89[k]
                   + f_2 * gh_s_108[k];

        t_103[k] = pa_x[k] * fh_90[k]
                   + f_2 * gh_s_109[k];

        t_104[k] = pa_y[k] * fh_52[k]
                   + f_2 * gh_s_110[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pa_y, fg_45, fg_46, fh_53, fh_54, \
                         fh_91, fh_92, gh_s_111, gh_s_112, gh_s_113, \
                         gh_s_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_y[k] * fh_53[k]
                   + f_2 * gh_s_111[k];

        t_106[k] = f_7 * fg_45[k]
                   + pa_x[k] * fh_91[k]
                   + f_2 * gh_s_112[k];

        t_107[k] = pa_y[k] * fh_54[k]
                   + f_2 * gh_s_113[k];

        t_108[k] = f_6 * fg_46[k]
                   + pa_x[k] * fh_92[k]
                   + f_2 * gh_s_114[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_y, fh_55, fh_95, fh_96, fh_97, \
                         gh_s_115, gh_s_116, gh_s_117, gh_s_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_y[k] * fh_55[k]
                   + f_2 * gh_s_115[k];

        t_110[k] = pa_x[k] * fh_95[k]
                   + f_2 * gh_s_116[k];

        t_111[k] = pa_x[k] * fh_96[k]
                   + f_2 * gh_s_117[k];

        t_112[k] = pa_x[k] * fh_97[k]
                   + f_2 * gh_s_118[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pb_z, fg_24, fg_51, fh_98, fh_99, \
                         fh_101, gh_s_119, gh_s_120, gh_s_121, gh_s_123, \
                         gg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_x[k] * fh_98[k]
                   + f_2 * gh_s_119[k];

        t_114[k] = pa_x[k] * fh_99[k]
                   + f_2 * gh_s_120[k];

        t_115[k] = f_11 * fg_51[k]
                   + pa_x[k] * fh_101[k]
                   + f_2 * gh_s_121[k];

        t_116[k] = f_7 * fg_24[k]
                   + f_2 * gh_s_123[k]
                   + pb_z[k] * gg_45[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_x, pb_x, fg_54, fg_57, fg_62, fh_106, fh_110, \
                         gh_s_126, gh_s_130, gh_s_131, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_7 * fg_54[k]
                   + pa_x[k] * fh_106[k]
                   + f_2 * gh_s_126[k];

        t_118[k] = f_6 * fg_57[k]
                   + pa_x[k] * fh_110[k]
                   + f_2 * gh_s_130[k];

        t_119[k] = f_4 * fg_62[k]
                   + f_2 * gh_s_131[k]
                   + pb_x[k] * gg_46[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pa_x, fh_115, fh_116, fh_117, \
                         fh_118, fh_120, gh_s_132, gh_s_133, gh_s_134, gh_s_135, \
                         gh_s_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * fh_115[k]
                   + f_2 * gh_s_132[k];

        t_121[k] = pa_x[k] * fh_116[k]
                   + f_2 * gh_s_133[k];

        t_122[k] = pa_x[k] * fh_117[k]
                   + f_2 * gh_s_134[k];

        t_123[k] = pa_x[k] * fh_118[k]
                   + f_2 * gh_s_135[k];

        t_124[k] = pa_x[k] * fh_120[k]
                   + f_2 * gh_s_136[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pb_x, gf_s_32, gf_s_33, gf_s_34, gh_s_137, \
                         gh_s_138, gh_s_139, gf_25, gf_26, gf_27, gg_47, gg_48, \
                         gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -f_1 * gf_s_32[k]
                   + f_2 * gh_s_137[k]
                   + f_0 * gf_25[k]
                   + pb_x[k] * gg_47[k];

        t_126[k] = -f_9 * gf_s_33[k]
                   + f_2 * gh_s_138[k]
                   + f_7 * gf_26[k]
                   + pb_x[k] * gg_48[k];

        t_127[k] = -f_5 * gf_s_34[k]
                   + f_2 * gh_s_139[k]
                   + f_6 * gf_27[k]
                   + pb_x[k] * gg_49[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, gf_s_35, gf_s_36, gf_s_38, gh_s_140, \
                         gh_s_141, gh_s_142, gf_28, gf_29, gf_31, gg_50, gg_51, \
                         gg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = -f_5 * gf_s_35[k]
                   + f_2 * gh_s_140[k]
                   + f_6 * gf_28[k]
                   + pb_x[k] * gg_50[k];

        t_129[k] = -f_3 * gf_s_36[k]
                   + f_2 * gh_s_141[k]
                   + f_4 * gf_29[k]
                   + pb_x[k] * gg_51[k];

        t_130[k] = -f_3 * gf_s_38[k]
                   + f_2 * gh_s_142[k]
                   + f_4 * gf_31[k]
                   + pb_x[k] * gg_52[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pb_x, gf_s_39, gh_s_143, gh_s_144, \
                         gh_s_145, gh_s_146, gf_32, gg_53, gg_54, gg_56, \
                         gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -f_3 * gf_s_39[k]
                   + f_2 * gh_s_143[k]
                   + f_4 * gf_32[k]
                   + pb_x[k] * gg_53[k];

        t_132[k] = f_2 * gh_s_144[k]
                   + pb_x[k] * gg_54[k];

        t_133[k] = f_2 * gh_s_145[k]
                   + pb_x[k] * gg_56[k];

        t_134[k] = f_2 * gh_s_146[k]
                   + pb_x[k] * gg_57[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pb_x, pb_y, pb_z, fg_35, gf_s_36, gh_s_147, \
                         gh_s_148, gh_s_149, gf_29, gg_54, gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_2 * gh_s_147[k]
                   + pb_x[k] * gg_58[k];

        t_136[k] = f_0 * fg_35[k]
                   - f_1 * gf_s_36[k]
                   + f_2 * gh_s_148[k]
                   + f_0 * gf_29[k]
                   + pb_y[k] * gg_54[k];

        t_137[k] = f_2 * gh_s_149[k]
                   + pb_z[k] * gg_54[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_y, pb_z, fg_39, gf_s_36, gf_s_37, gh_s_150, \
                         gh_s_151, gh_s_152, gf_29, gf_30, gg_55, gg_56, \
                         gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = -f_3 * gf_s_36[k]
                   + f_2 * gh_s_150[k]
                   + f_4 * gf_29[k]
                   + pb_z[k] * gg_55[k];

        t_139[k] = -f_5 * gf_s_37[k]
                   + f_2 * gh_s_151[k]
                   + f_6 * gf_30[k]
                   + pb_z[k] * gg_56[k];

        t_140[k] = f_0 * fg_39[k]
                   + f_2 * gh_s_152[k]
                   + pb_y[k] * gg_58[k];
    }

#pragma omp simd aligned(t_141, t_142, pb_x, pb_z, gf_s_39, gf_s_40, gh_s_153, gh_s_154, \
                         gf_32, gf_33, gg_58, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -f_1 * gf_s_39[k]
                   + f_2 * gh_s_153[k]
                   + f_0 * gf_32[k]
                   + pb_z[k] * gg_58[k];

        t_142[k] = -f_9 * gf_s_40[k]
                   + f_2 * gh_s_154[k]
                   + f_7 * gf_33[k]
                   + pb_x[k] * gg_59[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, gf_s_42, gf_s_46, gh_s_156, gh_s_159, \
                         gh_s_163, gf_34, gf_35, gg_60, gg_61, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = -f_5 * gf_s_42[k]
                   + f_2 * gh_s_156[k]
                   + f_6 * gf_34[k]
                   + pb_x[k] * gg_60[k];

        t_144[k] = -f_3 * gf_s_46[k]
                   + f_2 * gh_s_159[k]
                   + f_4 * gf_35[k]
                   + pb_x[k] * gg_61[k];

        t_145[k] = f_2 * gh_s_163[k]
                   + pb_x[k] * gg_63[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_z, pb_z, fg_35, fg_36, fh_75, fh_77, \
                         gh_s_164, gh_s_165, gh_s_166, gg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_z[k] * fh_75[k]
                   + f_2 * gh_s_164[k];

        t_147[k] = f_4 * fg_35[k]
                   + f_2 * gh_s_165[k]
                   + pb_z[k] * gg_62[k];

        t_148[k] = f_6 * fg_36[k]
                   + pa_z[k] * fh_77[k]
                   + f_2 * gh_s_166[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_y, pa_z, pb_y, dh_s_39, dh_39, fg_37, fg_44, \
                         fh_78, fh_90, gh_s_167, gh_s_168, gh_s_169, \
                         gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_7 * fg_37[k]
                   + pa_z[k] * fh_78[k]
                   + f_2 * gh_s_167[k];

        t_150[k] = f_7 * fg_44[k]
                   + f_2 * gh_s_168[k]
                   + pb_y[k] * gg_63[k];

        t_151[k] = -f_8 * dh_s_39[k]
                   + f_6 * dh_39[k]
                   + pa_y[k] * fh_90[k]
                   + f_2 * gh_s_169[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_x, gf_s_47, gf_s_48, gf_s_49, gh_s_170, \
                         gh_s_171, gh_s_172, gf_36, gf_37, gf_38, gg_64, gg_65, \
                         gg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -f_1 * gf_s_47[k]
                   + f_2 * gh_s_170[k]
                   + f_0 * gf_36[k]
                   + pb_x[k] * gg_64[k];

        t_153[k] = -f_9 * gf_s_48[k]
                   + f_2 * gh_s_171[k]
                   + f_7 * gf_37[k]
                   + pb_x[k] * gg_65[k];

        t_154[k] = -f_9 * gf_s_49[k]
                   + f_2 * gh_s_172[k]
                   + f_7 * gf_38[k]
                   + pb_x[k] * gg_66[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pb_x, gf_s_50, gf_s_51, gf_s_52, gh_s_173, \
                         gh_s_174, gh_s_175, gf_39, gf_40, gf_41, gg_67, gg_68, \
                         gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -f_5 * gf_s_50[k]
                   + f_2 * gh_s_173[k]
                   + f_6 * gf_39[k]
                   + pb_x[k] * gg_67[k];

        t_156[k] = -f_5 * gf_s_51[k]
                   + f_2 * gh_s_174[k]
                   + f_6 * gf_40[k]
                   + pb_x[k] * gg_68[k];

        t_157[k] = -f_5 * gf_s_52[k]
                   + f_2 * gh_s_175[k]
                   + f_6 * gf_41[k]
                   + pb_x[k] * gg_69[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pb_x, gf_s_53, gf_s_54, gf_s_55, gh_s_176, \
                         gh_s_177, gh_s_178, gf_42, gf_43, gf_44, gg_70, gg_71, \
                         gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = -f_3 * gf_s_53[k]
                   + f_2 * gh_s_176[k]
                   + f_4 * gf_42[k]
                   + pb_x[k] * gg_70[k];

        t_159[k] = -f_3 * gf_s_54[k]
                   + f_2 * gh_s_177[k]
                   + f_4 * gf_43[k]
                   + pb_x[k] * gg_71[k];

        t_160[k] = -f_3 * gf_s_55[k]
                   + f_2 * gh_s_178[k]
                   + f_4 * gf_44[k]
                   + pb_x[k] * gg_72[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_x, gf_s_56, gh_s_179, gh_s_180, \
                         gh_s_181, gh_s_182, gf_45, gg_73, gg_74, gg_75, \
                         gg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -f_3 * gf_s_56[k]
                   + f_2 * gh_s_179[k]
                   + f_4 * gf_45[k]
                   + pb_x[k] * gg_73[k];

        t_162[k] = f_2 * gh_s_180[k]
                   + pb_x[k] * gg_74[k];

        t_163[k] = f_2 * gh_s_181[k]
                   + pb_x[k] * gg_75[k];

        t_164[k] = f_2 * gh_s_182[k]
                   + pb_x[k] * gg_76[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_z, pb_x, dh_s_28, dh_28, fh_85, gh_s_183, \
                         gh_s_184, gh_s_185, gg_77, gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_2 * gh_s_183[k]
                   + pb_x[k] * gg_77[k];

        t_166[k] = f_2 * gh_s_184[k]
                   + pb_x[k] * gg_78[k];

        t_167[k] = -f_10 * dh_s_28[k]
                   + f_4 * dh_28[k]
                   + pa_z[k] * fh_85[k]
                   + f_2 * gh_s_185[k];
    }

#pragma omp simd aligned(t_168, t_169, pb_y, pb_z, fg_42, fg_48, gf_s_55, gh_s_186, gh_s_187, \
                         gf_44, gg_74, gg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_6 * fg_42[k]
                   + f_2 * gh_s_186[k]
                   + pb_z[k] * gg_74[k];

        t_169[k] = f_6 * fg_48[k]
                   - f_5 * gf_s_55[k]
                   + f_2 * gh_s_187[k]
                   + f_6 * gf_44[k]
                   + pb_y[k] * gg_76[k];
    }

#pragma omp simd aligned(t_170, t_171, pb_y, fg_49, fg_50, gf_s_56, gh_s_188, gh_s_189, gf_45, \
                         gg_77, gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_6 * fg_49[k]
                   - f_3 * gf_s_56[k]
                   + f_2 * gh_s_188[k]
                   + f_4 * gf_45[k]
                   + pb_y[k] * gg_77[k];

        t_171[k] = f_6 * fg_50[k]
                   + f_2 * gh_s_189[k]
                   + pb_y[k] * gg_78[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pa_y, pb_z, dh_s_55, dh_55, fg_47, fg_58, \
                         fh_100, fh_115, gh_s_190, gh_s_201, gh_s_202, \
                         gg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -f_10 * dh_s_55[k]
                   + f_4 * dh_55[k]
                   + pa_y[k] * fh_100[k]
                   + f_2 * gh_s_190[k];

        t_173[k] = f_11 * fg_58[k]
                   + pa_y[k] * fh_115[k]
                   + f_2 * gh_s_201[k];

        t_174[k] = f_7 * fg_47[k]
                   + f_2 * gh_s_202[k]
                   + pb_z[k] * gg_79[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_y, pb_y, fg_60, fg_61, fg_62, fh_117, fh_118, \
                         gh_s_203, gh_s_204, gh_s_205, gg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_7 * fg_60[k]
                   + pa_y[k] * fh_117[k]
                   + f_2 * gh_s_203[k];

        t_176[k] = f_6 * fg_61[k]
                   + pa_y[k] * fh_118[k]
                   + f_2 * gh_s_204[k];

        t_177[k] = f_4 * fg_62[k]
                   + f_2 * gh_s_205[k]
                   + pb_y[k] * gg_80[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pa_y, pb_x, fh_120, gf_s_64, gf_s_65, gh_s_206, \
                         gh_s_207, gh_s_208, gf_46, gf_47, gg_81, \
                         gg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = pa_y[k] * fh_120[k]
                   + f_2 * gh_s_206[k];

        t_179[k] = -f_1 * gf_s_64[k]
                   + f_2 * gh_s_207[k]
                   + f_0 * gf_46[k]
                   + pb_x[k] * gg_81[k];

        t_180[k] = -f_9 * gf_s_65[k]
                   + f_2 * gh_s_208[k]
                   + f_7 * gf_47[k]
                   + pb_x[k] * gg_82[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, pb_x, gf_s_66, gf_s_67, gf_s_68, gh_s_209, \
                         gh_s_210, gh_s_211, gf_48, gf_49, gf_50, gg_83, gg_84, \
                         gg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = -f_5 * gf_s_66[k]
                   + f_2 * gh_s_209[k]
                   + f_6 * gf_48[k]
                   + pb_x[k] * gg_83[k];

        t_182[k] = -f_5 * gf_s_67[k]
                   + f_2 * gh_s_210[k]
                   + f_6 * gf_49[k]
                   + pb_x[k] * gg_84[k];

        t_183[k] = -f_3 * gf_s_68[k]
                   + f_2 * gh_s_211[k]
                   + f_4 * gf_50[k]
                   + pb_x[k] * gg_85[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pb_x, gf_s_69, gf_s_71, gh_s_212, gh_s_213, \
                         gh_s_214, gf_51, gf_53, gg_86, gg_87, gg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -f_3 * gf_s_69[k]
                   + f_2 * gh_s_212[k]
                   + f_4 * gf_51[k]
                   + pb_x[k] * gg_86[k];

        t_185[k] = -f_3 * gf_s_71[k]
                   + f_2 * gh_s_213[k]
                   + f_4 * gf_53[k]
                   + pb_x[k] * gg_87[k];

        t_186[k] = f_2 * gh_s_214[k]
                   + pb_x[k] * gg_88[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_x, pb_y, gf_s_68, gh_s_215, gh_s_216, \
                         gh_s_217, gh_s_218, gf_50, gg_88, gg_89, gg_90, \
                         gg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_2 * gh_s_215[k]
                   + pb_x[k] * gg_89[k];

        t_188[k] = f_2 * gh_s_216[k]
                   + pb_x[k] * gg_90[k];

        t_189[k] = f_2 * gh_s_217[k]
                   + pb_x[k] * gg_92[k];

        t_190[k] = -f_1 * gf_s_68[k]
                   + f_2 * gh_s_218[k]
                   + f_0 * gf_50[k]
                   + pb_y[k] * gg_88[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pb_y, gf_s_69, gf_s_70, gf_s_71, gh_s_219, \
                         gh_s_220, gh_s_221, gf_51, gf_52, gf_53, gg_89, gg_90, \
                         gg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = -f_9 * gf_s_69[k]
                   + f_2 * gh_s_219[k]
                   + f_7 * gf_51[k]
                   + pb_y[k] * gg_89[k];

        t_192[k] = -f_5 * gf_s_70[k]
                   + f_2 * gh_s_220[k]
                   + f_6 * gf_52[k]
                   + pb_y[k] * gg_90[k];

        t_193[k] = -f_3 * gf_s_71[k]
                   + f_2 * gh_s_221[k]
                   + f_4 * gf_53[k]
                   + pb_y[k] * gg_91[k];
    }

#pragma omp simd aligned(t_194, t_195, pb_y, pb_z, fg_62, gf_s_71, gh_s_222, gh_s_223, gf_53, \
                         gg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_2 * gh_s_222[k]
                   + pb_y[k] * gg_92[k];

        t_195[k] = f_0 * fg_62[k]
                   - f_1 * gf_s_71[k]
                   + f_2 * gh_s_223[k]
                   + f_0 * gf_53[k]
                   + pb_z[k] * gg_92[k];
    }
}

auto
compute_prim_gh_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = beta / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_14 = buffer.data(dh_s + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_14 = buffer.data(dh + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);

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

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_3 = buffer.data(gf_s + 3);
    const auto *gf_s_5 = buffer.data(gf_s + 5);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_17 = buffer.data(gf_s + 17);
    const auto *gf_s_22 = buffer.data(gf_s + 22);
    const auto *gf_s_23 = buffer.data(gf_s + 23);
    const auto *gf_s_24 = buffer.data(gf_s + 24);
    const auto *gf_s_28 = buffer.data(gf_s + 28);
    const auto *gf_s_39 = buffer.data(gf_s + 39);
    const auto *gf_s_40 = buffer.data(gf_s + 40);
    const auto *gf_s_41 = buffer.data(gf_s + 41);
    const auto *gf_s_42 = buffer.data(gf_s + 42);
    const auto *gf_s_43 = buffer.data(gf_s + 43);
    const auto *gf_s_44 = buffer.data(gf_s + 44);
    const auto *gf_s_45 = buffer.data(gf_s + 45);
    const auto *gf_s_46 = buffer.data(gf_s + 46);
    const auto *gf_s_47 = buffer.data(gf_s + 47);
    const auto *gf_s_49 = buffer.data(gf_s + 49);
    const auto *gf_s_50 = buffer.data(gf_s + 50);
    const auto *gf_s_51 = buffer.data(gf_s + 51);
    const auto *gf_s_52 = buffer.data(gf_s + 52);
    const auto *gf_s_53 = buffer.data(gf_s + 53);
    const auto *gf_s_54 = buffer.data(gf_s + 54);
    const auto *gf_s_55 = buffer.data(gf_s + 55);
    const auto *gf_s_60 = buffer.data(gf_s + 60);
    const auto *gf_s_61 = buffer.data(gf_s + 61);
    const auto *gf_s_62 = buffer.data(gf_s + 62);
    const auto *gf_s_63 = buffer.data(gf_s + 63);
    const auto *gf_s_64 = buffer.data(gf_s + 64);
    const auto *gf_s_65 = buffer.data(gf_s + 65);
    const auto *gf_s_66 = buffer.data(gf_s + 66);
    const auto *gf_s_67 = buffer.data(gf_s + 67);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_20 = buffer.data(gh_s + 20);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_31 = buffer.data(gh_s + 31);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_33 = buffer.data(gh_s + 33);
    const auto *gh_s_34 = buffer.data(gh_s + 34);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);
    const auto *gh_s_58 = buffer.data(gh_s + 58);
    const auto *gh_s_59 = buffer.data(gh_s + 59);
    const auto *gh_s_60 = buffer.data(gh_s + 60);
    const auto *gh_s_61 = buffer.data(gh_s + 61);
    const auto *gh_s_62 = buffer.data(gh_s + 62);
    const auto *gh_s_63 = buffer.data(gh_s + 63);
    const auto *gh_s_64 = buffer.data(gh_s + 64);
    const auto *gh_s_65 = buffer.data(gh_s + 65);
    const auto *gh_s_66 = buffer.data(gh_s + 66);
    const auto *gh_s_67 = buffer.data(gh_s + 67);
    const auto *gh_s_68 = buffer.data(gh_s + 68);
    const auto *gh_s_69 = buffer.data(gh_s + 69);
    const auto *gh_s_70 = buffer.data(gh_s + 70);
    const auto *gh_s_71 = buffer.data(gh_s + 71);
    const auto *gh_s_72 = buffer.data(gh_s + 72);
    const auto *gh_s_73 = buffer.data(gh_s + 73);
    const auto *gh_s_74 = buffer.data(gh_s + 74);
    const auto *gh_s_75 = buffer.data(gh_s + 75);
    const auto *gh_s_76 = buffer.data(gh_s + 76);
    const auto *gh_s_77 = buffer.data(gh_s + 77);
    const auto *gh_s_78 = buffer.data(gh_s + 78);
    const auto *gh_s_79 = buffer.data(gh_s + 79);
    const auto *gh_s_80 = buffer.data(gh_s + 80);
    const auto *gh_s_81 = buffer.data(gh_s + 81);
    const auto *gh_s_82 = buffer.data(gh_s + 82);
    const auto *gh_s_83 = buffer.data(gh_s + 83);
    const auto *gh_s_84 = buffer.data(gh_s + 84);
    const auto *gh_s_85 = buffer.data(gh_s + 85);
    const auto *gh_s_86 = buffer.data(gh_s + 86);
    const auto *gh_s_87 = buffer.data(gh_s + 87);
    const auto *gh_s_88 = buffer.data(gh_s + 88);
    const auto *gh_s_89 = buffer.data(gh_s + 89);
    const auto *gh_s_90 = buffer.data(gh_s + 90);
    const auto *gh_s_91 = buffer.data(gh_s + 91);
    const auto *gh_s_92 = buffer.data(gh_s + 92);
    const auto *gh_s_93 = buffer.data(gh_s + 93);
    const auto *gh_s_94 = buffer.data(gh_s + 94);
    const auto *gh_s_95 = buffer.data(gh_s + 95);
    const auto *gh_s_96 = buffer.data(gh_s + 96);
    const auto *gh_s_97 = buffer.data(gh_s + 97);
    const auto *gh_s_98 = buffer.data(gh_s + 98);
    const auto *gh_s_99 = buffer.data(gh_s + 99);
    const auto *gh_s_100 = buffer.data(gh_s + 100);
    const auto *gh_s_101 = buffer.data(gh_s + 101);
    const auto *gh_s_102 = buffer.data(gh_s + 102);
    const auto *gh_s_103 = buffer.data(gh_s + 103);
    const auto *gh_s_104 = buffer.data(gh_s + 104);
    const auto *gh_s_105 = buffer.data(gh_s + 105);
    const auto *gh_s_106 = buffer.data(gh_s + 106);
    const auto *gh_s_107 = buffer.data(gh_s + 107);
    const auto *gh_s_108 = buffer.data(gh_s + 108);
    const auto *gh_s_109 = buffer.data(gh_s + 109);
    const auto *gh_s_110 = buffer.data(gh_s + 110);
    const auto *gh_s_111 = buffer.data(gh_s + 111);
    const auto *gh_s_112 = buffer.data(gh_s + 112);
    const auto *gh_s_113 = buffer.data(gh_s + 113);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gf_0, gg_0, gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_1[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];

        t_2[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_2[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, gf_s_1, gf_s_2, gh_s_3, gh_s_4, gh_s_5, \
                         gf_1, gf_2, gg_3, gg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_3[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_4[k] = f_2 * gh_s_4[k]
                 + pb_z[k] * gg_3[k];

        t_5[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_5[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, fg_5, fg_6, gf_s_3, gh_s_6, gh_s_7, \
                         gh_s_8, gf_3, gg_5, gg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * fg_5[k]
                 + f_2 * gh_s_6[k]
                 + pb_x[k] * gg_5[k];

        t_7[k] = f_0 * fg_6[k]
                 + f_2 * gh_s_7[k]
                 + pb_x[k] * gg_7[k];

        t_8[k] = -f_1 * gf_s_3[k]
                 + f_2 * gh_s_8[k]
                 + f_0 * gf_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_y, pb_z, fg_0, fh_0, gf_s_5, gh_s_9, \
                         gh_s_10, gh_s_11, gf_5, gg_7, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_1 * gf_s_5[k]
                 + f_2 * gh_s_9[k]
                 + f_0 * gf_5[k]
                 + pb_z[k] * gg_7[k];

        t_10[k] = pa_y[k] * fh_0[k]
                  + f_2 * gh_s_10[k];

        t_11[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_11[k]
                  + pb_y[k] * gg_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_x, fg_1, fg_3, fg_9, fh_1, fh_3, gh_s_12, \
                         gh_s_13, gh_s_14, gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_6 * fg_1[k]
                  + pa_y[k] * fh_1[k]
                  + f_2 * gh_s_12[k];

        t_13[k] = f_7 * fg_3[k]
                  + pa_y[k] * fh_3[k]
                  + f_2 * gh_s_13[k];

        t_14[k] = f_7 * fg_9[k]
                  + f_2 * gh_s_14[k]
                  + pb_x[k] * gg_11[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_z, pb_z, dh_s_3, dh_3, fg_0, fh_0, fh_6, \
                         gh_s_15, gh_s_16, gh_s_17, gg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_8 * dh_s_3[k]
                  + f_6 * dh_3[k]
                  + pa_x[k] * fh_6[k]
                  + f_2 * gh_s_15[k];

        t_16[k] = pa_z[k] * fh_0[k]
                  + f_2 * gh_s_16[k];

        t_17[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_17[k]
                  + pb_z[k] * gg_15[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_x, fg_2, fg_4, fg_13, fh_2, fh_4, gh_s_18, \
                         gh_s_19, gh_s_20, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_6 * fg_2[k]
                  + pa_z[k] * fh_2[k]
                  + f_2 * gh_s_18[k];

        t_19[k] = f_7 * fg_4[k]
                  + pa_z[k] * fh_4[k]
                  + f_2 * gh_s_19[k];

        t_20[k] = f_7 * fg_13[k]
                  + f_2 * gh_s_20[k]
                  + pb_x[k] * gg_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_x, pa_y, dh_s_0, dh_s_4, dh_0, dh_4, fh_5, fh_10, \
                         gh_s_21, gh_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_8 * dh_s_4[k]
                  + f_6 * dh_4[k]
                  + pa_x[k] * fh_10[k]
                  + f_2 * gh_s_21[k];

        t_22[k] = -f_9 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_5[k]
                  + f_2 * gh_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pb_y, fg_7, fg_15, gf_s_16, gh_s_23, gh_s_24, \
                         gf_16, gg_21, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * fg_7[k]
                  + f_2 * gh_s_23[k]
                  + pb_y[k] * gg_21[k];

        t_24[k] = f_6 * fg_15[k]
                  - f_5 * gf_s_16[k]
                  + f_2 * gh_s_24[k]
                  + f_6 * gf_16[k]
                  + pb_x[k] * gg_22[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, fg_16, fg_17, gf_s_17, gh_s_25, gh_s_26, gf_17, \
                         gg_24, gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * fg_16[k]
                  - f_3 * gf_s_17[k]
                  + f_2 * gh_s_25[k]
                  + f_4 * gf_17[k]
                  + pb_x[k] * gg_24[k];

        t_26[k] = f_6 * fg_17[k]
                  + f_2 * gh_s_26[k]
                  + pb_x[k] * gg_25[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, dh_s_5, dh_5, fh_8, fh_9, fh_11, \
                         gh_s_27, gh_s_28, gh_s_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_9 * dh_s_5[k]
                  + f_4 * dh_5[k]
                  + pa_x[k] * fh_11[k]
                  + f_2 * gh_s_27[k];

        t_28[k] = pa_y[k] * fh_8[k]
                  + f_2 * gh_s_28[k];

        t_29[k] = pa_y[k] * fh_9[k]
                  + f_2 * gh_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, pa_x, dh_s_8, dh_s_9, dh_8, dh_9, fh_12, fh_13, gh_s_30, \
                         gh_s_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_9 * dh_s_8[k]
                  + f_4 * dh_8[k]
                  + pa_x[k] * fh_12[k]
                  + f_2 * gh_s_30[k];

        t_31[k] = -f_9 * dh_s_9[k]
                  + f_4 * dh_9[k]
                  + pa_x[k] * fh_13[k]
                  + f_2 * gh_s_31[k];
    }

#pragma omp simd aligned(t_32, t_33, pa_z, pb_z, dh_s_0, dh_0, fg_10, fh_7, gh_s_32, gh_s_33, \
                         gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_9 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_z[k] * fh_7[k]
                  + f_2 * gh_s_32[k];

        t_33[k] = f_6 * fg_10[k]
                  + f_2 * gh_s_33[k]
                  + pb_z[k] * gg_37[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_x, pb_y, fg_21, gf_s_22, gf_s_24, gh_s_34, gh_s_35, \
                         gf_22, gf_24, gg_38, gg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -f_3 * gf_s_22[k]
                  + f_2 * gh_s_34[k]
                  + f_4 * gf_22[k]
                  + pb_y[k] * gg_38[k];

        t_35[k] = f_6 * fg_21[k]
                  - f_5 * gf_s_24[k]
                  + f_2 * gh_s_35[k]
                  + f_6 * gf_24[k]
                  + pb_x[k] * gg_41[k];
    }

#pragma omp simd aligned(t_36, t_37, pb_x, pb_y, fg_22, gf_s_23, gf_s_28, gh_s_36, gh_s_37, \
                         gf_23, gf_28, gg_40, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_5 * gf_s_23[k]
                  + f_2 * gh_s_36[k]
                  + f_6 * gf_23[k]
                  + pb_y[k] * gg_40[k];

        t_37[k] = f_6 * fg_22[k]
                  - f_3 * gf_s_28[k]
                  + f_2 * gh_s_37[k]
                  + f_4 * gf_28[k]
                  + pb_x[k] * gg_42[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_x, dh_s_14, dh_14, fg_23, fg_24, fh_14, \
                         fh_15, gh_s_38, gh_s_39, gh_s_40, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * fg_23[k]
                  + f_2 * gh_s_38[k]
                  + pb_x[k] * gg_46[k];

        t_39[k] = -f_9 * dh_s_14[k]
                  + f_4 * dh_14[k]
                  + pa_x[k] * fh_14[k]
                  + f_2 * gh_s_39[k];

        t_40[k] = f_10 * fg_24[k]
                  + pa_x[k] * fh_15[k]
                  + f_2 * gh_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_x, pb_y, fg_14, fg_25, fg_27, fh_16, fh_17, \
                         gh_s_41, gh_s_42, gh_s_43, gg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * fg_14[k]
                  + f_2 * gh_s_41[k]
                  + pb_y[k] * gg_47[k];

        t_42[k] = f_7 * fg_25[k]
                  + pa_x[k] * fh_16[k]
                  + f_2 * gh_s_42[k];

        t_43[k] = f_6 * fg_27[k]
                  + pa_x[k] * fh_17[k]
                  + f_2 * gh_s_43[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, fg_28, fh_18, fh_22, fh_23, \
                         gh_s_44, gh_s_45, gh_s_46, gh_s_47, gg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * fg_28[k]
                  + f_2 * gh_s_44[k]
                  + pb_x[k] * gg_51[k];

        t_45[k] = pa_x[k] * fh_18[k]
                  + f_2 * gh_s_45[k];

        t_46[k] = pa_x[k] * fh_22[k]
                  + f_2 * gh_s_46[k];

        t_47[k] = pa_x[k] * fh_23[k]
                  + f_2 * gh_s_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, fh_24, fh_25, fh_26, fh_27, gh_s_48, \
                         gh_s_49, gh_s_50, gh_s_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_x[k] * fh_24[k]
                  + f_2 * gh_s_48[k];

        t_49[k] = pa_x[k] * fh_25[k]
                  + f_2 * gh_s_49[k];

        t_50[k] = pa_x[k] * fh_26[k]
                  + f_2 * gh_s_50[k];

        t_51[k] = pa_x[k] * fh_27[k]
                  + f_2 * gh_s_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_z, fg_19, fg_45, fg_48, fh_29, fh_31, \
                         gh_s_52, gh_s_53, gh_s_54, gg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * fg_45[k]
                  + pa_x[k] * fh_29[k]
                  + f_2 * gh_s_52[k];

        t_53[k] = f_7 * fg_19[k]
                  + f_2 * gh_s_53[k]
                  + pb_z[k] * gg_71[k];

        t_54[k] = f_7 * fg_48[k]
                  + pa_x[k] * fh_31[k]
                  + f_2 * gh_s_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_x, pb_x, fg_50, fg_55, fh_33, fh_37, gh_s_55, \
                         gh_s_56, gh_s_57, gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_6 * fg_50[k]
                  + pa_x[k] * fh_33[k]
                  + f_2 * gh_s_55[k];

        t_56[k] = f_4 * fg_55[k]
                  + f_2 * gh_s_56[k]
                  + pb_x[k] * gg_78[k];

        t_57[k] = pa_x[k] * fh_37[k]
                  + f_2 * gh_s_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, gf_s_39, gf_s_40, gf_s_41, gh_s_58, gh_s_59, \
                         gh_s_60, gf_38, gf_39, gf_40, gg_79, gg_80, \
                         gg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -f_1 * gf_s_39[k]
                  + f_2 * gh_s_58[k]
                  + f_0 * gf_38[k]
                  + pb_x[k] * gg_79[k];

        t_59[k] = -f_11 * gf_s_40[k]
                  + f_2 * gh_s_59[k]
                  + f_7 * gf_39[k]
                  + pb_x[k] * gg_80[k];

        t_60[k] = -f_5 * gf_s_41[k]
                  + f_2 * gh_s_60[k]
                  + f_6 * gf_40[k]
                  + pb_x[k] * gg_81[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, pb_z, gf_s_42, gf_s_43, gh_s_61, gh_s_62, \
                         gh_s_63, gf_41, gf_42, gg_80, gg_83, gg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_2 * gh_s_61[k]
                  + pb_z[k] * gg_80[k];

        t_62[k] = -f_5 * gf_s_42[k]
                  + f_2 * gh_s_62[k]
                  + f_6 * gf_41[k]
                  + pb_x[k] * gg_83[k];

        t_63[k] = -f_3 * gf_s_43[k]
                  + f_2 * gh_s_63[k]
                  + f_4 * gf_42[k]
                  + pb_x[k] * gg_84[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_z, gf_s_45, gf_s_46, gh_s_64, gh_s_65, \
                         gh_s_66, gf_44, gf_45, gg_81, gg_85, gg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * gh_s_64[k]
                  + pb_z[k] * gg_81[k];

        t_65[k] = -f_3 * gf_s_45[k]
                  + f_2 * gh_s_65[k]
                  + f_4 * gf_44[k]
                  + pb_x[k] * gg_85[k];

        t_66[k] = -f_3 * gf_s_46[k]
                  + f_2 * gh_s_66[k]
                  + f_4 * gf_45[k]
                  + pb_x[k] * gg_86[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_y, pb_z, fg_28, gf_s_43, gf_s_44, gh_s_67, \
                         gh_s_68, gh_s_69, gf_42, gf_43, gg_87, gg_88, \
                         gg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * fg_28[k]
                  - f_1 * gf_s_43[k]
                  + f_2 * gh_s_67[k]
                  + f_0 * gf_42[k]
                  + pb_y[k] * gg_87[k];

        t_68[k] = -f_3 * gf_s_43[k]
                  + f_2 * gh_s_68[k]
                  + f_4 * gf_42[k]
                  + pb_z[k] * gg_88[k];

        t_69[k] = -f_5 * gf_s_44[k]
                  + f_2 * gh_s_69[k]
                  + f_6 * gf_43[k]
                  + pb_z[k] * gg_89[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_x, pb_y, pb_z, fg_32, gf_s_46, gf_s_47, gh_s_70, \
                         gh_s_71, gh_s_72, gf_45, gf_46, gg_91, gg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * fg_32[k]
                  + f_2 * gh_s_70[k]
                  + pb_y[k] * gg_91[k];

        t_71[k] = -f_1 * gf_s_46[k]
                  + f_2 * gh_s_71[k]
                  + f_0 * gf_45[k]
                  + pb_z[k] * gg_91[k];

        t_72[k] = -f_5 * gf_s_47[k]
                  + f_2 * gh_s_72[k]
                  + f_6 * gf_46[k]
                  + pb_x[k] * gg_92[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_z, pb_x, pb_z, fg_28, fh_18, gf_s_49, gh_s_73, \
                         gh_s_74, gh_s_75, gf_48, gg_93, gg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -f_3 * gf_s_49[k]
                  + f_2 * gh_s_73[k]
                  + f_4 * gf_48[k]
                  + pb_x[k] * gg_93[k];

        t_74[k] = pa_z[k] * fh_18[k]
                  + f_2 * gh_s_74[k];

        t_75[k] = f_4 * fg_28[k]
                  + f_2 * gh_s_75[k]
                  + pb_z[k] * gg_94[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pb_y, fg_29, fg_30, fg_38, fh_19, fh_20, \
                         gh_s_76, gh_s_77, gh_s_78, gg_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_6 * fg_29[k]
                  + pa_z[k] * fh_19[k]
                  + f_2 * gh_s_76[k];

        t_77[k] = f_7 * fg_30[k]
                  + pa_z[k] * fh_20[k]
                  + f_2 * gh_s_77[k];

        t_78[k] = f_7 * fg_38[k]
                  + f_2 * gh_s_78[k]
                  + pb_y[k] * gg_98[k];
    }

#pragma omp simd aligned(t_79, t_80, pa_y, pb_x, dh_s_10, dh_10, fh_24, gf_s_50, gh_s_79, \
                         gh_s_80, gf_49, gg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -f_8 * dh_s_10[k]
                  + f_6 * dh_10[k]
                  + pa_y[k] * fh_24[k]
                  + f_2 * gh_s_79[k];

        t_80[k] = -f_1 * gf_s_50[k]
                  + f_2 * gh_s_80[k]
                  + f_0 * gf_49[k]
                  + pb_x[k] * gg_99[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, gf_s_51, gf_s_52, gf_s_53, gh_s_81, gh_s_82, \
                         gh_s_83, gf_50, gf_51, gf_52, gg_100, gg_101, \
                         gg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_5 * gf_s_51[k]
                  + f_2 * gh_s_81[k]
                  + f_6 * gf_50[k]
                  + pb_x[k] * gg_100[k];

        t_82[k] = -f_5 * gf_s_52[k]
                  + f_2 * gh_s_82[k]
                  + f_6 * gf_51[k]
                  + pb_x[k] * gg_101[k];

        t_83[k] = -f_3 * gf_s_53[k]
                  + f_2 * gh_s_83[k]
                  + f_4 * gf_52[k]
                  + pb_x[k] * gg_102[k];
    }

#pragma omp simd aligned(t_84, t_85, pa_z, pb_x, dh_s_5, dh_5, fh_21, gf_s_55, gh_s_84, \
                         gh_s_85, gf_54, gg_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_3 * gf_s_55[k]
                  + f_2 * gh_s_84[k]
                  + f_4 * gf_54[k]
                  + pb_x[k] * gg_103[k];

        t_85[k] = -f_9 * dh_s_5[k]
                  + f_4 * dh_5[k]
                  + pa_z[k] * fh_21[k]
                  + f_2 * gh_s_85[k];
    }

#pragma omp simd aligned(t_86, t_87, pb_y, pb_z, fg_34, fg_42, gf_s_54, gh_s_86, gh_s_87, \
                         gf_53, gg_104, gg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_6 * fg_34[k]
                  + f_2 * gh_s_86[k]
                  + pb_z[k] * gg_104[k];

        t_87[k] = f_6 * fg_42[k]
                  - f_5 * gf_s_54[k]
                  + f_2 * gh_s_87[k]
                  + f_6 * gf_53[k]
                  + pb_y[k] * gg_106[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, fg_43, fg_44, gf_s_55, gh_s_88, gh_s_89, gf_54, \
                         gg_107, gg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_6 * fg_43[k]
                  - f_3 * gf_s_55[k]
                  + f_2 * gh_s_88[k]
                  + f_4 * gf_54[k]
                  + pb_y[k] * gg_107[k];

        t_89[k] = f_6 * fg_44[k]
                  + f_2 * gh_s_89[k]
                  + pb_y[k] * gg_108[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_y, dh_s_14, dh_14, fg_46, fg_47, fh_28, fh_30, \
                         fh_32, gh_s_90, gh_s_91, gh_s_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -f_9 * dh_s_14[k]
                  + f_4 * dh_14[k]
                  + pa_y[k] * fh_28[k]
                  + f_2 * gh_s_90[k];

        t_91[k] = f_6 * fg_46[k]
                  + pa_y[k] * fh_30[k]
                  + f_2 * gh_s_91[k];

        t_92[k] = f_7 * fg_47[k]
                  + pa_y[k] * fh_32[k]
                  + f_2 * gh_s_92[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_y, pb_z, fg_40, fg_51, fg_53, fh_34, fh_35, \
                         gh_s_93, gh_s_94, gh_s_95, gg_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_10 * fg_51[k]
                  + pa_y[k] * fh_34[k]
                  + f_2 * gh_s_93[k];

        t_94[k] = f_7 * fg_40[k]
                  + f_2 * gh_s_94[k]
                  + pb_z[k] * gg_111[k];

        t_95[k] = f_7 * fg_53[k]
                  + pa_y[k] * fh_35[k]
                  + f_2 * gh_s_95[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_y, pb_y, fg_54, fg_55, fh_36, fh_37, gh_s_96, \
                         gh_s_97, gh_s_98, gg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_6 * fg_54[k]
                  + pa_y[k] * fh_36[k]
                  + f_2 * gh_s_96[k];

        t_97[k] = f_4 * fg_55[k]
                  + f_2 * gh_s_97[k]
                  + pb_y[k] * gg_115[k];

        t_98[k] = pa_y[k] * fh_37[k]
                  + f_2 * gh_s_98[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_y, gf_s_60, gf_s_61, gh_s_99, gh_s_100, \
                         gh_s_101, gf_59, gf_60, gg_116, gg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -f_1 * gf_s_60[k]
                  + f_2 * gh_s_99[k]
                  + f_0 * gf_59[k]
                  + pb_x[k] * gg_116[k];

        t_100[k] = f_2 * gh_s_100[k]
                   + pb_y[k] * gg_116[k];

        t_101[k] = -f_11 * gf_s_61[k]
                   + f_2 * gh_s_101[k]
                   + f_7 * gf_60[k]
                   + pb_x[k] * gg_118[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_x, pb_y, gf_s_62, gf_s_63, gh_s_102, \
                         gh_s_103, gh_s_104, gf_61, gf_62, gg_118, gg_119, \
                         gg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -f_5 * gf_s_62[k]
                   + f_2 * gh_s_102[k]
                   + f_6 * gf_61[k]
                   + pb_x[k] * gg_119[k];

        t_103[k] = f_2 * gh_s_103[k]
                   + pb_y[k] * gg_118[k];

        t_104[k] = -f_5 * gf_s_63[k]
                   + f_2 * gh_s_104[k]
                   + f_6 * gf_62[k]
                   + pb_x[k] * gg_121[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pb_x, pb_y, gf_s_64, gf_s_65, gh_s_105, \
                         gh_s_106, gh_s_107, gf_63, gf_64, gg_121, gg_122, \
                         gg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -f_3 * gf_s_64[k]
                   + f_2 * gh_s_105[k]
                   + f_4 * gf_63[k]
                   + pb_x[k] * gg_122[k];

        t_106[k] = -f_3 * gf_s_65[k]
                   + f_2 * gh_s_106[k]
                   + f_4 * gf_64[k]
                   + pb_x[k] * gg_123[k];

        t_107[k] = f_2 * gh_s_107[k]
                   + pb_y[k] * gg_121[k];
    }

#pragma omp simd aligned(t_108, t_109, pb_x, pb_y, gf_s_64, gf_s_67, gh_s_108, gh_s_109, \
                         gf_63, gf_66, gg_124, gg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = -f_3 * gf_s_67[k]
                   + f_2 * gh_s_108[k]
                   + f_4 * gf_66[k]
                   + pb_x[k] * gg_124[k];

        t_109[k] = -f_1 * gf_s_64[k]
                   + f_2 * gh_s_109[k]
                   + f_0 * gf_63[k]
                   + pb_y[k] * gg_125[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_y, gf_s_65, gf_s_66, gf_s_67, gh_s_110, \
                         gh_s_111, gh_s_112, gf_64, gf_65, gf_66, gg_126, gg_127, \
                         gg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -f_11 * gf_s_65[k]
                   + f_2 * gh_s_110[k]
                   + f_7 * gf_64[k]
                   + pb_y[k] * gg_126[k];

        t_111[k] = -f_5 * gf_s_66[k]
                   + f_2 * gh_s_111[k]
                   + f_6 * gf_65[k]
                   + pb_y[k] * gg_127[k];

        t_112[k] = -f_3 * gf_s_67[k]
                   + f_2 * gh_s_112[k]
                   + f_4 * gf_66[k]
                   + pb_y[k] * gg_128[k];
    }

#pragma omp simd aligned(t_113, pb_z, fg_55, gf_s_67, gh_s_113, gf_66, \
                         gg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_0 * fg_55[k]
                   - f_1 * gf_s_67[k]
                   + f_2 * gh_s_113[k]
                   + f_0 * gf_66[k]
                   + pb_z[k] * gg_129[k];
    }
}

auto
compute_prim_gh_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_27 = buffer.data(dh_s + 27);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_27 = buffer.data(dh + 27);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_6 = buffer.data(fg + 6);
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
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
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
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_38 = buffer.data(fh + 38);
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
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_78 = buffer.data(fh + 78);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_3 = buffer.data(gf_s + 3);
    const auto *gf_s_4 = buffer.data(gf_s + 4);
    const auto *gf_s_5 = buffer.data(gf_s + 5);
    const auto *gf_s_7 = buffer.data(gf_s + 7);
    const auto *gf_s_8 = buffer.data(gf_s + 8);
    const auto *gf_s_10 = buffer.data(gf_s + 10);
    const auto *gf_s_11 = buffer.data(gf_s + 11);
    const auto *gf_s_12 = buffer.data(gf_s + 12);
    const auto *gf_s_13 = buffer.data(gf_s + 13);
    const auto *gf_s_14 = buffer.data(gf_s + 14);
    const auto *gf_s_15 = buffer.data(gf_s + 15);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_17 = buffer.data(gf_s + 17);
    const auto *gf_s_18 = buffer.data(gf_s + 18);
    const auto *gf_s_19 = buffer.data(gf_s + 19);
    const auto *gf_s_20 = buffer.data(gf_s + 20);
    const auto *gf_s_21 = buffer.data(gf_s + 21);
    const auto *gf_s_22 = buffer.data(gf_s + 22);
    const auto *gf_s_23 = buffer.data(gf_s + 23);
    const auto *gf_s_24 = buffer.data(gf_s + 24);
    const auto *gf_s_25 = buffer.data(gf_s + 25);
    const auto *gf_s_33 = buffer.data(gf_s + 33);
    const auto *gf_s_34 = buffer.data(gf_s + 34);
    const auto *gf_s_35 = buffer.data(gf_s + 35);
    const auto *gf_s_36 = buffer.data(gf_s + 36);
    const auto *gf_s_37 = buffer.data(gf_s + 37);
    const auto *gf_s_38 = buffer.data(gf_s + 38);
    const auto *gf_s_39 = buffer.data(gf_s + 39);
    const auto *gf_s_40 = buffer.data(gf_s + 40);
    const auto *gf_s_41 = buffer.data(gf_s + 41);
    const auto *gf_s_43 = buffer.data(gf_s + 43);
    const auto *gf_s_44 = buffer.data(gf_s + 44);
    const auto *gf_s_45 = buffer.data(gf_s + 45);
    const auto *gf_s_46 = buffer.data(gf_s + 46);
    const auto *gf_s_47 = buffer.data(gf_s + 47);
    const auto *gf_s_48 = buffer.data(gf_s + 48);
    const auto *gf_s_49 = buffer.data(gf_s + 49);
    const auto *gf_s_54 = buffer.data(gf_s + 54);
    const auto *gf_s_55 = buffer.data(gf_s + 55);
    const auto *gf_s_56 = buffer.data(gf_s + 56);
    const auto *gf_s_57 = buffer.data(gf_s + 57);
    const auto *gf_s_58 = buffer.data(gf_s + 58);
    const auto *gf_s_59 = buffer.data(gf_s + 59);
    const auto *gf_s_60 = buffer.data(gf_s + 60);
    const auto *gf_s_61 = buffer.data(gf_s + 61);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_20 = buffer.data(gh_s + 20);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_31 = buffer.data(gh_s + 31);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_33 = buffer.data(gh_s + 33);
    const auto *gh_s_34 = buffer.data(gh_s + 34);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);
    const auto *gh_s_58 = buffer.data(gh_s + 58);
    const auto *gh_s_59 = buffer.data(gh_s + 59);
    const auto *gh_s_60 = buffer.data(gh_s + 60);
    const auto *gh_s_61 = buffer.data(gh_s + 61);
    const auto *gh_s_62 = buffer.data(gh_s + 62);
    const auto *gh_s_63 = buffer.data(gh_s + 63);
    const auto *gh_s_64 = buffer.data(gh_s + 64);
    const auto *gh_s_65 = buffer.data(gh_s + 65);
    const auto *gh_s_66 = buffer.data(gh_s + 66);
    const auto *gh_s_67 = buffer.data(gh_s + 67);
    const auto *gh_s_68 = buffer.data(gh_s + 68);
    const auto *gh_s_69 = buffer.data(gh_s + 69);
    const auto *gh_s_70 = buffer.data(gh_s + 70);
    const auto *gh_s_71 = buffer.data(gh_s + 71);
    const auto *gh_s_72 = buffer.data(gh_s + 72);
    const auto *gh_s_73 = buffer.data(gh_s + 73);
    const auto *gh_s_74 = buffer.data(gh_s + 74);
    const auto *gh_s_75 = buffer.data(gh_s + 75);
    const auto *gh_s_76 = buffer.data(gh_s + 76);
    const auto *gh_s_77 = buffer.data(gh_s + 77);
    const auto *gh_s_78 = buffer.data(gh_s + 78);
    const auto *gh_s_79 = buffer.data(gh_s + 79);
    const auto *gh_s_80 = buffer.data(gh_s + 80);
    const auto *gh_s_81 = buffer.data(gh_s + 81);
    const auto *gh_s_82 = buffer.data(gh_s + 82);
    const auto *gh_s_83 = buffer.data(gh_s + 83);
    const auto *gh_s_84 = buffer.data(gh_s + 84);
    const auto *gh_s_85 = buffer.data(gh_s + 85);
    const auto *gh_s_86 = buffer.data(gh_s + 86);
    const auto *gh_s_87 = buffer.data(gh_s + 87);
    const auto *gh_s_88 = buffer.data(gh_s + 88);
    const auto *gh_s_89 = buffer.data(gh_s + 89);
    const auto *gh_s_90 = buffer.data(gh_s + 90);
    const auto *gh_s_91 = buffer.data(gh_s + 91);
    const auto *gh_s_92 = buffer.data(gh_s + 92);
    const auto *gh_s_93 = buffer.data(gh_s + 93);
    const auto *gh_s_94 = buffer.data(gh_s + 94);
    const auto *gh_s_95 = buffer.data(gh_s + 95);
    const auto *gh_s_96 = buffer.data(gh_s + 96);
    const auto *gh_s_97 = buffer.data(gh_s + 97);
    const auto *gh_s_98 = buffer.data(gh_s + 98);
    const auto *gh_s_99 = buffer.data(gh_s + 99);
    const auto *gh_s_100 = buffer.data(gh_s + 100);
    const auto *gh_s_101 = buffer.data(gh_s + 101);
    const auto *gh_s_102 = buffer.data(gh_s + 102);
    const auto *gh_s_103 = buffer.data(gh_s + 103);
    const auto *gh_s_104 = buffer.data(gh_s + 104);
    const auto *gh_s_105 = buffer.data(gh_s + 105);
    const auto *gh_s_106 = buffer.data(gh_s + 106);
    const auto *gh_s_107 = buffer.data(gh_s + 107);
    const auto *gh_s_108 = buffer.data(gh_s + 108);
    const auto *gh_s_109 = buffer.data(gh_s + 109);
    const auto *gh_s_110 = buffer.data(gh_s + 110);
    const auto *gh_s_111 = buffer.data(gh_s + 111);
    const auto *gh_s_112 = buffer.data(gh_s + 112);
    const auto *gh_s_113 = buffer.data(gh_s + 113);
    const auto *gh_s_114 = buffer.data(gh_s + 114);
    const auto *gh_s_115 = buffer.data(gh_s + 115);
    const auto *gh_s_116 = buffer.data(gh_s + 116);
    const auto *gh_s_117 = buffer.data(gh_s + 117);
    const auto *gh_s_118 = buffer.data(gh_s + 118);
    const auto *gh_s_119 = buffer.data(gh_s + 119);
    const auto *gh_s_120 = buffer.data(gh_s + 120);
    const auto *gh_s_121 = buffer.data(gh_s + 121);
    const auto *gh_s_122 = buffer.data(gh_s + 122);
    const auto *gh_s_123 = buffer.data(gh_s + 123);
    const auto *gh_s_124 = buffer.data(gh_s + 124);
    const auto *gh_s_125 = buffer.data(gh_s + 125);
    const auto *gh_s_126 = buffer.data(gh_s + 126);
    const auto *gh_s_127 = buffer.data(gh_s + 127);
    const auto *gh_s_128 = buffer.data(gh_s + 128);
    const auto *gh_s_129 = buffer.data(gh_s + 129);
    const auto *gh_s_130 = buffer.data(gh_s + 130);
    const auto *gh_s_131 = buffer.data(gh_s + 131);
    const auto *gh_s_132 = buffer.data(gh_s + 132);
    const auto *gh_s_133 = buffer.data(gh_s + 133);
    const auto *gh_s_134 = buffer.data(gh_s + 134);
    const auto *gh_s_135 = buffer.data(gh_s + 135);
    const auto *gh_s_136 = buffer.data(gh_s + 136);
    const auto *gh_s_137 = buffer.data(gh_s + 137);
    const auto *gh_s_138 = buffer.data(gh_s + 138);
    const auto *gh_s_139 = buffer.data(gh_s + 139);
    const auto *gh_s_140 = buffer.data(gh_s + 140);
    const auto *gh_s_141 = buffer.data(gh_s + 141);
    const auto *gh_s_142 = buffer.data(gh_s + 142);
    const auto *gh_s_143 = buffer.data(gh_s + 143);
    const auto *gh_s_144 = buffer.data(gh_s + 144);
    const auto *gh_s_145 = buffer.data(gh_s + 145);
    const auto *gh_s_146 = buffer.data(gh_s + 146);
    const auto *gh_s_147 = buffer.data(gh_s + 147);
    const auto *gh_s_148 = buffer.data(gh_s + 148);
    const auto *gh_s_149 = buffer.data(gh_s + 149);
    const auto *gh_s_150 = buffer.data(gh_s + 150);
    const auto *gh_s_151 = buffer.data(gh_s + 151);
    const auto *gh_s_152 = buffer.data(gh_s + 152);
    const auto *gh_s_153 = buffer.data(gh_s + 153);
    const auto *gh_s_154 = buffer.data(gh_s + 154);
    const auto *gh_s_155 = buffer.data(gh_s + 155);
    const auto *gh_s_156 = buffer.data(gh_s + 156);
    const auto *gh_s_157 = buffer.data(gh_s + 157);
    const auto *gh_s_158 = buffer.data(gh_s + 158);
    const auto *gh_s_159 = buffer.data(gh_s + 159);
    const auto *gh_s_160 = buffer.data(gh_s + 160);
    const auto *gh_s_161 = buffer.data(gh_s + 161);
    const auto *gh_s_162 = buffer.data(gh_s + 162);
    const auto *gh_s_163 = buffer.data(gh_s + 163);
    const auto *gh_s_164 = buffer.data(gh_s + 164);
    const auto *gh_s_165 = buffer.data(gh_s + 165);
    const auto *gh_s_166 = buffer.data(gh_s + 166);
    const auto *gh_s_167 = buffer.data(gh_s + 167);
    const auto *gh_s_168 = buffer.data(gh_s + 168);
    const auto *gh_s_169 = buffer.data(gh_s + 169);
    const auto *gh_s_170 = buffer.data(gh_s + 170);
    const auto *gh_s_171 = buffer.data(gh_s + 171);
    const auto *gh_s_172 = buffer.data(gh_s + 172);
    const auto *gh_s_173 = buffer.data(gh_s + 173);
    const auto *gh_s_174 = buffer.data(gh_s + 174);
    const auto *gh_s_175 = buffer.data(gh_s + 175);
    const auto *gh_s_176 = buffer.data(gh_s + 176);
    const auto *gh_s_177 = buffer.data(gh_s + 177);
    const auto *gh_s_178 = buffer.data(gh_s + 178);
    const auto *gh_s_179 = buffer.data(gh_s + 179);
    const auto *gh_s_180 = buffer.data(gh_s + 180);
    const auto *gh_s_181 = buffer.data(gh_s + 181);
    const auto *gh_s_182 = buffer.data(gh_s + 182);
    const auto *gh_s_183 = buffer.data(gh_s + 183);
    const auto *gh_s_184 = buffer.data(gh_s + 184);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);

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
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
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
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_48 = buffer.data(gg + 48);
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
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gh_s_3, gf_0, gg_0, gg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_2 * gh_s_1[k]
                 + pb_y[k] * gg_0[k];

        t_2[k] = f_2 * gh_s_2[k]
                 + pb_z[k] * gg_0[k];

        t_3[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_3[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, gf_s_0, gf_s_1, gh_s_4, gh_s_5, gh_s_6, \
                         gf_0, gf_1, gg_2, gg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_4[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];

        t_5[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_5[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = f_2 * gh_s_6[k]
                 + pb_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, gf_s_2, gf_s_3, gh_s_7, gh_s_8, gh_s_9, \
                         gf_2, gf_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * gh_s_7[k]
                 + pb_y[k] * gg_4[k];

        t_8[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_8[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_4[k];

        t_9[k] = -f_1 * gf_s_3[k]
                 + f_2 * gh_s_9[k]
                 + f_0 * gf_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, gf_s_4, gf_s_5, gh_s_10, gh_s_11, \
                         gh_s_12, gf_4, gf_5, gg_6, gg_7, gg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_5 * gf_s_4[k]
                  + f_2 * gh_s_10[k]
                  + f_6 * gf_4[k]
                  + pb_y[k] * gg_6[k];

        t_11[k] = -f_3 * gf_s_5[k]
                  + f_2 * gh_s_11[k]
                  + f_4 * gf_5[k]
                  + pb_y[k] * gg_7[k];

        t_12[k] = -f_1 * gf_s_5[k]
                  + f_2 * gh_s_12[k]
                  + f_0 * gf_5[k]
                  + pb_z[k] * gg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, fg_1, fg_3, fh_0, fh_3, fh_4, fh_5, \
                         gh_s_13, gh_s_14, gh_s_15, gh_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fh_0[k]
                  + f_2 * gh_s_13[k];

        t_14[k] = f_6 * fg_1[k]
                  + pa_y[k] * fh_3[k]
                  + f_2 * gh_s_14[k];

        t_15[k] = pa_y[k] * fh_4[k]
                  + f_2 * gh_s_15[k];

        t_16[k] = f_7 * fg_3[k]
                  + pa_y[k] * fh_5[k]
                  + f_2 * gh_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, dh_s_6, dh_6, fh_8, fh_14, \
                         gf_s_7, gh_s_17, gh_s_18, gh_s_19, gf_7, \
                         gg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_y[k] * fh_8[k]
                  + f_2 * gh_s_17[k];

        t_18[k] = -f_8 * dh_s_6[k]
                  + f_6 * dh_6[k]
                  + pa_x[k] * fh_14[k]
                  + f_2 * gh_s_18[k];

        t_19[k] = -f_3 * gf_s_7[k]
                  + f_2 * gh_s_19[k]
                  + f_4 * gf_7[k]
                  + pb_z[k] * gg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, fg_6, fh_10, gf_s_8, gh_s_20, \
                         gh_s_21, gh_s_22, gf_8, gg_12, gg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_5 * gf_s_8[k]
                  + f_2 * gh_s_20[k]
                  + f_6 * gf_8[k]
                  + pb_z[k] * gg_12[k];

        t_21[k] = f_4 * fg_6[k]
                  + f_2 * gh_s_21[k]
                  + pb_y[k] * gg_13[k];

        t_22[k] = pa_y[k] * fh_10[k]
                  + f_2 * gh_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_z, pb_z, fg_0, fg_2, fh_0, fh_4, gh_s_23, \
                         gh_s_24, gh_s_25, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_z[k] * fh_0[k]
                  + f_2 * gh_s_23[k];

        t_24[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_24[k]
                  + pb_z[k] * gg_14[k];

        t_25[k] = f_6 * fg_2[k]
                  + pa_z[k] * fh_4[k]
                  + f_2 * gh_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, fg_4, fh_8, gf_s_10, gh_s_26, gh_s_27, \
                         gh_s_28, gf_10, gg_15, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * gh_s_26[k]
                  + pb_y[k] * gg_15[k];

        t_27[k] = f_7 * fg_4[k]
                  + pa_z[k] * fh_8[k]
                  + f_2 * gh_s_27[k];

        t_28[k] = -f_9 * gf_s_10[k]
                  + f_2 * gh_s_28[k]
                  + f_7 * gf_10[k]
                  + pb_y[k] * gg_16[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_y, gf_s_11, gf_s_12, gh_s_29, gh_s_30, gf_11, gf_12, \
                         gg_17, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_5 * gf_s_11[k]
                  + f_2 * gh_s_29[k]
                  + f_6 * gf_11[k]
                  + pb_y[k] * gg_17[k];

        t_30[k] = -f_3 * gf_s_12[k]
                  + f_2 * gh_s_30[k]
                  + f_4 * gf_12[k]
                  + pb_y[k] * gg_18[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_x, pa_y, dh_s_0, dh_s_7, dh_0, dh_7, fh_11, fh_19, \
                         gh_s_31, gh_s_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_8 * dh_s_7[k]
                  + f_6 * dh_7[k]
                  + pa_x[k] * fh_19[k]
                  + f_2 * gh_s_31[k];

        t_32[k] = -f_10 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_11[k]
                  + f_2 * gh_s_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_z, fg_13, gf_s_13, gf_s_15, gh_s_33, gh_s_34, \
                         gf_13, gf_15, gg_21, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * fg_13[k]
                  - f_5 * gf_s_15[k]
                  + f_2 * gh_s_33[k]
                  + f_6 * gf_15[k]
                  + pb_x[k] * gg_22[k];

        t_34[k] = -f_3 * gf_s_13[k]
                  + f_2 * gh_s_34[k]
                  + f_4 * gf_13[k]
                  + pb_z[k] * gg_21[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pb_z, fg_14, gf_s_14, gf_s_16, gh_s_35, \
                         gh_s_36, gh_s_37, gf_14, gf_16, gg_22, gg_23, \
                         gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * fg_14[k]
                  - f_3 * gf_s_16[k]
                  + f_2 * gh_s_35[k]
                  + f_4 * gf_16[k]
                  + pb_x[k] * gg_24[k];

        t_36[k] = f_2 * gh_s_36[k]
                  + pb_z[k] * gg_22[k];

        t_37[k] = -f_5 * gf_s_14[k]
                  + f_2 * gh_s_37[k]
                  + f_6 * gf_14[k]
                  + pb_z[k] * gg_23[k];
    }

#pragma omp simd aligned(t_38, t_39, pa_x, pb_x, dh_s_11, dh_11, fg_15, fh_24, gh_s_38, \
                         gh_s_39, gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * fg_15[k]
                  + f_2 * gh_s_38[k]
                  + pb_x[k] * gg_25[k];

        t_39[k] = -f_10 * dh_s_11[k]
                  + f_4 * dh_11[k]
                  + pa_x[k] * fh_24[k]
                  + f_2 * gh_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pb_z, fg_9, gf_s_16, gf_s_17, gh_s_40, \
                         gh_s_41, gh_s_42, gf_16, gf_17, gg_26, gg_27, \
                         gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_3 * gf_s_16[k]
                  + f_2 * gh_s_40[k]
                  + f_4 * gf_16[k]
                  + pb_z[k] * gg_26[k];

        t_41[k] = -f_5 * gf_s_17[k]
                  + f_2 * gh_s_41[k]
                  + f_6 * gf_17[k]
                  + pb_z[k] * gg_27[k];

        t_42[k] = f_6 * fg_9[k]
                  + f_2 * gh_s_42[k]
                  + pb_y[k] * gg_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_z, fh_12, fh_16, gf_s_18, gh_s_43, \
                         gh_s_44, gh_s_45, gf_18, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_1 * gf_s_18[k]
                  + f_2 * gh_s_43[k]
                  + f_0 * gf_18[k]
                  + pb_z[k] * gg_28[k];

        t_44[k] = pa_y[k] * fh_16[k]
                  + f_2 * gh_s_44[k];

        t_45[k] = pa_z[k] * fh_12[k]
                  + f_2 * gh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, fh_13, fh_14, fh_17, fh_18, \
                         gh_s_46, gh_s_47, gh_s_48, gh_s_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * fh_17[k]
                  + f_2 * gh_s_46[k];

        t_47[k] = pa_z[k] * fh_13[k]
                  + f_2 * gh_s_47[k];

        t_48[k] = pa_y[k] * fh_18[k]
                  + f_2 * gh_s_48[k];

        t_49[k] = pa_z[k] * fh_14[k]
                  + f_2 * gh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_z, dh_s_15, dh_s_16, dh_15, dh_16, fg_8, \
                         fh_25, fh_26, gh_s_50, gh_s_51, gh_s_52, \
                         gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_4 * fg_8[k]
                  + f_2 * gh_s_50[k]
                  + pb_z[k] * gg_29[k];

        t_51[k] = -f_10 * dh_s_15[k]
                  + f_4 * dh_15[k]
                  + pa_x[k] * fh_25[k]
                  + f_2 * gh_s_51[k];

        t_52[k] = -f_10 * dh_s_16[k]
                  + f_4 * dh_16[k]
                  + pa_x[k] * fh_26[k]
                  + f_2 * gh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pa_z, pb_y, dh_s_0, dh_0, fg_11, fh_15, \
                         fh_19, gh_s_53, gh_s_54, gh_s_55, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * fg_11[k]
                  + f_2 * gh_s_53[k]
                  + pb_y[k] * gg_30[k];

        t_54[k] = pa_y[k] * fh_19[k]
                  + f_2 * gh_s_54[k];

        t_55[k] = -f_10 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_z[k] * fh_15[k]
                  + f_2 * gh_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_y, pb_z, fg_10, gf_s_19, gh_s_56, gh_s_57, \
                         gh_s_58, gf_19, gg_31, gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * gh_s_56[k]
                  + pb_y[k] * gg_31[k];

        t_57[k] = f_6 * fg_10[k]
                  + f_2 * gh_s_57[k]
                  + pb_z[k] * gg_31[k];

        t_58[k] = -f_3 * gf_s_19[k]
                  + f_2 * gh_s_58[k]
                  + f_4 * gf_19[k]
                  + pb_y[k] * gg_32[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_x, pb_y, fg_17, gf_s_20, gf_s_21, gh_s_59, \
                         gh_s_60, gh_s_61, gf_20, gf_21, gg_33, gg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_6 * fg_17[k]
                  - f_5 * gf_s_21[k]
                  + f_2 * gh_s_59[k]
                  + f_6 * gf_21[k]
                  + pb_x[k] * gg_34[k];

        t_60[k] = -f_5 * gf_s_20[k]
                  + f_2 * gh_s_60[k]
                  + f_6 * gf_20[k]
                  + pb_y[k] * gg_33[k];

        t_61[k] = f_2 * gh_s_61[k]
                  + pb_y[k] * gg_34[k];
    }

#pragma omp simd aligned(t_62, t_63, pb_x, fg_18, fg_19, gf_s_25, gh_s_62, gh_s_63, gf_25, \
                         gg_35, gg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_6 * fg_18[k]
                  - f_3 * gf_s_25[k]
                  + f_2 * gh_s_62[k]
                  + f_4 * gf_25[k]
                  + pb_x[k] * gg_35[k];

        t_63[k] = f_6 * fg_19[k]
                  + f_2 * gh_s_63[k]
                  + pb_x[k] * gg_40[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_y, gf_s_22, gf_s_23, gf_s_24, gh_s_64, gh_s_65, \
                         gh_s_66, gf_22, gf_23, gf_24, gg_36, gg_37, \
                         gg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = -f_1 * gf_s_22[k]
                  + f_2 * gh_s_64[k]
                  + f_0 * gf_22[k]
                  + pb_y[k] * gg_36[k];

        t_65[k] = -f_9 * gf_s_23[k]
                  + f_2 * gh_s_65[k]
                  + f_7 * gf_23[k]
                  + pb_y[k] * gg_37[k];

        t_66[k] = -f_5 * gf_s_24[k]
                  + f_2 * gh_s_66[k]
                  + f_6 * gf_24[k]
                  + pb_y[k] * gg_38[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_x, pb_y, dh_s_27, dh_27, fg_20, fh_32, fh_33, \
                         gf_s_25, gh_s_67, gh_s_68, gh_s_69, gf_25, \
                         gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -f_3 * gf_s_25[k]
                  + f_2 * gh_s_67[k]
                  + f_4 * gf_25[k]
                  + pb_y[k] * gg_39[k];

        t_68[k] = -f_10 * dh_s_27[k]
                  + f_4 * dh_27[k]
                  + pa_x[k] * fh_32[k]
                  + f_2 * gh_s_68[k];

        t_69[k] = f_11 * fg_20[k]
                  + pa_x[k] * fh_33[k]
                  + f_2 * gh_s_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_x, fg_21, fg_22, fg_23, fh_34, fh_35, fh_36, \
                         gh_s_70, gh_s_71, gh_s_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_7 * fg_21[k]
                  + pa_x[k] * fh_34[k]
                  + f_2 * gh_s_70[k];

        t_71[k] = f_7 * fg_22[k]
                  + pa_x[k] * fh_35[k]
                  + f_2 * gh_s_71[k];

        t_72[k] = f_6 * fg_23[k]
                  + pa_x[k] * fh_36[k]
                  + f_2 * gh_s_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, fg_24, fg_25, fh_38, fh_41, \
                         fh_43, gh_s_73, gh_s_74, gh_s_75, gh_s_76, \
                         gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * fg_24[k]
                  + pa_x[k] * fh_38[k]
                  + f_2 * gh_s_73[k];

        t_74[k] = f_4 * fg_25[k]
                  + f_2 * gh_s_74[k]
                  + pb_x[k] * gg_44[k];

        t_75[k] = pa_x[k] * fh_41[k]
                  + f_2 * gh_s_75[k];

        t_76[k] = pa_x[k] * fh_43[k]
                  + f_2 * gh_s_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_z, fh_20, fh_44, fh_45, fh_46, \
                         gh_s_77, gh_s_78, gh_s_79, gh_s_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_x[k] * fh_44[k]
                  + f_2 * gh_s_77[k];

        t_78[k] = pa_x[k] * fh_45[k]
                  + f_2 * gh_s_78[k];

        t_79[k] = pa_x[k] * fh_46[k]
                  + f_2 * gh_s_79[k];

        t_80[k] = pa_z[k] * fh_20[k]
                  + f_2 * gh_s_80[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, pa_z, pb_z, fg_12, fg_29, fh_21, fh_47, \
                         gh_s_81, gh_s_82, gh_s_83, gg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * fg_12[k]
                  + f_2 * gh_s_81[k]
                  + pb_z[k] * gg_45[k];

        t_82[k] = pa_z[k] * fh_21[k]
                  + f_2 * gh_s_82[k];

        t_83[k] = f_7 * fg_29[k]
                  + pa_x[k] * fh_47[k]
                  + f_2 * gh_s_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pa_z, fg_30, fh_22, fh_48, fh_50, \
                         fh_51, gh_s_84, gh_s_85, gh_s_86, gh_s_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_z[k] * fh_22[k]
                  + f_2 * gh_s_84[k];

        t_85[k] = f_6 * fg_30[k]
                  + pa_x[k] * fh_48[k]
                  + f_2 * gh_s_85[k];

        t_86[k] = pa_x[k] * fh_50[k]
                  + f_2 * gh_s_86[k];

        t_87[k] = pa_x[k] * fh_51[k]
                  + f_2 * gh_s_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, fh_27, fh_52, fh_53, fh_54, \
                         gh_s_88, gh_s_89, gh_s_90, gh_s_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pa_x[k] * fh_52[k]
                  + f_2 * gh_s_88[k];

        t_89[k] = pa_x[k] * fh_53[k]
                  + f_2 * gh_s_89[k];

        t_90[k] = pa_x[k] * fh_54[k]
                  + f_2 * gh_s_90[k];

        t_91[k] = pa_y[k] * fh_27[k]
                  + f_2 * gh_s_91[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pa_y, fg_33, fg_34, fh_28, fh_29, \
                         fh_55, fh_56, gh_s_92, gh_s_93, gh_s_94, \
                         gh_s_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_y[k] * fh_28[k]
                  + f_2 * gh_s_92[k];

        t_93[k] = f_7 * fg_33[k]
                  + pa_x[k] * fh_55[k]
                  + f_2 * gh_s_93[k];

        t_94[k] = pa_y[k] * fh_29[k]
                  + f_2 * gh_s_94[k];

        t_95[k] = f_6 * fg_34[k]
                  + pa_x[k] * fh_56[k]
                  + f_2 * gh_s_95[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pa_y, fh_30, fh_57, fh_58, fh_59, \
                         gh_s_96, gh_s_97, gh_s_98, gh_s_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_y[k] * fh_30[k]
                  + f_2 * gh_s_96[k];

        t_97[k] = pa_x[k] * fh_57[k]
                  + f_2 * gh_s_97[k];

        t_98[k] = pa_x[k] * fh_58[k]
                  + f_2 * gh_s_98[k];

        t_99[k] = pa_x[k] * fh_59[k]
                  + f_2 * gh_s_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_z, fg_16, fg_39, fh_60, fh_61, \
                         fh_63, gh_s_100, gh_s_101, gh_s_102, gh_s_103, \
                         gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * fh_60[k]
                   + f_2 * gh_s_100[k];

        t_101[k] = pa_x[k] * fh_61[k]
                   + f_2 * gh_s_101[k];

        t_102[k] = f_11 * fg_39[k]
                   + pa_x[k] * fh_63[k]
                   + f_2 * gh_s_102[k];

        t_103[k] = f_7 * fg_16[k]
                   + f_2 * gh_s_103[k]
                   + pb_z[k] * gg_48[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_x, pb_x, fg_42, fg_44, fg_49, fh_66, fh_69, \
                         gh_s_104, gh_s_105, gh_s_106, gg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_7 * fg_42[k]
                   + pa_x[k] * fh_66[k]
                   + f_2 * gh_s_104[k];

        t_105[k] = f_6 * fg_44[k]
                   + pa_x[k] * fh_69[k]
                   + f_2 * gh_s_105[k];

        t_106[k] = f_4 * fg_49[k]
                   + f_2 * gh_s_106[k]
                   + pb_x[k] * gg_51[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, pa_x, fh_73, fh_74, fh_75, fh_76, \
                         fh_78, gh_s_107, gh_s_108, gh_s_109, gh_s_110, \
                         gh_s_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * fh_73[k]
                   + f_2 * gh_s_107[k];

        t_108[k] = pa_x[k] * fh_74[k]
                   + f_2 * gh_s_108[k];

        t_109[k] = pa_x[k] * fh_75[k]
                   + f_2 * gh_s_109[k];

        t_110[k] = pa_x[k] * fh_76[k]
                   + f_2 * gh_s_110[k];

        t_111[k] = pa_x[k] * fh_78[k]
                   + f_2 * gh_s_111[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, gf_s_33, gf_s_34, gf_s_35, gh_s_112, \
                         gh_s_113, gh_s_114, gf_30, gf_31, gf_32, gg_52, gg_53, \
                         gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -f_1 * gf_s_33[k]
                   + f_2 * gh_s_112[k]
                   + f_0 * gf_30[k]
                   + pb_x[k] * gg_52[k];

        t_113[k] = -f_9 * gf_s_34[k]
                   + f_2 * gh_s_113[k]
                   + f_7 * gf_31[k]
                   + pb_x[k] * gg_53[k];

        t_114[k] = -f_5 * gf_s_35[k]
                   + f_2 * gh_s_114[k]
                   + f_6 * gf_32[k]
                   + pb_x[k] * gg_54[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_x, pb_z, gf_s_36, gf_s_37, gh_s_115, \
                         gh_s_116, gh_s_117, gf_33, gf_34, gg_53, gg_55, \
                         gg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_2 * gh_s_115[k]
                   + pb_z[k] * gg_53[k];

        t_116[k] = -f_5 * gf_s_36[k]
                   + f_2 * gh_s_116[k]
                   + f_6 * gf_33[k]
                   + pb_x[k] * gg_55[k];

        t_117[k] = -f_3 * gf_s_37[k]
                   + f_2 * gh_s_117[k]
                   + f_4 * gf_34[k]
                   + pb_x[k] * gg_56[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_z, gf_s_39, gf_s_40, gh_s_118, \
                         gh_s_119, gh_s_120, gf_36, gf_37, gg_54, gg_57, \
                         gg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_2 * gh_s_118[k]
                   + pb_z[k] * gg_54[k];

        t_119[k] = -f_3 * gf_s_39[k]
                   + f_2 * gh_s_119[k]
                   + f_4 * gf_36[k]
                   + pb_x[k] * gg_57[k];

        t_120[k] = -f_3 * gf_s_40[k]
                   + f_2 * gh_s_120[k]
                   + f_4 * gf_37[k]
                   + pb_x[k] * gg_58[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, gh_s_121, gh_s_122, gh_s_123, \
                         gh_s_124, gg_59, gg_61, gg_62, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_2 * gh_s_121[k]
                   + pb_x[k] * gg_59[k];

        t_122[k] = f_2 * gh_s_122[k]
                   + pb_x[k] * gg_61[k];

        t_123[k] = f_2 * gh_s_123[k]
                   + pb_x[k] * gg_62[k];

        t_124[k] = f_2 * gh_s_124[k]
                   + pb_x[k] * gg_63[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pb_y, pb_z, fg_25, gf_s_37, gh_s_125, gh_s_126, \
                         gh_s_127, gf_34, gg_59, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * fg_25[k]
                   - f_1 * gf_s_37[k]
                   + f_2 * gh_s_125[k]
                   + f_0 * gf_34[k]
                   + pb_y[k] * gg_59[k];

        t_126[k] = f_2 * gh_s_126[k]
                   + pb_z[k] * gg_59[k];

        t_127[k] = -f_3 * gf_s_37[k]
                   + f_2 * gh_s_127[k]
                   + f_4 * gf_34[k]
                   + pb_z[k] * gg_60[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_y, pb_z, fg_28, gf_s_38, gf_s_40, gh_s_128, \
                         gh_s_129, gh_s_130, gf_35, gf_37, gg_61, \
                         gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = -f_5 * gf_s_38[k]
                   + f_2 * gh_s_128[k]
                   + f_6 * gf_35[k]
                   + pb_z[k] * gg_61[k];

        t_129[k] = f_0 * fg_28[k]
                   + f_2 * gh_s_129[k]
                   + pb_y[k] * gg_63[k];

        t_130[k] = -f_1 * gf_s_40[k]
                   + f_2 * gh_s_130[k]
                   + f_0 * gf_37[k]
                   + pb_z[k] * gg_63[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pb_x, gf_s_41, gf_s_43, gh_s_131, gh_s_132, \
                         gh_s_133, gf_38, gf_40, gg_64, gg_65, gg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -f_5 * gf_s_41[k]
                   + f_2 * gh_s_131[k]
                   + f_6 * gf_38[k]
                   + pb_x[k] * gg_64[k];

        t_132[k] = -f_3 * gf_s_43[k]
                   + f_2 * gh_s_132[k]
                   + f_4 * gf_40[k]
                   + pb_x[k] * gg_65[k];

        t_133[k] = f_2 * gh_s_133[k]
                   + pb_x[k] * gg_67[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_z, pb_x, pb_z, fg_25, fh_41, gh_s_134, \
                         gh_s_135, gh_s_136, gg_66, gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_2 * gh_s_134[k]
                   + pb_x[k] * gg_68[k];

        t_135[k] = pa_z[k] * fh_41[k]
                   + f_2 * gh_s_135[k];

        t_136[k] = f_4 * fg_25[k]
                   + f_2 * gh_s_136[k]
                   + pb_z[k] * gg_66[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_z, pb_y, fg_26, fg_27, fg_32, fh_43, fh_44, \
                         gh_s_137, gh_s_138, gh_s_139, gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * fg_26[k]
                   + pa_z[k] * fh_43[k]
                   + f_2 * gh_s_137[k];

        t_138[k] = f_7 * fg_27[k]
                   + pa_z[k] * fh_44[k]
                   + f_2 * gh_s_138[k];

        t_139[k] = f_7 * fg_32[k]
                   + f_2 * gh_s_139[k]
                   + pb_y[k] * gg_68[k];
    }

#pragma omp simd aligned(t_140, t_141, pa_y, pb_x, dh_s_17, dh_17, fh_54, gf_s_44, gh_s_140, \
                         gh_s_141, gf_41, gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -f_8 * dh_s_17[k]
                   + f_6 * dh_17[k]
                   + pa_y[k] * fh_54[k]
                   + f_2 * gh_s_140[k];

        t_141[k] = -f_1 * gf_s_44[k]
                   + f_2 * gh_s_141[k]
                   + f_0 * gf_41[k]
                   + pb_x[k] * gg_69[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_x, gf_s_45, gf_s_46, gf_s_47, gh_s_142, \
                         gh_s_143, gh_s_144, gf_42, gf_43, gf_44, gg_70, gg_71, \
                         gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -f_5 * gf_s_45[k]
                   + f_2 * gh_s_142[k]
                   + f_6 * gf_42[k]
                   + pb_x[k] * gg_70[k];

        t_143[k] = -f_5 * gf_s_46[k]
                   + f_2 * gh_s_143[k]
                   + f_6 * gf_43[k]
                   + pb_x[k] * gg_71[k];

        t_144[k] = -f_3 * gf_s_47[k]
                   + f_2 * gh_s_144[k]
                   + f_4 * gf_44[k]
                   + pb_x[k] * gg_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pb_x, gf_s_49, gh_s_145, gh_s_146, \
                         gh_s_147, gh_s_148, gf_46, gg_73, gg_74, gg_75, \
                         gg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -f_3 * gf_s_49[k]
                   + f_2 * gh_s_145[k]
                   + f_4 * gf_46[k]
                   + pb_x[k] * gg_73[k];

        t_146[k] = f_2 * gh_s_146[k]
                   + pb_x[k] * gg_74[k];

        t_147[k] = f_2 * gh_s_147[k]
                   + pb_x[k] * gg_75[k];

        t_148[k] = f_2 * gh_s_148[k]
                   + pb_x[k] * gg_77[k];
    }

#pragma omp simd aligned(t_149, t_150, pa_z, pb_z, dh_s_11, dh_11, fg_31, fh_49, gh_s_149, \
                         gh_s_150, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -f_10 * dh_s_11[k]
                   + f_4 * dh_11[k]
                   + pa_z[k] * fh_49[k]
                   + f_2 * gh_s_149[k];

        t_150[k] = f_6 * fg_31[k]
                   + f_2 * gh_s_150[k]
                   + pb_z[k] * gg_74[k];
    }

#pragma omp simd aligned(t_151, t_152, pb_y, fg_36, fg_37, gf_s_48, gf_s_49, gh_s_151, \
                         gh_s_152, gf_45, gf_46, gg_75, gg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_6 * fg_36[k]
                   - f_5 * gf_s_48[k]
                   + f_2 * gh_s_151[k]
                   + f_6 * gf_45[k]
                   + pb_y[k] * gg_75[k];

        t_152[k] = f_6 * fg_37[k]
                   - f_3 * gf_s_49[k]
                   + f_2 * gh_s_152[k]
                   + f_4 * gf_46[k]
                   + pb_y[k] * gg_76[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_y, pb_y, dh_s_27, dh_27, fg_38, fg_40, fh_62, \
                         fh_65, gh_s_153, gh_s_154, gh_s_155, gg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * fg_38[k]
                   + f_2 * gh_s_153[k]
                   + pb_y[k] * gg_77[k];

        t_154[k] = -f_10 * dh_s_27[k]
                   + f_4 * dh_27[k]
                   + pa_y[k] * fh_62[k]
                   + f_2 * gh_s_154[k];

        t_155[k] = f_6 * fg_40[k]
                   + pa_y[k] * fh_65[k]
                   + f_2 * gh_s_155[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pa_y, pb_x, fg_41, fg_45, fh_67, fh_73, \
                         gh_s_156, gh_s_157, gh_s_158, gh_s_159, gg_80, \
                         gg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_7 * fg_41[k]
                   + pa_y[k] * fh_67[k]
                   + f_2 * gh_s_156[k];

        t_157[k] = f_2 * gh_s_157[k]
                   + pb_x[k] * gg_80[k];

        t_158[k] = f_2 * gh_s_158[k]
                   + pb_x[k] * gg_81[k];

        t_159[k] = f_11 * fg_45[k]
                   + pa_y[k] * fh_73[k]
                   + f_2 * gh_s_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_y, pb_z, fg_35, fg_47, fg_48, fh_75, fh_76, \
                         gh_s_160, gh_s_161, gh_s_162, gg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * fg_35[k]
                   + f_2 * gh_s_160[k]
                   + pb_z[k] * gg_80[k];

        t_161[k] = f_7 * fg_47[k]
                   + pa_y[k] * fh_75[k]
                   + f_2 * gh_s_161[k];

        t_162[k] = f_6 * fg_48[k]
                   + pa_y[k] * fh_76[k]
                   + f_2 * gh_s_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_y, pb_x, pb_y, fg_49, fh_78, gf_s_54, \
                         gh_s_163, gh_s_164, gh_s_165, gf_50, gg_83, \
                         gg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_4 * fg_49[k]
                   + f_2 * gh_s_163[k]
                   + pb_y[k] * gg_83[k];

        t_164[k] = pa_y[k] * fh_78[k]
                   + f_2 * gh_s_164[k];

        t_165[k] = -f_1 * gf_s_54[k]
                   + f_2 * gh_s_165[k]
                   + f_0 * gf_50[k]
                   + pb_x[k] * gg_84[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pb_x, pb_y, gf_s_55, gf_s_56, gh_s_166, \
                         gh_s_167, gh_s_168, gf_51, gf_52, gg_84, gg_85, \
                         gg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_2 * gh_s_166[k]
                   + pb_y[k] * gg_84[k];

        t_167[k] = -f_9 * gf_s_55[k]
                   + f_2 * gh_s_167[k]
                   + f_7 * gf_51[k]
                   + pb_x[k] * gg_85[k];

        t_168[k] = -f_5 * gf_s_56[k]
                   + f_2 * gh_s_168[k]
                   + f_6 * gf_52[k]
                   + pb_x[k] * gg_86[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pb_x, pb_y, gf_s_57, gf_s_58, gh_s_169, \
                         gh_s_170, gh_s_171, gf_53, gf_54, gg_85, gg_87, \
                         gg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_2 * gh_s_169[k]
                   + pb_y[k] * gg_85[k];

        t_170[k] = -f_5 * gf_s_57[k]
                   + f_2 * gh_s_170[k]
                   + f_6 * gf_53[k]
                   + pb_x[k] * gg_87[k];

        t_171[k] = -f_3 * gf_s_58[k]
                   + f_2 * gh_s_171[k]
                   + f_4 * gf_54[k]
                   + pb_x[k] * gg_88[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_y, gf_s_59, gf_s_61, gh_s_172, \
                         gh_s_173, gh_s_174, gf_55, gf_57, gg_87, gg_89, \
                         gg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -f_3 * gf_s_59[k]
                   + f_2 * gh_s_172[k]
                   + f_4 * gf_55[k]
                   + pb_x[k] * gg_89[k];

        t_173[k] = f_2 * gh_s_173[k]
                   + pb_y[k] * gg_87[k];

        t_174[k] = -f_3 * gf_s_61[k]
                   + f_2 * gh_s_174[k]
                   + f_4 * gf_57[k]
                   + pb_x[k] * gg_90[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pb_x, gh_s_175, gh_s_176, gh_s_177, \
                         gh_s_178, gg_91, gg_92, gg_93, gg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_2 * gh_s_175[k]
                   + pb_x[k] * gg_91[k];

        t_176[k] = f_2 * gh_s_176[k]
                   + pb_x[k] * gg_92[k];

        t_177[k] = f_2 * gh_s_177[k]
                   + pb_x[k] * gg_93[k];

        t_178[k] = f_2 * gh_s_178[k]
                   + pb_x[k] * gg_95[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pb_y, gf_s_58, gf_s_59, gf_s_60, gh_s_179, \
                         gh_s_180, gh_s_181, gf_54, gf_55, gf_56, gg_91, gg_92, \
                         gg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -f_1 * gf_s_58[k]
                   + f_2 * gh_s_179[k]
                   + f_0 * gf_54[k]
                   + pb_y[k] * gg_91[k];

        t_180[k] = -f_9 * gf_s_59[k]
                   + f_2 * gh_s_180[k]
                   + f_7 * gf_55[k]
                   + pb_y[k] * gg_92[k];

        t_181[k] = -f_5 * gf_s_60[k]
                   + f_2 * gh_s_181[k]
                   + f_6 * gf_56[k]
                   + pb_y[k] * gg_93[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pb_y, pb_z, fg_49, gf_s_61, gh_s_182, gh_s_183, \
                         gh_s_184, gf_57, gg_94, gg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -f_3 * gf_s_61[k]
                   + f_2 * gh_s_182[k]
                   + f_4 * gf_57[k]
                   + pb_y[k] * gg_94[k];

        t_183[k] = f_2 * gh_s_183[k]
                   + pb_y[k] * gg_95[k];

        t_184[k] = f_0 * fg_49[k]
                   - f_1 * gf_s_61[k]
                   + f_2 * gh_s_184[k]
                   + f_0 * gf_57[k]
                   + pb_z[k] * gg_95[k];
    }
}

auto
compute_prim_gh_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = beta / p;
    const auto f_10 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_14 = buffer.data(dh_s + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_14 = buffer.data(dh + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);

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

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_11 = buffer.data(gf_s + 11);
    const auto *gf_s_12 = buffer.data(gf_s + 12);
    const auto *gf_s_18 = buffer.data(gf_s + 18);
    const auto *gf_s_19 = buffer.data(gf_s + 19);
    const auto *gf_s_29 = buffer.data(gf_s + 29);
    const auto *gf_s_30 = buffer.data(gf_s + 30);
    const auto *gf_s_32 = buffer.data(gf_s + 32);
    const auto *gf_s_33 = buffer.data(gf_s + 33);
    const auto *gf_s_42 = buffer.data(gf_s + 42);
    const auto *gf_s_43 = buffer.data(gf_s + 43);
    const auto *gf_s_48 = buffer.data(gf_s + 48);
    const auto *gf_s_50 = buffer.data(gf_s + 50);
    const auto *gf_s_51 = buffer.data(gf_s + 51);
    const auto *gf_s_52 = buffer.data(gf_s + 52);
    const auto *gf_s_54 = buffer.data(gf_s + 54);
    const auto *gf_s_55 = buffer.data(gf_s + 55);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_20 = buffer.data(gh_s + 20);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_31 = buffer.data(gh_s + 31);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_33 = buffer.data(gh_s + 33);
    const auto *gh_s_34 = buffer.data(gh_s + 34);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gf_0, gg_0, gg_1, gg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_1[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];

        t_2[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_2[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, fh_0, gf_s_1, gf_s_2, gh_s_3, \
                         gh_s_4, gh_s_5, gf_1, gf_2, gg_3, gg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_3[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_4[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_4[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_4[k];

        t_5[k] = pa_y[k] * fh_0[k]
                 + f_2 * gh_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, dh_s_3, dh_3, fg_1, fh_0, fh_1, fh_4, \
                         gh_s_6, gh_s_7, gh_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_7 * dh_s_3[k]
                 + f_6 * dh_3[k]
                 + pa_x[k] * fh_4[k]
                 + f_2 * gh_s_6[k];

        t_7[k] = pa_z[k] * fh_0[k]
                 + f_2 * gh_s_7[k];

        t_8[k] = f_6 * fg_1[k]
                 + pa_z[k] * fh_1[k]
                 + f_2 * gh_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pa_z, dh_s_4, dh_4, fg_3, fh_2, fh_8, gh_s_9, \
                         gh_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * fg_3[k]
                 + pa_z[k] * fh_2[k]
                 + f_2 * gh_s_9[k];

        t_10[k] = -f_7 * dh_s_4[k]
                  + f_6 * dh_4[k]
                  + pa_x[k] * fh_8[k]
                  + f_2 * gh_s_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_y, pb_x, dh_s_0, dh_0, fg_9, fh_3, gf_s_11, gh_s_11, \
                         gh_s_12, gf_11, gg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_9 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_3[k]
                  + f_2 * gh_s_11[k];

        t_12[k] = f_6 * fg_9[k]
                  - f_5 * gf_s_11[k]
                  + f_2 * gh_s_12[k]
                  + f_6 * gf_11[k]
                  + pb_x[k] * gg_15[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pb_x, dh_s_5, dh_5, fg_10, fh_9, gf_s_12, gh_s_13, \
                         gh_s_14, gf_12, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * fg_10[k]
                  - f_3 * gf_s_12[k]
                  + f_2 * gh_s_13[k]
                  + f_4 * gf_12[k]
                  + pb_x[k] * gg_16[k];

        t_14[k] = -f_9 * dh_s_5[k]
                  + f_4 * dh_5[k]
                  + pa_x[k] * fh_9[k]
                  + f_2 * gh_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, dh_s_8, dh_8, fh_6, fh_7, fh_10, \
                         gh_s_15, gh_s_16, gh_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * fh_6[k]
                  + f_2 * gh_s_15[k];

        t_16[k] = pa_y[k] * fh_7[k]
                  + f_2 * gh_s_16[k];

        t_17[k] = -f_9 * dh_s_8[k]
                  + f_4 * dh_8[k]
                  + pa_x[k] * fh_10[k]
                  + f_2 * gh_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, pa_x, pa_z, dh_s_0, dh_s_9, dh_0, dh_9, fh_5, fh_11, \
                         gh_s_18, gh_s_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_9 * dh_s_9[k]
                  + f_4 * dh_9[k]
                  + pa_x[k] * fh_11[k]
                  + f_2 * gh_s_18[k];

        t_19[k] = -f_9 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_z[k] * fh_5[k]
                  + f_2 * gh_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, pb_x, fg_13, fg_14, gf_s_18, gf_s_19, gh_s_20, gh_s_21, \
                         gf_18, gf_19, gg_23, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * fg_13[k]
                  - f_5 * gf_s_18[k]
                  + f_2 * gh_s_20[k]
                  + f_6 * gf_18[k]
                  + pb_x[k] * gg_23[k];

        t_21[k] = f_6 * fg_14[k]
                  - f_3 * gf_s_19[k]
                  + f_2 * gh_s_21[k]
                  + f_4 * gf_19[k]
                  + pb_x[k] * gg_24[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, dh_s_14, dh_14, fh_12, fh_13, fh_17, \
                         fh_18, gh_s_22, gh_s_23, gh_s_24, gh_s_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_9 * dh_s_14[k]
                  + f_4 * dh_14[k]
                  + pa_x[k] * fh_12[k]
                  + f_2 * gh_s_22[k];

        t_23[k] = pa_x[k] * fh_13[k]
                  + f_2 * gh_s_23[k];

        t_24[k] = pa_x[k] * fh_17[k]
                  + f_2 * gh_s_24[k];

        t_25[k] = pa_x[k] * fh_18[k]
                  + f_2 * gh_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, fh_19, fh_20, fh_21, fh_22, \
                         fh_27, gh_s_26, gh_s_27, gh_s_28, gh_s_29, \
                         gh_s_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * fh_19[k]
                  + f_2 * gh_s_26[k];

        t_27[k] = pa_x[k] * fh_20[k]
                  + f_2 * gh_s_27[k];

        t_28[k] = pa_x[k] * fh_21[k]
                  + f_2 * gh_s_28[k];

        t_29[k] = pa_x[k] * fh_22[k]
                  + f_2 * gh_s_29[k];

        t_30[k] = pa_x[k] * fh_27[k]
                  + f_2 * gh_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, gf_s_29, gf_s_30, gf_s_32, gh_s_31, gh_s_32, \
                         gh_s_33, gf_29, gf_30, gf_32, gg_39, gg_40, \
                         gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_1 * gf_s_29[k]
                  + f_2 * gh_s_31[k]
                  + f_0 * gf_29[k]
                  + pb_x[k] * gg_39[k];

        t_32[k] = -f_5 * gf_s_30[k]
                  + f_2 * gh_s_32[k]
                  + f_6 * gf_30[k]
                  + pb_x[k] * gg_40[k];

        t_33[k] = -f_3 * gf_s_32[k]
                  + f_2 * gh_s_33[k]
                  + f_4 * gf_32[k]
                  + pb_x[k] * gg_42[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_y, pb_z, fg_18, gf_s_32, gf_s_33, gh_s_34, \
                         gh_s_35, gh_s_36, gf_32, gf_33, gg_43, gg_44, \
                         gg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * fg_18[k]
                  - f_1 * gf_s_32[k]
                  + f_2 * gh_s_34[k]
                  + f_0 * gf_32[k]
                  + pb_y[k] * gg_43[k];

        t_35[k] = -f_3 * gf_s_32[k]
                  + f_2 * gh_s_35[k]
                  + f_4 * gf_32[k]
                  + pb_z[k] * gg_44[k];

        t_36[k] = -f_5 * gf_s_33[k]
                  + f_2 * gh_s_36[k]
                  + f_6 * gf_33[k]
                  + pb_z[k] * gg_45[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_z, fg_19, fg_20, fh_13, fh_14, fh_15, gh_s_37, \
                         gh_s_38, gh_s_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_z[k] * fh_13[k]
                  + f_2 * gh_s_37[k];

        t_38[k] = f_6 * fg_19[k]
                  + pa_z[k] * fh_14[k]
                  + f_2 * gh_s_38[k];

        t_39[k] = f_8 * fg_20[k]
                  + pa_z[k] * fh_15[k]
                  + f_2 * gh_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, pa_y, pa_z, dh_s_5, dh_s_10, dh_5, dh_10, fh_16, fh_19, \
                         gh_s_40, gh_s_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_7 * dh_s_10[k]
                  + f_6 * dh_10[k]
                  + pa_y[k] * fh_19[k]
                  + f_2 * gh_s_40[k];

        t_41[k] = -f_9 * dh_s_5[k]
                  + f_4 * dh_5[k]
                  + pa_z[k] * fh_16[k]
                  + f_2 * gh_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, pb_y, fg_25, fg_26, gf_s_42, gf_s_43, gh_s_42, gh_s_43, \
                         gf_42, gf_43, gg_59, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_6 * fg_25[k]
                  - f_5 * gf_s_42[k]
                  + f_2 * gh_s_42[k]
                  + f_6 * gf_42[k]
                  + pb_y[k] * gg_59[k];

        t_43[k] = f_6 * fg_26[k]
                  - f_3 * gf_s_43[k]
                  + f_2 * gh_s_43[k]
                  + f_4 * gf_43[k]
                  + pb_y[k] * gg_60[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, dh_s_14, dh_14, fg_31, fg_32, fh_23, fh_24, \
                         fh_25, gh_s_44, gh_s_45, gh_s_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = -f_9 * dh_s_14[k]
                  + f_4 * dh_14[k]
                  + pa_y[k] * fh_23[k]
                  + f_2 * gh_s_44[k];

        t_45[k] = f_10 * fg_31[k]
                  + pa_y[k] * fh_24[k]
                  + f_2 * gh_s_45[k];

        t_46[k] = f_8 * fg_32[k]
                  + pa_y[k] * fh_25[k]
                  + f_2 * gh_s_46[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pb_x, fg_33, fh_26, fh_27, gf_s_48, gh_s_47, \
                         gh_s_48, gh_s_49, gf_48, gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * fg_33[k]
                  + pa_y[k] * fh_26[k]
                  + f_2 * gh_s_47[k];

        t_48[k] = pa_y[k] * fh_27[k]
                  + f_2 * gh_s_48[k];

        t_49[k] = -f_1 * gf_s_48[k]
                  + f_2 * gh_s_49[k]
                  + f_0 * gf_48[k]
                  + pb_x[k] * gg_68[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, gf_s_50, gf_s_51, gf_s_52, gh_s_50, gh_s_51, \
                         gh_s_52, gf_50, gf_51, gf_52, gg_70, gg_71, \
                         gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_5 * gf_s_50[k]
                  + f_2 * gh_s_50[k]
                  + f_6 * gf_50[k]
                  + pb_x[k] * gg_70[k];

        t_51[k] = -f_5 * gf_s_51[k]
                  + f_2 * gh_s_51[k]
                  + f_6 * gf_51[k]
                  + pb_x[k] * gg_71[k];

        t_52[k] = -f_3 * gf_s_52[k]
                  + f_2 * gh_s_52[k]
                  + f_4 * gf_52[k]
                  + pb_x[k] * gg_72[k];
    }

#pragma omp simd aligned(t_53, t_54, pb_x, pb_y, gf_s_52, gf_s_55, gh_s_53, gh_s_54, gf_52, \
                         gf_55, gg_73, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_3 * gf_s_55[k]
                  + f_2 * gh_s_53[k]
                  + f_4 * gf_55[k]
                  + pb_x[k] * gg_73[k];

        t_54[k] = -f_1 * gf_s_52[k]
                  + f_2 * gh_s_54[k]
                  + f_0 * gf_52[k]
                  + pb_y[k] * gg_74[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_y, pb_z, fg_34, gf_s_54, gf_s_55, gh_s_55, \
                         gh_s_56, gh_s_57, gf_54, gf_55, gg_76, gg_77, \
                         gg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_5 * gf_s_54[k]
                  + f_2 * gh_s_55[k]
                  + f_6 * gf_54[k]
                  + pb_y[k] * gg_76[k];

        t_56[k] = -f_3 * gf_s_55[k]
                  + f_2 * gh_s_56[k]
                  + f_4 * gf_55[k]
                  + pb_y[k] * gg_77[k];

        t_57[k] = f_0 * fg_34[k]
                  - f_1 * gf_s_55[k]
                  + f_2 * gh_s_57[k]
                  + f_0 * gf_55[k]
                  + pb_z[k] * gg_78[k];
    }
}

auto
compute_prim_gh_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = beta / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_21 = buffer.data(dh_s + 21);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_21 = buffer.data(dh + 21);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_8 = buffer.data(fg + 8);
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
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);

    const auto *fh_0 = buffer.data(fh + 0);
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
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
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
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_47 = buffer.data(fh + 47);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_3 = buffer.data(gf_s + 3);
    const auto *gf_s_4 = buffer.data(gf_s + 4);
    const auto *gf_s_10 = buffer.data(gf_s + 10);
    const auto *gf_s_11 = buffer.data(gf_s + 11);
    const auto *gf_s_13 = buffer.data(gf_s + 13);
    const auto *gf_s_14 = buffer.data(gf_s + 14);
    const auto *gf_s_15 = buffer.data(gf_s + 15);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_21 = buffer.data(gf_s + 21);
    const auto *gf_s_22 = buffer.data(gf_s + 22);
    const auto *gf_s_23 = buffer.data(gf_s + 23);
    const auto *gf_s_24 = buffer.data(gf_s + 24);
    const auto *gf_s_25 = buffer.data(gf_s + 25);
    const auto *gf_s_26 = buffer.data(gf_s + 26);
    const auto *gf_s_27 = buffer.data(gf_s + 27);
    const auto *gf_s_29 = buffer.data(gf_s + 29);
    const auto *gf_s_30 = buffer.data(gf_s + 30);
    const auto *gf_s_31 = buffer.data(gf_s + 31);
    const auto *gf_s_32 = buffer.data(gf_s + 32);
    const auto *gf_s_33 = buffer.data(gf_s + 33);
    const auto *gf_s_34 = buffer.data(gf_s + 34);
    const auto *gf_s_35 = buffer.data(gf_s + 35);
    const auto *gf_s_40 = buffer.data(gf_s + 40);
    const auto *gf_s_42 = buffer.data(gf_s + 42);
    const auto *gf_s_43 = buffer.data(gf_s + 43);
    const auto *gf_s_44 = buffer.data(gf_s + 44);
    const auto *gf_s_45 = buffer.data(gf_s + 45);
    const auto *gf_s_46 = buffer.data(gf_s + 46);
    const auto *gf_s_47 = buffer.data(gf_s + 47);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_20 = buffer.data(gh_s + 20);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_31 = buffer.data(gh_s + 31);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_33 = buffer.data(gh_s + 33);
    const auto *gh_s_34 = buffer.data(gh_s + 34);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);
    const auto *gh_s_58 = buffer.data(gh_s + 58);
    const auto *gh_s_59 = buffer.data(gh_s + 59);
    const auto *gh_s_60 = buffer.data(gh_s + 60);
    const auto *gh_s_61 = buffer.data(gh_s + 61);
    const auto *gh_s_62 = buffer.data(gh_s + 62);
    const auto *gh_s_63 = buffer.data(gh_s + 63);
    const auto *gh_s_64 = buffer.data(gh_s + 64);
    const auto *gh_s_65 = buffer.data(gh_s + 65);
    const auto *gh_s_66 = buffer.data(gh_s + 66);
    const auto *gh_s_67 = buffer.data(gh_s + 67);
    const auto *gh_s_68 = buffer.data(gh_s + 68);
    const auto *gh_s_69 = buffer.data(gh_s + 69);
    const auto *gh_s_70 = buffer.data(gh_s + 70);
    const auto *gh_s_71 = buffer.data(gh_s + 71);
    const auto *gh_s_72 = buffer.data(gh_s + 72);
    const auto *gh_s_73 = buffer.data(gh_s + 73);
    const auto *gh_s_74 = buffer.data(gh_s + 74);
    const auto *gh_s_75 = buffer.data(gh_s + 75);
    const auto *gh_s_76 = buffer.data(gh_s + 76);
    const auto *gh_s_77 = buffer.data(gh_s + 77);
    const auto *gh_s_78 = buffer.data(gh_s + 78);
    const auto *gh_s_79 = buffer.data(gh_s + 79);
    const auto *gh_s_80 = buffer.data(gh_s + 80);
    const auto *gh_s_81 = buffer.data(gh_s + 81);
    const auto *gh_s_82 = buffer.data(gh_s + 82);
    const auto *gh_s_83 = buffer.data(gh_s + 83);
    const auto *gh_s_84 = buffer.data(gh_s + 84);
    const auto *gh_s_85 = buffer.data(gh_s + 85);
    const auto *gh_s_86 = buffer.data(gh_s + 86);
    const auto *gh_s_87 = buffer.data(gh_s + 87);
    const auto *gh_s_88 = buffer.data(gh_s + 88);
    const auto *gh_s_89 = buffer.data(gh_s + 89);
    const auto *gh_s_90 = buffer.data(gh_s + 90);
    const auto *gh_s_91 = buffer.data(gh_s + 91);
    const auto *gh_s_92 = buffer.data(gh_s + 92);
    const auto *gh_s_93 = buffer.data(gh_s + 93);
    const auto *gh_s_94 = buffer.data(gh_s + 94);
    const auto *gh_s_95 = buffer.data(gh_s + 95);
    const auto *gh_s_96 = buffer.data(gh_s + 96);
    const auto *gh_s_97 = buffer.data(gh_s + 97);
    const auto *gh_s_98 = buffer.data(gh_s + 98);
    const auto *gh_s_99 = buffer.data(gh_s + 99);
    const auto *gh_s_100 = buffer.data(gh_s + 100);
    const auto *gh_s_101 = buffer.data(gh_s + 101);
    const auto *gh_s_102 = buffer.data(gh_s + 102);
    const auto *gh_s_103 = buffer.data(gh_s + 103);
    const auto *gh_s_104 = buffer.data(gh_s + 104);
    const auto *gh_s_105 = buffer.data(gh_s + 105);
    const auto *gh_s_106 = buffer.data(gh_s + 106);
    const auto *gh_s_107 = buffer.data(gh_s + 107);
    const auto *gh_s_108 = buffer.data(gh_s + 108);
    const auto *gh_s_109 = buffer.data(gh_s + 109);
    const auto *gh_s_110 = buffer.data(gh_s + 110);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_32 = buffer.data(gg + 32);
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
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gh_s_3, gf_0, gg_0, gg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_2 * gh_s_1[k]
                 + pb_y[k] * gg_0[k];

        t_2[k] = f_2 * gh_s_2[k]
                 + pb_z[k] * gg_0[k];

        t_3[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_3[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, gf_s_0, gf_s_1, gh_s_4, gh_s_5, gh_s_6, \
                         gf_0, gf_1, gg_2, gg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_4[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];

        t_5[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_5[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_6[k] = f_2 * gh_s_6[k]
                 + pb_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, gf_s_2, gf_s_3, gh_s_7, gh_s_8, gh_s_9, \
                         gf_2, gf_3, gg_4, gg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * gh_s_7[k]
                 + pb_y[k] * gg_4[k];

        t_8[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_8[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_4[k];

        t_9[k] = -f_1 * gf_s_3[k]
                 + f_2 * gh_s_9[k]
                 + f_0 * gf_3[k]
                 + pb_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_z, fg_1, fh_0, fh_2, gf_s_4, gh_s_10, \
                         gh_s_11, gh_s_12, gf_4, gg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * gf_s_4[k]
                  + f_2 * gh_s_10[k]
                  + f_0 * gf_4[k]
                  + pb_z[k] * gg_6[k];

        t_11[k] = pa_y[k] * fh_0[k]
                  + f_2 * gh_s_11[k];

        t_12[k] = f_6 * fg_1[k]
                  + pa_y[k] * fh_2[k]
                  + f_2 * gh_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pa_z, dh_s_4, dh_4, fg_3, fh_0, fh_4, \
                         fh_7, gh_s_13, gh_s_14, gh_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * fg_3[k]
                  + pa_y[k] * fh_4[k]
                  + f_2 * gh_s_13[k];

        t_14[k] = -f_8 * dh_s_4[k]
                  + f_6 * dh_4[k]
                  + pa_x[k] * fh_7[k]
                  + f_2 * gh_s_14[k];

        t_15[k] = pa_z[k] * fh_0[k]
                  + f_2 * gh_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, fg_0, fg_2, fg_4, fh_3, fh_5, gh_s_16, \
                         gh_s_17, gh_s_18, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_16[k]
                  + pb_z[k] * gg_10[k];

        t_17[k] = f_6 * fg_2[k]
                  + pa_z[k] * fh_3[k]
                  + f_2 * gh_s_17[k];

        t_18[k] = f_7 * fg_4[k]
                  + pa_z[k] * fh_5[k]
                  + f_2 * gh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_x, pa_y, dh_s_0, dh_s_5, dh_0, dh_5, fh_6, fh_11, \
                         gh_s_19, gh_s_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_8 * dh_s_5[k]
                  + f_6 * dh_5[k]
                  + pa_x[k] * fh_11[k]
                  + f_2 * gh_s_19[k];

        t_20[k] = -f_9 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_6[k]
                  + f_2 * gh_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, fg_11, fg_12, gf_s_10, gf_s_11, gh_s_21, gh_s_22, \
                         gf_10, gf_11, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_6 * fg_11[k]
                  - f_5 * gf_s_10[k]
                  + f_2 * gh_s_21[k]
                  + f_6 * gf_10[k]
                  + pb_x[k] * gg_13[k];

        t_22[k] = f_6 * fg_12[k]
                  - f_3 * gf_s_11[k]
                  + f_2 * gh_s_22[k]
                  + f_4 * gf_11[k]
                  + pb_x[k] * gg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, dh_s_8, dh_8, fg_13, fh_9, fh_15, \
                         gh_s_23, gh_s_24, gh_s_25, gg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * fg_13[k]
                  + f_2 * gh_s_23[k]
                  + pb_x[k] * gg_15[k];

        t_24[k] = -f_9 * dh_s_8[k]
                  + f_4 * dh_8[k]
                  + pa_x[k] * fh_15[k]
                  + f_2 * gh_s_24[k];

        t_25[k] = pa_y[k] * fh_9[k]
                  + f_2 * gh_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, pa_y, dh_s_12, dh_s_13, dh_12, dh_13, fh_10, \
                         fh_16, fh_17, gh_s_26, gh_s_27, gh_s_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * fh_10[k]
                  + f_2 * gh_s_26[k];

        t_27[k] = -f_9 * dh_s_12[k]
                  + f_4 * dh_12[k]
                  + pa_x[k] * fh_16[k]
                  + f_2 * gh_s_27[k];

        t_28[k] = -f_9 * dh_s_13[k]
                  + f_4 * dh_13[k]
                  + pa_x[k] * fh_17[k]
                  + f_2 * gh_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_z, pb_z, dh_s_0, dh_0, fg_8, fh_8, gh_s_29, gh_s_30, \
                         gg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_9 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_z[k] * fh_8[k]
                  + f_2 * gh_s_29[k];

        t_30[k] = f_6 * fg_8[k]
                  + f_2 * gh_s_30[k]
                  + pb_z[k] * gg_19[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pb_y, fg_15, gf_s_13, gf_s_15, gh_s_31, gh_s_32, \
                         gf_13, gf_15, gg_20, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_3 * gf_s_13[k]
                  + f_2 * gh_s_31[k]
                  + f_4 * gf_13[k]
                  + pb_y[k] * gg_20[k];

        t_32[k] = f_6 * fg_15[k]
                  - f_5 * gf_s_15[k]
                  + f_2 * gh_s_32[k]
                  + f_6 * gf_15[k]
                  + pb_x[k] * gg_22[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, fg_16, gf_s_14, gf_s_16, gh_s_33, gh_s_34, \
                         gf_14, gf_16, gg_21, gg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_5 * gf_s_14[k]
                  + f_2 * gh_s_33[k]
                  + f_6 * gf_14[k]
                  + pb_y[k] * gg_21[k];

        t_34[k] = f_6 * fg_16[k]
                  - f_3 * gf_s_16[k]
                  + f_2 * gh_s_34[k]
                  + f_4 * gf_16[k]
                  + pb_x[k] * gg_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_x, dh_s_21, dh_21, fg_17, fg_18, fh_21, \
                         fh_22, gh_s_35, gh_s_36, gh_s_37, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * fg_17[k]
                  + f_2 * gh_s_35[k]
                  + pb_x[k] * gg_24[k];

        t_36[k] = -f_9 * dh_s_21[k]
                  + f_4 * dh_21[k]
                  + pa_x[k] * fh_21[k]
                  + f_2 * gh_s_36[k];

        t_37[k] = f_10 * fg_18[k]
                  + pa_x[k] * fh_22[k]
                  + f_2 * gh_s_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_x, fg_19, fg_20, fg_21, fh_23, fh_24, \
                         gh_s_38, gh_s_39, gh_s_40, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_7 * fg_19[k]
                  + pa_x[k] * fh_23[k]
                  + f_2 * gh_s_38[k];

        t_39[k] = f_6 * fg_20[k]
                  + pa_x[k] * fh_24[k]
                  + f_2 * gh_s_39[k];

        t_40[k] = f_4 * fg_21[k]
                  + f_2 * gh_s_40[k]
                  + pb_x[k] * gg_28[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_x, fh_25, fh_30, fh_31, fh_32, \
                         fh_33, gh_s_41, gh_s_42, gh_s_43, gh_s_44, \
                         gh_s_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_x[k] * fh_25[k]
                  + f_2 * gh_s_41[k];

        t_42[k] = pa_x[k] * fh_30[k]
                  + f_2 * gh_s_42[k];

        t_43[k] = pa_x[k] * fh_31[k]
                  + f_2 * gh_s_43[k];

        t_44[k] = pa_x[k] * fh_32[k]
                  + f_2 * gh_s_44[k];

        t_45[k] = pa_x[k] * fh_33[k]
                  + f_2 * gh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_z, fg_14, fg_31, fh_34, fh_35, \
                         fh_38, gh_s_46, gh_s_47, gh_s_48, gh_s_49, \
                         gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * fh_34[k]
                  + f_2 * gh_s_46[k];

        t_47[k] = pa_x[k] * fh_35[k]
                  + f_2 * gh_s_47[k];

        t_48[k] = f_10 * fg_31[k]
                  + pa_x[k] * fh_38[k]
                  + f_2 * gh_s_48[k];

        t_49[k] = f_7 * fg_14[k]
                  + f_2 * gh_s_49[k]
                  + pb_z[k] * gg_32[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, fg_34, fg_36, fg_40, fh_40, fh_42, \
                         gh_s_50, gh_s_51, gh_s_52, gg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * fg_34[k]
                  + pa_x[k] * fh_40[k]
                  + f_2 * gh_s_50[k];

        t_51[k] = f_6 * fg_36[k]
                  + pa_x[k] * fh_42[k]
                  + f_2 * gh_s_51[k];

        t_52[k] = f_4 * fg_40[k]
                  + f_2 * gh_s_52[k]
                  + pb_x[k] * gg_35[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pb_x, fh_47, gf_s_21, gf_s_22, gh_s_53, \
                         gh_s_54, gh_s_55, gf_21, gf_22, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * fh_47[k]
                  + f_2 * gh_s_53[k];

        t_54[k] = -f_1 * gf_s_21[k]
                  + f_2 * gh_s_54[k]
                  + f_0 * gf_21[k]
                  + pb_x[k] * gg_36[k];

        t_55[k] = -f_5 * gf_s_22[k]
                  + f_2 * gh_s_55[k]
                  + f_6 * gf_22[k]
                  + pb_x[k] * gg_37[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_x, pb_z, gf_s_23, gf_s_24, gh_s_56, gh_s_57, \
                         gh_s_58, gf_23, gf_24, gg_37, gg_38, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -f_5 * gf_s_23[k]
                  + f_2 * gh_s_56[k]
                  + f_6 * gf_23[k]
                  + pb_x[k] * gg_38[k];

        t_57[k] = -f_3 * gf_s_24[k]
                  + f_2 * gh_s_57[k]
                  + f_4 * gf_24[k]
                  + pb_x[k] * gg_39[k];

        t_58[k] = f_2 * gh_s_58[k]
                  + pb_z[k] * gg_37[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_x, gf_s_26, gh_s_59, gh_s_60, gh_s_61, gf_26, \
                         gg_40, gg_41, gg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_3 * gf_s_26[k]
                  + f_2 * gh_s_59[k]
                  + f_4 * gf_26[k]
                  + pb_x[k] * gg_40[k];

        t_60[k] = f_2 * gh_s_60[k]
                  + pb_x[k] * gg_41[k];

        t_61[k] = f_2 * gh_s_61[k]
                  + pb_x[k] * gg_43[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pb_y, pb_z, fg_21, gf_s_24, gh_s_62, gh_s_63, \
                         gh_s_64, gf_24, gg_41, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * fg_21[k]
                  - f_1 * gf_s_24[k]
                  + f_2 * gh_s_62[k]
                  + f_0 * gf_24[k]
                  + pb_y[k] * gg_41[k];

        t_63[k] = f_2 * gh_s_63[k]
                  + pb_z[k] * gg_41[k];

        t_64[k] = -f_3 * gf_s_24[k]
                  + f_2 * gh_s_64[k]
                  + f_4 * gf_24[k]
                  + pb_z[k] * gg_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_y, pb_z, fg_24, gf_s_25, gf_s_26, gh_s_65, \
                         gh_s_66, gh_s_67, gf_25, gf_26, gg_43, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -f_5 * gf_s_25[k]
                  + f_2 * gh_s_65[k]
                  + f_6 * gf_25[k]
                  + pb_z[k] * gg_43[k];

        t_66[k] = f_0 * fg_24[k]
                  + f_2 * gh_s_66[k]
                  + pb_y[k] * gg_44[k];

        t_67[k] = -f_1 * gf_s_26[k]
                  + f_2 * gh_s_67[k]
                  + f_0 * gf_26[k]
                  + pb_z[k] * gg_44[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_z, pb_x, fh_25, gf_s_27, gf_s_29, gh_s_68, \
                         gh_s_69, gh_s_70, gf_27, gf_29, gg_45, gg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -f_5 * gf_s_27[k]
                  + f_2 * gh_s_68[k]
                  + f_6 * gf_27[k]
                  + pb_x[k] * gg_45[k];

        t_69[k] = -f_3 * gf_s_29[k]
                  + f_2 * gh_s_69[k]
                  + f_4 * gf_29[k]
                  + pb_x[k] * gg_46[k];

        t_70[k] = pa_z[k] * fh_25[k]
                  + f_2 * gh_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_z, pb_z, fg_21, fg_22, fg_23, fh_27, fh_28, \
                         gh_s_71, gh_s_72, gh_s_73, gg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_4 * fg_21[k]
                  + f_2 * gh_s_71[k]
                  + pb_z[k] * gg_47[k];

        t_72[k] = f_6 * fg_22[k]
                  + pa_z[k] * fh_27[k]
                  + f_2 * gh_s_72[k];

        t_73[k] = f_7 * fg_23[k]
                  + pa_z[k] * fh_28[k]
                  + f_2 * gh_s_73[k];
    }

#pragma omp simd aligned(t_74, t_75, pa_y, pb_y, dh_s_14, dh_14, fg_26, fh_32, gh_s_74, \
                         gh_s_75, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_7 * fg_26[k]
                  + f_2 * gh_s_74[k]
                  + pb_y[k] * gg_48[k];

        t_75[k] = -f_8 * dh_s_14[k]
                  + f_6 * dh_14[k]
                  + pa_y[k] * fh_32[k]
                  + f_2 * gh_s_75[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_x, gf_s_30, gf_s_31, gf_s_32, gh_s_76, gh_s_77, \
                         gh_s_78, gf_30, gf_31, gf_32, gg_49, gg_50, \
                         gg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -f_1 * gf_s_30[k]
                  + f_2 * gh_s_76[k]
                  + f_0 * gf_30[k]
                  + pb_x[k] * gg_49[k];

        t_77[k] = -f_5 * gf_s_31[k]
                  + f_2 * gh_s_77[k]
                  + f_6 * gf_31[k]
                  + pb_x[k] * gg_50[k];

        t_78[k] = -f_5 * gf_s_32[k]
                  + f_2 * gh_s_78[k]
                  + f_6 * gf_32[k]
                  + pb_x[k] * gg_51[k];
    }

#pragma omp simd aligned(t_79, t_80, pb_x, gf_s_33, gf_s_35, gh_s_79, gh_s_80, gf_33, gf_35, \
                         gg_52, gg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -f_3 * gf_s_33[k]
                  + f_2 * gh_s_79[k]
                  + f_4 * gf_33[k]
                  + pb_x[k] * gg_52[k];

        t_80[k] = -f_3 * gf_s_35[k]
                  + f_2 * gh_s_80[k]
                  + f_4 * gf_35[k]
                  + pb_x[k] * gg_53[k];
    }

#pragma omp simd aligned(t_81, t_82, pa_z, pb_z, dh_s_8, dh_8, fg_25, fh_29, gh_s_81, gh_s_82, \
                         gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_9 * dh_s_8[k]
                  + f_4 * dh_8[k]
                  + pa_z[k] * fh_29[k]
                  + f_2 * gh_s_81[k];

        t_82[k] = f_6 * fg_25[k]
                  + f_2 * gh_s_82[k]
                  + pb_z[k] * gg_54[k];
    }

#pragma omp simd aligned(t_83, t_84, pb_y, fg_28, fg_29, gf_s_34, gf_s_35, gh_s_83, gh_s_84, \
                         gf_34, gf_35, gg_55, gg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_6 * fg_28[k]
                  - f_5 * gf_s_34[k]
                  + f_2 * gh_s_83[k]
                  + f_6 * gf_34[k]
                  + pb_y[k] * gg_55[k];

        t_84[k] = f_6 * fg_29[k]
                  - f_3 * gf_s_35[k]
                  + f_2 * gh_s_84[k]
                  + f_4 * gf_35[k]
                  + pb_y[k] * gg_56[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pa_y, pb_y, dh_s_21, dh_21, fg_30, fg_32, fh_37, \
                         fh_39, gh_s_85, gh_s_86, gh_s_87, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_6 * fg_30[k]
                  + f_2 * gh_s_85[k]
                  + pb_y[k] * gg_57[k];

        t_86[k] = -f_9 * dh_s_21[k]
                  + f_4 * dh_21[k]
                  + pa_y[k] * fh_37[k]
                  + f_2 * gh_s_86[k];

        t_87[k] = f_6 * fg_32[k]
                  + pa_y[k] * fh_39[k]
                  + f_2 * gh_s_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pa_y, pb_z, fg_27, fg_33, fg_37, fh_41, fh_43, \
                         gh_s_88, gh_s_89, gh_s_90, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_7 * fg_33[k]
                  + pa_y[k] * fh_41[k]
                  + f_2 * gh_s_88[k];

        t_89[k] = f_10 * fg_37[k]
                  + pa_y[k] * fh_43[k]
                  + f_2 * gh_s_89[k];

        t_90[k] = f_7 * fg_27[k]
                  + f_2 * gh_s_90[k]
                  + pb_z[k] * gg_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_y, pb_y, fg_38, fg_39, fg_40, fh_44, fh_45, \
                         gh_s_91, gh_s_92, gh_s_93, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_7 * fg_38[k]
                  + pa_y[k] * fh_44[k]
                  + f_2 * gh_s_91[k];

        t_92[k] = f_6 * fg_39[k]
                  + pa_y[k] * fh_45[k]
                  + f_2 * gh_s_92[k];

        t_93[k] = f_4 * fg_40[k]
                  + f_2 * gh_s_93[k]
                  + pb_y[k] * gg_63[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pb_x, pb_y, fh_47, gf_s_40, gh_s_94, gh_s_95, \
                         gh_s_96, gf_39, gg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pa_y[k] * fh_47[k]
                  + f_2 * gh_s_94[k];

        t_95[k] = -f_1 * gf_s_40[k]
                  + f_2 * gh_s_95[k]
                  + f_0 * gf_39[k]
                  + pb_x[k] * gg_64[k];

        t_96[k] = f_2 * gh_s_96[k]
                  + pb_y[k] * gg_64[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, gf_s_42, gf_s_43, gf_s_44, gh_s_97, gh_s_98, \
                         gh_s_99, gf_41, gf_42, gf_43, gg_66, gg_67, \
                         gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -f_5 * gf_s_42[k]
                  + f_2 * gh_s_97[k]
                  + f_6 * gf_41[k]
                  + pb_x[k] * gg_66[k];

        t_98[k] = -f_5 * gf_s_43[k]
                  + f_2 * gh_s_98[k]
                  + f_6 * gf_42[k]
                  + pb_x[k] * gg_67[k];

        t_99[k] = -f_3 * gf_s_44[k]
                  + f_2 * gh_s_99[k]
                  + f_4 * gf_43[k]
                  + pb_x[k] * gg_68[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_x, pb_y, gf_s_47, gh_s_100, gh_s_101, \
                         gh_s_102, gh_s_103, gf_46, gg_67, gg_69, gg_70, \
                         gg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_2 * gh_s_100[k]
                   + pb_y[k] * gg_67[k];

        t_101[k] = -f_3 * gf_s_47[k]
                   + f_2 * gh_s_101[k]
                   + f_4 * gf_46[k]
                   + pb_x[k] * gg_69[k];

        t_102[k] = f_2 * gh_s_102[k]
                   + pb_x[k] * gg_70[k];

        t_103[k] = f_2 * gh_s_103[k]
                   + pb_x[k] * gg_72[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_x, pb_y, gf_s_44, gf_s_45, gh_s_104, \
                         gh_s_105, gh_s_106, gf_43, gf_44, gg_70, gg_71, \
                         gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_2 * gh_s_104[k]
                   + pb_x[k] * gg_74[k];

        t_105[k] = -f_1 * gf_s_44[k]
                   + f_2 * gh_s_105[k]
                   + f_0 * gf_43[k]
                   + pb_y[k] * gg_70[k];

        t_106[k] = -f_11 * gf_s_45[k]
                   + f_2 * gh_s_106[k]
                   + f_7 * gf_44[k]
                   + pb_y[k] * gg_71[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_y, gf_s_46, gf_s_47, gh_s_107, gh_s_108, \
                         gh_s_109, gf_45, gf_46, gg_72, gg_73, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_5 * gf_s_46[k]
                   + f_2 * gh_s_107[k]
                   + f_6 * gf_45[k]
                   + pb_y[k] * gg_72[k];

        t_108[k] = -f_3 * gf_s_47[k]
                   + f_2 * gh_s_108[k]
                   + f_4 * gf_46[k]
                   + pb_y[k] * gg_73[k];

        t_109[k] = f_2 * gh_s_109[k]
                   + pb_y[k] * gg_74[k];
    }

#pragma omp simd aligned(t_110, pb_z, fg_40, gf_s_47, gh_s_110, gf_46, \
                         gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_0 * fg_40[k]
                   - f_1 * gf_s_47[k]
                   + f_2 * gh_s_110[k]
                   + f_0 * gf_46[k]
                   + pb_z[k] * gg_74[k];
    }
}

}  // namespace simdkin
