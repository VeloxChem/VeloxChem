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


#include "SimdElectronRepulsionVrrRecHG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_hg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);

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

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_3 = buffer.data(fg0 + 3);
    const auto *fg0_4 = buffer.data(fg0 + 4);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_6 = buffer.data(fg0 + 6);
    const auto *fg0_7 = buffer.data(fg0 + 7);
    const auto *fg0_8 = buffer.data(fg0 + 8);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_1 = buffer.data(fg1 + 1);
    const auto *fg1_2 = buffer.data(fg1 + 2);
    const auto *fg1_3 = buffer.data(fg1 + 3);
    const auto *fg1_4 = buffer.data(fg1 + 4);
    const auto *fg1_5 = buffer.data(fg1 + 5);
    const auto *fg1_6 = buffer.data(fg1 + 6);
    const auto *fg1_7 = buffer.data(fg1 + 7);
    const auto *fg1_8 = buffer.data(fg1 + 8);

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
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_111 = buffer.data(hf + 111);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_114 = buffer.data(hf + 114);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_123 = buffer.data(hf + 123);
    const auto *hf_124 = buffer.data(hf + 124);
    const auto *hf_125 = buffer.data(hf + 125);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_131 = buffer.data(hf + 131);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_133 = buffer.data(hf + 133);
    const auto *hf_134 = buffer.data(hf + 134);
    const auto *hf_135 = buffer.data(hf + 135);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_139 = buffer.data(hf + 139);
    const auto *hf_140 = buffer.data(hf + 140);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, \
                         hf_0, hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = pb_y[k] * hf_2[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, gf_3, gf_6, hd0_1, hd1_1, \
                         hf_3, hf_4, hf_5, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gf_3[k]
                 + pb_x[k] * hf_5[k];

        t_7[k] = pb_z[k] * hf_3[k];

        t_8[k] = pb_y[k] * hf_4[k];

        t_9[k] = f_0 * gf_6[k]
                 + pb_x[k] * hf_7[k];

        t_10[k] = f_1 * hd0_1[k]
                  - f_2 * hd1_1[k]
                  + pb_y[k] * hf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, gg_0, hd0_2, hd1_2, \
                         hf_5, hf_6, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * hf_5[k];

        t_12[k] = f_3 * hd0_2[k]
                  - f_4 * hd1_2[k]
                  + pb_y[k] * hf_6[k];

        t_13[k] = pb_y[k] * hf_7[k];

        t_14[k] = f_1 * hd0_2[k]
                  - f_2 * hd1_2[k]
                  + pb_z[k] * hf_7[k];

        t_15[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, gf_0, gf_1, gg_1, \
                         gg_2, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * gf_0[k]
                  + pb_y[k] * hf_8[k];

        t_17[k] = pb_z[k] * hf_8[k];

        t_18[k] = f_6 * gf_1[k]
                  + pa_y[k] * gg_1[k];

        t_19[k] = pb_z[k] * hf_9[k];

        t_20[k] = pa_y[k] * gg_2[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, gf_3, gf_8, gf_9, \
                         gg_4, gg_5, hf_10, hf_11, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * gf_8[k]
                  + pb_x[k] * hf_11[k];

        t_22[k] = pb_z[k] * hf_10[k];

        t_23[k] = f_7 * gf_9[k]
                  + pb_x[k] * hf_12[k];

        t_24[k] = pa_y[k] * gg_4[k];

        t_25[k] = f_7 * gf_3[k]
                  + pa_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, gf_5, gf_6, \
                         gg_0, gg_6, gg_7, hf_11, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * hf_11[k];

        t_27[k] = f_6 * gf_5[k]
                  + pa_y[k] * gg_6[k];

        t_28[k] = f_5 * gf_6[k]
                  + pb_y[k] * hf_13[k];

        t_29[k] = pa_y[k] * gg_7[k];

        t_30[k] = pa_z[k] * gg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, gf_0, gf_2, \
                         gg_1, gg_2, gg_3, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * hf_14[k];

        t_32[k] = f_5 * gf_0[k]
                  + pb_z[k] * hf_14[k];

        t_33[k] = pa_z[k] * gg_1[k];

        t_34[k] = pb_y[k] * hf_15[k];

        t_35[k] = f_6 * gf_2[k]
                  + pa_z[k] * gg_2[k];

        t_36[k] = pa_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, gf_14, gf_16, gg_5, hf_16, \
                         hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * gf_14[k]
                  + pb_x[k] * hf_18[k];

        t_38[k] = pb_y[k] * hf_16[k];

        t_39[k] = f_7 * gf_16[k]
                  + pb_x[k] * hf_19[k];

        t_40[k] = pa_z[k] * gg_5[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, gf_3, gf_4, gf_6, gg_6, \
                         gg_7, hf_17, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * gf_3[k]
                  + pb_z[k] * hf_17[k];

        t_42[k] = f_6 * gf_4[k]
                  + pa_z[k] * gg_6[k];

        t_43[k] = pb_y[k] * hf_19[k];

        t_44[k] = f_7 * gf_6[k]
                  + pa_z[k] * gg_7[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, fg0_0, fg1_0, gf_7, gg_8, \
                         hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_y[k] * gg_8[k];

        t_46[k] = f_6 * gf_7[k]
                  + pb_y[k] * hf_20[k];

        t_47[k] = pb_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, gf_19, gf_20, hd0_3, hd0_4, \
                         hd1_3, hd1_4, hf_21, hf_22, hf_23, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * gf_19[k]
                  + f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_x[k] * hf_23[k];

        t_49[k] = pb_z[k] * hf_21[k];

        t_50[k] = f_3 * hd0_3[k]
                  - f_4 * hd1_3[k]
                  + pb_z[k] * hf_22[k];

        t_51[k] = f_10 * gf_20[k]
                  + pb_x[k] * hf_24[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, fg0_3, fg1_3, gf_22, gf_23, \
                         gg_24, hf_23, hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * hf_23[k];

        t_53[k] = f_10 * gf_22[k]
                  + pb_x[k] * hf_26[k];

        t_54[k] = f_10 * gf_23[k]
                  + pb_x[k] * hf_27[k];

        t_55[k] = f_11 * fg0_3[k]
                  - f_12 * fg1_3[k]
                  + pa_x[k] * gg_24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, gf_10, hd0_4, hd0_5, hd1_4, \
                         hd1_5, hf_24, hf_25, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * hf_24[k];

        t_57[k] = f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_z[k] * hf_25[k];

        t_58[k] = f_6 * gf_10[k]
                  + pb_y[k] * hf_27[k];

        t_59[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_5[k]
                  + pb_z[k] * hf_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, gf_12, gg_9, \
                         gg_10, gg_13, gg_14, gg_15, hf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * gg_13[k];

        t_61[k] = pa_z[k] * gg_9[k];

        t_62[k] = pa_y[k] * gg_14[k];

        t_63[k] = pa_z[k] * gg_10[k];

        t_64[k] = f_5 * gf_12[k]
                  + pb_y[k] * hf_28[k];

        t_65[k] = pa_y[k] * gg_15[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, gf_26, gf_27, gg_11, \
                         gg_12, gg_16, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * gg_11[k];

        t_67[k] = f_10 * gf_26[k]
                  + pb_x[k] * hf_30[k];

        t_68[k] = f_10 * gf_27[k]
                  + pb_x[k] * hf_31[k];

        t_69[k] = pa_y[k] * gg_16[k];

        t_70[k] = pa_z[k] * gg_12[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, gf_8, gf_15, gf_16, gg_17, \
                         gg_18, hf_29, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * gf_8[k]
                  + pb_z[k] * hf_29[k];

        t_72[k] = f_6 * gf_15[k]
                  + pa_y[k] * gg_17[k];

        t_73[k] = f_5 * gf_16[k]
                  + pb_y[k] * hf_32[k];

        t_74[k] = pa_y[k] * gg_18[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, fg0_0, fg1_0, gf_11, gg_13, \
                         hd0_6, hd1_6, hf_33, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_z[k] * gg_13[k];

        t_76[k] = pb_y[k] * hf_33[k];

        t_77[k] = f_6 * gf_11[k]
                  + pb_z[k] * hf_33[k];

        t_78[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_y[k] * hf_34[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, gf_32, gf_33, gf_34, hd0_8, \
                         hd1_8, hf_35, hf_36, hf_37, hf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * hf_35[k];

        t_80[k] = f_10 * gf_32[k]
                  + f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_x[k] * hf_36[k];

        t_81[k] = f_10 * gf_33[k]
                  + pb_x[k] * hf_37[k];

        t_82[k] = f_10 * gf_34[k]
                  + pb_x[k] * hf_38[k];

        t_83[k] = pb_y[k] * hf_36[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, gf_13, gf_36, hd0_7, hd0_8, \
                         hd1_7, hd1_8, hf_37, hf_39, hf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_10 * gf_36[k]
                  + pb_x[k] * hf_40[k];

        t_85[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_y[k] * hf_37[k];

        t_86[k] = f_6 * gf_13[k]
                  + pb_z[k] * hf_37[k];

        t_87[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_39[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, fg0_1, fg0_4, fg1_1, fg1_4, \
                         gf_17, gg_19, gg_34, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * hf_40[k];

        t_89[k] = f_11 * fg0_4[k]
                  - f_12 * fg1_4[k]
                  + pa_x[k] * gg_34[k];

        t_90[k] = f_11 * fg0_1[k]
                  - f_12 * fg1_1[k]
                  + pa_y[k] * gg_19[k];

        t_91[k] = f_10 * gf_17[k]
                  + pb_y[k] * hf_41[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, gf_38, hd0_9, hd0_10, hd1_9, \
                         hd1_10, hf_41, hf_42, hf_43, hf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * hf_41[k];

        t_93[k] = f_6 * gf_38[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_44[k];

        t_94[k] = pb_z[k] * hf_42[k];

        t_95[k] = f_3 * hd0_9[k]
                  - f_4 * hd1_9[k]
                  + pb_z[k] * hf_43[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, gf_39, gf_40, gf_41, hf_44, \
                         hf_45, hf_47, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_6 * gf_39[k]
                  + pb_x[k] * hf_45[k];

        t_97[k] = pb_z[k] * hf_44[k];

        t_98[k] = f_6 * gf_40[k]
                  + pb_x[k] * hf_47[k];

        t_99[k] = f_6 * gf_41[k]
                  + pb_x[k] * hf_48[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, fg0_5, fg1_5, gf_23, \
                         gg_39, hd0_10, hd1_10, hf_45, hf_46, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_8 * fg0_5[k]
                   - f_9 * fg1_5[k]
                   + pa_x[k] * gg_39[k];

        t_101[k] = pb_z[k] * hf_45[k];

        t_102[k] = f_3 * hd0_10[k]
                   - f_4 * hd1_10[k]
                   + pb_z[k] * hf_46[k];

        t_103[k] = f_10 * gf_23[k]
                   + pb_y[k] * hf_48[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, gf_17, gg_19, gg_20, \
                         gg_21, hd0_11, hd1_11, hf_48, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * hd0_11[k]
                   - f_2 * hd1_11[k]
                   + pb_z[k] * hf_48[k];

        t_105[k] = pa_z[k] * gg_19[k];

        t_106[k] = pa_z[k] * gg_20[k];

        t_107[k] = f_5 * gf_17[k]
                   + pb_z[k] * hf_49[k];

        t_108[k] = pa_z[k] * gg_21[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, gf_18, gf_24, gf_44, \
                         gg_22, gg_23, hf_50, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * gf_24[k]
                   + pb_y[k] * hf_50[k];

        t_110[k] = f_6 * gf_18[k]
                   + pa_z[k] * gg_22[k];

        t_111[k] = pa_z[k] * gg_23[k];

        t_112[k] = f_6 * gf_44[k]
                   + pb_x[k] * hf_52[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, gf_20, gf_45, gf_46, \
                         gg_24, hf_51, hf_53, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_6 * gf_45[k]
                   + pb_x[k] * hf_53[k];

        t_114[k] = f_6 * gf_46[k]
                   + pb_x[k] * hf_54[k];

        t_115[k] = pa_z[k] * gg_24[k];

        t_116[k] = f_5 * gf_20[k]
                   + pb_z[k] * hf_51[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, gf_21, gf_23, gf_28, \
                         gg_25, gg_26, gg_27, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * gf_21[k]
                   + pa_z[k] * gg_25[k];

        t_118[k] = f_6 * gf_28[k]
                   + pb_y[k] * hf_54[k];

        t_119[k] = f_7 * gf_23[k]
                   + pa_z[k] * gg_26[k];

        t_120[k] = pa_y[k] * gg_27[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, gf_29, gf_30, gf_31, \
                         gg_28, gg_29, gg_30, hf_55, hf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * gf_29[k]
                   + pb_y[k] * hf_55[k];

        t_122[k] = pa_y[k] * gg_28[k];

        t_123[k] = f_6 * gf_30[k]
                   + pa_y[k] * gg_29[k];

        t_124[k] = f_5 * gf_31[k]
                   + pb_y[k] * hf_56[k];

        t_125[k] = pa_y[k] * gg_30[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, gf_33, gf_49, gf_50, \
                         gf_51, gg_31, gg_32, hf_57, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_6 * gf_49[k]
                   + pb_x[k] * hf_57[k];

        t_127[k] = f_6 * gf_50[k]
                   + pb_x[k] * hf_58[k];

        t_128[k] = f_6 * gf_51[k]
                   + pb_x[k] * hf_59[k];

        t_129[k] = pa_y[k] * gg_31[k];

        t_130[k] = f_7 * gf_33[k]
                   + pa_y[k] * gg_32[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, gf_25, gf_35, gf_36, \
                         gg_33, gg_34, hf_57, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * gf_25[k]
                   + pb_z[k] * hf_57[k];

        t_132[k] = f_6 * gf_35[k]
                   + pa_y[k] * gg_33[k];

        t_133[k] = f_5 * gf_36[k]
                   + pb_y[k] * hf_60[k];

        t_134[k] = pa_y[k] * gg_34[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, fg0_2, fg1_2, gf_29, \
                         gg_27, hd0_12, hd1_12, hf_61, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * fg0_2[k]
                   - f_12 * fg1_2[k]
                   + pa_z[k] * gg_27[k];

        t_136[k] = pb_y[k] * hf_61[k];

        t_137[k] = f_10 * gf_29[k]
                   + pb_z[k] * hf_61[k];

        t_138[k] = f_3 * hd0_12[k]
                   - f_4 * hd1_12[k]
                   + pb_y[k] * hf_62[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, gf_54, gf_55, gf_56, \
                         hd0_14, hd1_14, hf_63, hf_64, hf_65, hf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * hf_63[k];

        t_140[k] = f_6 * gf_54[k]
                   + f_3 * hd0_14[k]
                   - f_4 * hd1_14[k]
                   + pb_x[k] * hf_64[k];

        t_141[k] = f_6 * gf_55[k]
                   + pb_x[k] * hf_65[k];

        t_142[k] = f_6 * gf_56[k]
                   + pb_x[k] * hf_66[k];

        t_143[k] = pb_y[k] * hf_64[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, gf_33, gf_57, hd0_13, \
                         hd0_14, hd1_13, hd1_14, hf_65, hf_67, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_6 * gf_57[k]
                   + pb_x[k] * hf_68[k];

        t_145[k] = f_1 * hd0_13[k]
                   - f_2 * hd1_13[k]
                   + pb_y[k] * hf_65[k];

        t_146[k] = f_10 * gf_33[k]
                   + pb_z[k] * hf_65[k];

        t_147[k] = f_3 * hd0_14[k]
                   - f_4 * hd1_14[k]
                   + pb_y[k] * hf_67[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, pa_x, pb_y, pb_z, fg0_8, fg1_8, \
                         gf_37, gf_58, gg_44, gg_45, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * hf_68[k];

        t_149[k] = f_8 * fg0_8[k]
                   - f_9 * fg1_8[k]
                   + pa_x[k] * gg_44[k];

        t_150[k] = f_7 * gf_58[k]
                   + pa_x[k] * gg_45[k];

        t_151[k] = f_7 * gf_37[k]
                   + pb_y[k] * hf_69[k];

        t_152[k] = pb_z[k] * hf_69[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, gf_60, gf_61, \
                         gf_62, gg_47, gg_48, hf_70, hf_71, hf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * gf_60[k]
                   + pa_x[k] * gg_47[k];

        t_154[k] = pb_z[k] * hf_70[k];

        t_155[k] = f_6 * gf_61[k]
                   + pa_x[k] * gg_48[k];

        t_156[k] = f_5 * gf_62[k]
                   + pb_x[k] * hf_72[k];

        t_157[k] = pb_z[k] * hf_71[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_x, pb_x, pb_z, gf_64, gf_65, \
                         gg_49, gg_50, hf_72, hf_73, hf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * gf_64[k]
                   + pb_x[k] * hf_73[k];

        t_159[k] = f_5 * gf_65[k]
                   + pb_x[k] * hf_74[k];

        t_160[k] = pa_x[k] * gg_49[k];

        t_161[k] = pb_z[k] * hf_72[k];

        t_162[k] = pa_x[k] * gg_50[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pa_z, pb_z, gf_37, \
                         gg_35, gg_36, gg_37, gg_51, gg_52, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * gg_51[k];

        t_164[k] = pa_x[k] * gg_52[k];

        t_165[k] = pa_z[k] * gg_35[k];

        t_166[k] = pa_z[k] * gg_36[k];

        t_167[k] = f_5 * gf_37[k]
                   + pb_z[k] * hf_75[k];

        t_168[k] = pa_z[k] * gg_37[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_x, pa_z, pb_x, pb_y, gf_43, gf_68, \
                         gf_70, gg_38, gg_53, hf_76, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_10 * gf_43[k]
                   + pb_y[k] * hf_76[k];

        t_170[k] = f_6 * gf_68[k]
                   + pa_x[k] * gg_53[k];

        t_171[k] = pa_z[k] * gg_38[k];

        t_172[k] = f_5 * gf_70[k]
                   + pb_x[k] * hf_77[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, pa_x, pb_x, gf_71, gf_72, \
                         gg_54, gg_55, gg_56, gg_57, hf_78, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_5 * gf_71[k]
                   + pb_x[k] * hf_78[k];

        t_174[k] = f_5 * gf_72[k]
                   + pb_x[k] * hf_79[k];

        t_175[k] = pa_x[k] * gg_54[k];

        t_176[k] = pa_x[k] * gg_55[k];

        t_177[k] = pa_x[k] * gg_56[k];

        t_178[k] = pa_x[k] * gg_57[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pa_x, pb_y, pb_z, gf_42, gf_47, \
                         gf_73, gf_75, gg_58, gg_59, gg_60, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_x[k] * gg_58[k];

        t_180[k] = f_7 * gf_73[k]
                   + pa_x[k] * gg_59[k];

        t_181[k] = f_6 * gf_47[k]
                   + pb_y[k] * hf_80[k];

        t_182[k] = f_6 * gf_42[k]
                   + pb_z[k] * hf_80[k];

        t_183[k] = f_6 * gf_75[k]
                   + pa_x[k] * gg_60[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pb_x, pb_y, gf_48, gf_76, gf_77, \
                         gf_78, gg_61, hf_81, hf_82, hf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_6 * gf_48[k]
                   + pb_y[k] * hf_81[k];

        t_185[k] = f_6 * gf_76[k]
                   + pa_x[k] * gg_61[k];

        t_186[k] = f_5 * gf_77[k]
                   + pb_x[k] * hf_82[k];

        t_187[k] = f_5 * gf_78[k]
                   + pb_x[k] * hf_83[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, t_193, pa_x, pb_x, gf_79, gf_80, \
                         gg_62, gg_63, gg_64, gg_65, hf_84, hf_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_5 * gf_79[k]
                   + pb_x[k] * hf_84[k];

        t_189[k] = f_5 * gf_80[k]
                   + pb_x[k] * hf_85[k];

        t_190[k] = pa_x[k] * gg_62[k];

        t_191[k] = pa_x[k] * gg_63[k];

        t_192[k] = pa_x[k] * gg_64[k];

        t_193[k] = pa_x[k] * gg_65[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, pa_x, pa_y, pb_y, gf_52, gf_83, \
                         gg_40, gg_41, gg_66, gg_67, hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * gg_66[k];

        t_195[k] = pa_y[k] * gg_40[k];

        t_196[k] = f_5 * gf_52[k]
                   + pb_y[k] * hf_86[k];

        t_197[k] = pa_y[k] * gg_41[k];

        t_198[k] = f_6 * gf_83[k]
                   + pa_x[k] * gg_67[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pb_x, pb_y, gf_53, gf_84, gf_85, \
                         gg_42, hf_87, hf_88, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_5 * gf_53[k]
                   + pb_y[k] * hf_87[k];

        t_200[k] = pa_y[k] * gg_42[k];

        t_201[k] = f_5 * gf_84[k]
                   + pb_x[k] * hf_88[k];

        t_202[k] = f_5 * gf_85[k]
                   + pb_x[k] * hf_89[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, pa_x, pa_y, pb_x, gf_86, \
                         gg_43, gg_68, gg_69, gg_70, gg_71, hf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_5 * gf_86[k]
                   + pb_x[k] * hf_90[k];

        t_204[k] = pa_y[k] * gg_43[k];

        t_205[k] = pa_x[k] * gg_68[k];

        t_206[k] = pa_x[k] * gg_69[k];

        t_207[k] = pa_x[k] * gg_70[k];

        t_208[k] = pa_x[k] * gg_71[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pa_x, pb_y, pb_z, gf_52, gf_88, \
                         gf_91, gg_72, gg_73, gg_75, hf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_x[k] * gg_72[k];

        t_210[k] = f_7 * gf_88[k]
                   + pa_x[k] * gg_73[k];

        t_211[k] = pb_y[k] * hf_91[k];

        t_212[k] = f_7 * gf_52[k]
                   + pb_z[k] * hf_91[k];

        t_213[k] = f_6 * gf_91[k]
                   + pa_x[k] * gg_75[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pa_x, pb_x, pb_y, gf_92, gf_93, \
                         gf_94, gg_76, hf_92, hf_93, hf_94, hf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * hf_92[k];

        t_215[k] = f_6 * gf_92[k]
                   + pa_x[k] * gg_76[k];

        t_216[k] = f_5 * gf_93[k]
                   + pb_x[k] * hf_94[k];

        t_217[k] = f_5 * gf_94[k]
                   + pb_x[k] * hf_95[k];

        t_218[k] = pb_y[k] * hf_93[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, pa_x, pb_x, pb_y, gf_96, \
                         gg_77, gg_78, gg_79, gg_80, hf_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_5 * gf_96[k]
                   + pb_x[k] * hf_96[k];

        t_220[k] = pa_x[k] * gg_77[k];

        t_221[k] = pa_x[k] * gg_78[k];

        t_222[k] = pa_x[k] * gg_79[k];

        t_223[k] = pb_y[k] * hf_96[k];

        t_224[k] = pa_x[k] * gg_80[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, pb_x, pb_y, pb_z, gf_58, hd0_15, \
                         hd0_16, hd1_15, hd1_16, hf_97, hf_98, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_1 * hd0_15[k]
                   - f_2 * hd1_15[k]
                   + pb_x[k] * hf_97[k];

        t_226[k] = f_0 * gf_58[k]
                   + pb_y[k] * hf_97[k];

        t_227[k] = pb_z[k] * hf_97[k];

        t_228[k] = f_3 * hd0_16[k]
                   - f_4 * hd1_16[k]
                   + pb_x[k] * hf_99[k];

        t_229[k] = pb_z[k] * hf_98[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, pb_x, hd0_17, hd1_17, hf_100, \
                         hf_101, hf_102, hf_103, hf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_3 * hd0_17[k]
                   - f_4 * hd1_17[k]
                   + pb_x[k] * hf_100[k];

        t_231[k] = pb_x[k] * hf_101[k];

        t_232[k] = pb_x[k] * hf_102[k];

        t_233[k] = pb_x[k] * hf_103[k];

        t_234[k] = pb_x[k] * hf_104[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pb_y, pb_z, gf_62, gf_65, hd0_16, \
                         hd0_17, hd1_16, hd1_17, hf_101, hf_102, \
                         hf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_0 * gf_62[k]
                   + f_1 * hd0_16[k]
                   - f_2 * hd1_16[k]
                   + pb_y[k] * hf_101[k];

        t_236[k] = pb_z[k] * hf_101[k];

        t_237[k] = f_3 * hd0_16[k]
                   - f_4 * hd1_16[k]
                   + pb_z[k] * hf_102[k];

        t_238[k] = f_0 * gf_65[k]
                   + pb_y[k] * hf_104[k];

        t_239[k] = f_1 * hd0_17[k]
                   - f_2 * hd1_17[k]
                   + pb_z[k] * hf_104[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pa_z, pb_y, pb_z, gf_58, gf_67, \
                         gg_45, gg_46, gg_47, hf_105, hf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_z[k] * gg_45[k];

        t_241[k] = pa_z[k] * gg_46[k];

        t_242[k] = f_5 * gf_58[k]
                   + pb_z[k] * hf_105[k];

        t_243[k] = pa_z[k] * gg_47[k];

        t_244[k] = f_7 * gf_67[k]
                   + pb_y[k] * hf_106[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, t_250, pa_z, pb_x, gf_59, gg_48, \
                         gg_49, hf_107, hf_108, hf_109, hf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_6 * gf_59[k]
                   + pa_z[k] * gg_48[k];

        t_246[k] = pb_x[k] * hf_107[k];

        t_247[k] = pb_x[k] * hf_108[k];

        t_248[k] = pb_x[k] * hf_109[k];

        t_249[k] = pb_x[k] * hf_110[k];

        t_250[k] = pa_z[k] * gg_49[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_z, pb_y, pb_z, gf_62, gf_63, gf_65, \
                         gf_72, gg_50, gg_52, hf_107, hf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_5 * gf_62[k]
                   + pb_z[k] * hf_107[k];

        t_252[k] = f_6 * gf_63[k]
                   + pa_z[k] * gg_50[k];

        t_253[k] = f_7 * gf_72[k]
                   + pb_y[k] * hf_110[k];

        t_254[k] = f_7 * gf_65[k]
                   + pa_z[k] * gg_52[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pb_x, pb_y, pb_z, gf_66, gf_73, hd0_18, \
                         hd0_19, hd1_18, hd1_19, hf_111, hf_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_1 * hd0_18[k]
                   - f_2 * hd1_18[k]
                   + pb_x[k] * hf_111[k];

        t_256[k] = f_10 * gf_73[k]
                   + pb_y[k] * hf_111[k];

        t_257[k] = f_6 * gf_66[k]
                   + pb_z[k] * hf_111[k];

        t_258[k] = f_3 * hd0_19[k]
                   - f_4 * hd1_19[k]
                   + pb_x[k] * hf_113[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, pb_x, pb_y, gf_74, hd0_20, hd1_20, \
                         hf_112, hf_114, hf_115, hf_116, hf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_10 * gf_74[k]
                   + pb_y[k] * hf_112[k];

        t_260[k] = f_3 * hd0_20[k]
                   - f_4 * hd1_20[k]
                   + pb_x[k] * hf_114[k];

        t_261[k] = pb_x[k] * hf_115[k];

        t_262[k] = pb_x[k] * hf_116[k];

        t_263[k] = pb_x[k] * hf_117[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_z, pb_x, pb_z, fg0_5, fg1_5, gf_69, gg_54, \
                         hf_115, hf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = pb_x[k] * hf_118[k];

        t_265[k] = f_8 * fg0_5[k]
                   - f_9 * fg1_5[k]
                   + pa_z[k] * gg_54[k];

        t_266[k] = f_6 * gf_69[k]
                   + pb_z[k] * hf_115[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_y, pb_y, fg0_7, fg1_7, gf_79, gf_80, gg_66, \
                         hd0_20, hd1_20, hf_117, hf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_10 * gf_79[k]
                   + f_3 * hd0_20[k]
                   - f_4 * hd1_20[k]
                   + pb_y[k] * hf_117[k];

        t_268[k] = f_10 * gf_80[k]
                   + pb_y[k] * hf_118[k];

        t_269[k] = f_11 * fg0_7[k]
                   - f_12 * fg1_7[k]
                   + pa_y[k] * gg_66[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_x, pb_y, pb_z, gf_73, gf_81, hd0_21, \
                         hd0_22, hd1_21, hd1_22, hf_119, hf_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * hd0_21[k]
                   - f_2 * hd1_21[k]
                   + pb_x[k] * hf_119[k];

        t_271[k] = f_6 * gf_81[k]
                   + pb_y[k] * hf_119[k];

        t_272[k] = f_10 * gf_73[k]
                   + pb_z[k] * hf_119[k];

        t_273[k] = f_3 * hd0_22[k]
                   - f_4 * hd1_22[k]
                   + pb_x[k] * hf_121[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_x, pb_y, gf_82, hd0_23, hd1_23, \
                         hf_120, hf_122, hf_123, hf_124, hf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_6 * gf_82[k]
                   + pb_y[k] * hf_120[k];

        t_275[k] = f_3 * hd0_23[k]
                   - f_4 * hd1_23[k]
                   + pb_x[k] * hf_122[k];

        t_276[k] = pb_x[k] * hf_123[k];

        t_277[k] = pb_x[k] * hf_124[k];

        t_278[k] = pb_x[k] * hf_125[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pa_z, pb_x, pb_z, fg0_6, fg1_6, gf_77, gg_62, \
                         hf_123, hf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_x[k] * hf_126[k];

        t_280[k] = f_11 * fg0_6[k]
                   - f_12 * fg1_6[k]
                   + pa_z[k] * gg_62[k];

        t_281[k] = f_10 * gf_77[k]
                   + pb_z[k] * hf_123[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_y, pb_y, fg0_8, fg1_8, gf_86, gf_87, \
                         gg_72, gg_73, hd0_23, hd1_23, hf_125, hf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_6 * gf_86[k]
                   + f_3 * hd0_23[k]
                   - f_4 * hd1_23[k]
                   + pb_y[k] * hf_125[k];

        t_283[k] = f_6 * gf_87[k]
                   + pb_y[k] * hf_126[k];

        t_284[k] = f_8 * fg0_8[k]
                   - f_9 * fg1_8[k]
                   + pa_y[k] * gg_72[k];

        t_285[k] = pa_y[k] * gg_73[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, pa_y, pb_y, gf_88, gf_89, gf_90, \
                         gg_74, gg_75, gg_76, hf_127, hf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_5 * gf_88[k]
                   + pb_y[k] * hf_127[k];

        t_287[k] = pa_y[k] * gg_74[k];

        t_288[k] = f_6 * gf_89[k]
                   + pa_y[k] * gg_75[k];

        t_289[k] = f_5 * gf_90[k]
                   + pb_y[k] * hf_128[k];

        t_290[k] = pa_y[k] * gg_76[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, t_295, t_296, pa_y, pb_x, pb_z, gf_84, \
                         gf_93, gg_77, hf_129, hf_130, hf_131, hf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = pb_x[k] * hf_129[k];

        t_292[k] = pb_x[k] * hf_130[k];

        t_293[k] = pb_x[k] * hf_131[k];

        t_294[k] = pb_x[k] * hf_132[k];

        t_295[k] = f_7 * gf_93[k]
                   + pa_y[k] * gg_77[k];

        t_296[k] = f_7 * gf_84[k]
                   + pb_z[k] * hf_129[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, pa_y, pb_x, pb_y, gf_95, gf_96, \
                         gg_79, gg_80, hd0_24, hd1_24, hf_132, hf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_6 * gf_95[k]
                   + pa_y[k] * gg_79[k];

        t_298[k] = f_5 * gf_96[k]
                   + pb_y[k] * hf_132[k];

        t_299[k] = pa_y[k] * gg_80[k];

        t_300[k] = f_1 * hd0_24[k]
                   - f_2 * hd1_24[k]
                   + pb_x[k] * hf_133[k];

        t_301[k] = pb_y[k] * hf_133[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pb_x, pb_y, pb_z, gf_88, hd0_25, hd0_26, \
                         hd1_25, hd1_26, hf_133, hf_134, hf_135, \
                         hf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * gf_88[k]
                   + pb_z[k] * hf_133[k];

        t_303[k] = f_3 * hd0_25[k]
                   - f_4 * hd1_25[k]
                   + pb_x[k] * hf_135[k];

        t_304[k] = pb_y[k] * hf_134[k];

        t_305[k] = f_3 * hd0_26[k]
                   - f_4 * hd1_26[k]
                   + pb_x[k] * hf_136[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, pb_x, pb_y, pb_z, gf_93, \
                         hd0_25, hd1_25, hf_137, hf_138, hf_139, \
                         hf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = pb_x[k] * hf_137[k];

        t_307[k] = pb_x[k] * hf_138[k];

        t_308[k] = pb_x[k] * hf_139[k];

        t_309[k] = pb_x[k] * hf_140[k];

        t_310[k] = f_1 * hd0_25[k]
                   - f_2 * hd1_25[k]
                   + pb_y[k] * hf_137[k];

        t_311[k] = f_0 * gf_93[k]
                   + pb_z[k] * hf_137[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_y, pb_z, gf_96, hd0_26, hd1_26, hf_139, \
                         hf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_3 * hd0_26[k]
                   - f_4 * hd1_26[k]
                   + pb_y[k] * hf_139[k];

        t_313[k] = pb_y[k] * hf_140[k];

        t_314[k] = f_0 * gf_96[k]
                   + f_1 * hd0_26[k]
                   - f_2 * hd1_26[k]
                   + pb_z[k] * hf_140[k];
    }
}

auto
compute_prim_hg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_3 = buffer.data(fg0 + 3);
    const auto *fg0_4 = buffer.data(fg0 + 4);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_6 = buffer.data(fg0 + 6);
    const auto *fg0_7 = buffer.data(fg0 + 7);
    const auto *fg0_8 = buffer.data(fg0 + 8);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_9 = buffer.data(fg1 + 9);
    const auto *fg1_12 = buffer.data(fg1 + 12);
    const auto *fg1_20 = buffer.data(fg1 + 20);
    const auto *fg1_25 = buffer.data(fg1 + 25);
    const auto *fg1_31 = buffer.data(fg1 + 31);
    const auto *fg1_37 = buffer.data(fg1 + 37);
    const auto *fg1_47 = buffer.data(fg1 + 47);
    const auto *fg1_58 = buffer.data(fg1 + 58);

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
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
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
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_32 = buffer.data(gg + 32);
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
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
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
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_85 = buffer.data(gg + 85);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, gf_3, gf_6, hd0_1, hd0_2, hd1_1, \
                         hd1_2, hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * gf_3[k]
                 + pb_x[k] * hf_3[k];

        t_6[k] = f_0 * gf_6[k]
                 + pb_x[k] * hf_5[k];

        t_7[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_8[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_9[k] = pb_y[k] * hf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, gf_0, gf_1, gg_0, gg_3, \
                         hd0_2, hd1_2, hf_5, hf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * hd0_2[k]
                  - f_2 * hd1_2[k]
                  + pb_z[k] * hf_5[k];

        t_11[k] = pa_y[k] * gg_0[k];

        t_12[k] = f_5 * gf_0[k]
                  + pb_y[k] * hf_6[k];

        t_13[k] = f_6 * gf_1[k]
                  + pa_y[k] * gg_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pb_x, gf_3, gf_5, gf_8, gg_4, gg_5, \
                         gg_6, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * gg_4[k];

        t_15[k] = f_7 * gf_8[k]
                  + pb_x[k] * hf_7[k];

        t_16[k] = f_7 * gf_3[k]
                  + pa_y[k] * gg_5[k];

        t_17[k] = f_6 * gf_5[k]
                  + pa_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, gf_0, gf_6, \
                         gg_0, gg_3, gg_8, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gf_6[k]
                  + pb_y[k] * hf_8[k];

        t_19[k] = pa_y[k] * gg_8[k];

        t_20[k] = pa_z[k] * gg_0[k];

        t_21[k] = f_5 * gf_0[k]
                  + pb_z[k] * hf_9[k];

        t_22[k] = pa_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_z, gf_2, gf_3, gf_13, gg_4, \
                         gg_5, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * gf_2[k]
                  + pa_z[k] * gg_4[k];

        t_24[k] = f_7 * gf_13[k]
                  + pb_x[k] * hf_11[k];

        t_25[k] = pa_z[k] * gg_5[k];

        t_26[k] = f_5 * gf_3[k]
                  + pb_z[k] * hf_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, fg0_0, fg1_0, gf_4, gf_6, \
                         gf_7, gg_6, gg_8, gg_9, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gf_4[k]
                  + pa_z[k] * gg_6[k];

        t_28[k] = f_7 * gf_6[k]
                  + pa_z[k] * gg_8[k];

        t_29[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_y[k] * gg_9[k];

        t_30[k] = f_6 * gf_7[k]
                  + pb_y[k] * hf_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, pb_z, gf_16, gf_17, hd0_3, hd0_4, \
                         hd1_3, hd1_4, hf_12, hf_13, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_z[k] * hf_12[k];

        t_32[k] = f_10 * gf_16[k]
                  + f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_x[k] * hf_14[k];

        t_33[k] = f_3 * hd0_3[k]
                  - f_4 * hd1_3[k]
                  + pb_z[k] * hf_13[k];

        t_34[k] = f_10 * gf_17[k]
                  + pb_x[k] * hf_15[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, pb_z, fg0_3, fg1_20, gf_9, gg_22, \
                         hd0_4, hd1_4, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_11 * fg0_3[k]
                  - f_12 * fg1_20[k]
                  + pa_x[k] * gg_22[k];

        t_36[k] = pb_z[k] * hf_15[k];

        t_37[k] = f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_z[k] * hf_16[k];

        t_38[k] = f_6 * gf_9[k]
                  + pb_y[k] * hf_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pa_z, pb_z, gg_10, gg_11, gg_13, \
                         gg_14, hd0_5, hd1_5, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_5[k]
                  + pb_z[k] * hf_17[k];

        t_40[k] = pa_y[k] * gg_13[k];

        t_41[k] = pa_z[k] * gg_10[k];

        t_42[k] = pa_y[k] * gg_14[k];

        t_43[k] = pa_z[k] * gg_11[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, gf_8, gf_12, gf_13, gg_15, \
                         gg_16, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * gf_8[k]
                  + pb_z[k] * hf_18[k];

        t_45[k] = f_6 * gf_12[k]
                  + pa_y[k] * gg_15[k];

        t_46[k] = f_5 * gf_13[k]
                  + pb_y[k] * hf_19[k];

        t_47[k] = pa_y[k] * gg_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, fg0_0, fg1_0, gf_10, gg_12, \
                         hd0_6, hd1_6, hf_20, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_z[k] * gg_12[k];

        t_49[k] = pb_y[k] * hf_20[k];

        t_50[k] = f_6 * gf_10[k]
                  + pb_z[k] * hf_20[k];

        t_51[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_y[k] * hf_21[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, gf_24, gf_27, hd0_7, hd0_8, hd1_7, \
                         hd1_8, hf_22, hf_23, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * gf_24[k]
                  + f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_x[k] * hf_22[k];

        t_53[k] = f_10 * gf_27[k]
                  + pb_x[k] * hf_25[k];

        t_54[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pb_y, pb_z, fg0_4, fg1_25, gf_11, \
                         gg_35, hd0_8, hd1_8, hf_23, hf_24, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_6 * gf_11[k]
                  + pb_z[k] * hf_23[k];

        t_56[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_24[k];

        t_57[k] = pb_y[k] * hf_25[k];

        t_58[k] = f_11 * fg0_4[k]
                  - f_12 * fg1_25[k]
                  + pa_x[k] * gg_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_y, pb_y, pb_z, fg0_1, fg1_9, gf_14, gg_17, \
                         hf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_11 * fg0_1[k]
                  - f_12 * fg1_9[k]
                  + pa_y[k] * gg_17[k];

        t_60[k] = f_10 * gf_14[k]
                  + pb_y[k] * hf_26[k];

        t_61[k] = pb_z[k] * hf_26[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pb_x, pb_z, gf_29, gf_30, hd0_9, hd0_10, hd1_9, \
                         hd1_10, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_6 * gf_29[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_28[k];

        t_63[k] = f_3 * hd0_9[k]
                  - f_4 * hd1_9[k]
                  + pb_z[k] * hf_27[k];

        t_64[k] = f_6 * gf_30[k]
                  + pb_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_y, pb_z, fg0_5, fg1_31, gf_19, \
                         gg_38, hd0_10, hd1_10, hf_29, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_8 * fg0_5[k]
                  - f_9 * fg1_31[k]
                  + pa_x[k] * gg_38[k];

        t_66[k] = pb_z[k] * hf_29[k];

        t_67[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_z[k] * hf_30[k];

        t_68[k] = f_10 * gf_19[k]
                  + pb_y[k] * hf_31[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_z, pb_z, gf_14, gf_15, gg_17, gg_19, \
                         gg_20, hd0_11, hd1_11, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * hd0_11[k]
                  - f_2 * hd1_11[k]
                  + pb_z[k] * hf_31[k];

        t_70[k] = pa_z[k] * gg_17[k];

        t_71[k] = f_5 * gf_14[k]
                  + pb_z[k] * hf_32[k];

        t_72[k] = pa_z[k] * gg_19[k];

        t_73[k] = f_6 * gf_15[k]
                  + pa_z[k] * gg_20[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pb_y, pb_z, gf_17, gf_18, gf_21, gg_22, \
                         gg_24, hf_33, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_z[k] * gg_22[k];

        t_75[k] = f_5 * gf_17[k]
                  + pb_z[k] * hf_33[k];

        t_76[k] = f_6 * gf_18[k]
                  + pa_z[k] * gg_24[k];

        t_77[k] = f_6 * gf_21[k]
                  + pb_y[k] * hf_34[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, gf_19, gf_23, gg_25, gg_26, \
                         gg_28, gg_29, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_7 * gf_19[k]
                  + pa_z[k] * gg_25[k];

        t_79[k] = pa_y[k] * gg_26[k];

        t_80[k] = pa_y[k] * gg_28[k];

        t_81[k] = f_6 * gf_23[k]
                  + pa_y[k] * gg_29[k];

        t_82[k] = pa_y[k] * gg_30[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_y, pb_z, gf_20, gf_25, gf_26, gf_27, \
                         gg_32, gg_33, hf_35, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_7 * gf_25[k]
                  + pa_y[k] * gg_32[k];

        t_84[k] = f_6 * gf_20[k]
                  + pb_z[k] * hf_35[k];

        t_85[k] = f_6 * gf_26[k]
                  + pa_y[k] * gg_33[k];

        t_86[k] = f_5 * gf_27[k]
                  + pb_y[k] * hf_36[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_y, pb_z, fg0_2, fg1_12, gf_22, \
                         gg_26, gg_35, hf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_y[k] * gg_35[k];

        t_88[k] = f_11 * fg0_2[k]
                  - f_12 * fg1_12[k]
                  + pa_z[k] * gg_26[k];

        t_89[k] = pb_y[k] * hf_37[k];

        t_90[k] = f_10 * gf_22[k]
                  + pb_z[k] * hf_37[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_y, gf_33, gf_34, hd0_12, hd0_14, hd1_12, \
                         hd1_14, hf_38, hf_39, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_y[k] * hf_38[k];

        t_92[k] = f_6 * gf_33[k]
                  + f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_x[k] * hf_39[k];

        t_93[k] = f_6 * gf_34[k]
                  + pb_x[k] * hf_42[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_y, pb_z, gf_25, hd0_13, hd0_14, hd1_13, \
                         hd1_14, hf_40, hf_41, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_y[k] * hf_40[k];

        t_95[k] = f_10 * gf_25[k]
                  + pb_z[k] * hf_40[k];

        t_96[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_41[k];

        t_97[k] = pb_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_x, pb_y, fg0_8, fg1_58, gf_28, gf_35, \
                         gf_37, gg_42, gg_43, gg_44, hf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_8 * fg0_8[k]
                  - f_9 * fg1_58[k]
                  + pa_x[k] * gg_42[k];

        t_99[k] = f_7 * gf_35[k]
                  + pa_x[k] * gg_43[k];

        t_100[k] = f_7 * gf_28[k]
                   + pb_y[k] * hf_43[k];

        t_101[k] = f_6 * gf_37[k]
                   + pa_x[k] * gg_44[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pa_x, pb_x, gf_38, gf_39, \
                         gg_45, gg_48, gg_50, gg_51, gg_52, hf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_6 * gf_38[k]
                   + pa_x[k] * gg_45[k];

        t_103[k] = f_5 * gf_39[k]
                   + pb_x[k] * hf_44[k];

        t_104[k] = pa_x[k] * gg_48[k];

        t_105[k] = pa_x[k] * gg_50[k];

        t_106[k] = pa_x[k] * gg_51[k];

        t_107[k] = pa_x[k] * gg_52[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_x, pa_z, pb_z, gf_28, gf_43, \
                         gg_36, gg_37, gg_53, gg_55, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * gg_36[k];

        t_109[k] = f_5 * gf_28[k]
                   + pb_z[k] * hf_45[k];

        t_110[k] = pa_z[k] * gg_37[k];

        t_111[k] = f_6 * gf_43[k]
                   + pa_x[k] * gg_53[k];

        t_112[k] = pa_x[k] * gg_55[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_x, pb_z, gf_31, gf_46, gg_56, \
                         gg_57, gg_58, gg_59, hf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_x[k] * gg_56[k];

        t_114[k] = pa_x[k] * gg_57[k];

        t_115[k] = pa_x[k] * gg_58[k];

        t_116[k] = f_7 * gf_46[k]
                   + pa_x[k] * gg_59[k];

        t_117[k] = f_6 * gf_31[k]
                   + pb_z[k] * hf_46[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, t_123, pa_x, gf_47, gf_48, gg_60, \
                         gg_61, gg_64, gg_65, gg_66, gg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_6 * gf_47[k]
                   + pa_x[k] * gg_60[k];

        t_119[k] = f_6 * gf_48[k]
                   + pa_x[k] * gg_61[k];

        t_120[k] = pa_x[k] * gg_64[k];

        t_121[k] = pa_x[k] * gg_65[k];

        t_122[k] = pa_x[k] * gg_66[k];

        t_123[k] = pa_x[k] * gg_67[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, pa_x, pa_y, gf_52, gg_39, \
                         gg_40, gg_41, gg_68, gg_69, gg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_x[k] * gg_68[k];

        t_125[k] = pa_y[k] * gg_39[k];

        t_126[k] = pa_y[k] * gg_40[k];

        t_127[k] = f_6 * gf_52[k]
                   + pa_x[k] * gg_69[k];

        t_128[k] = pa_y[k] * gg_41[k];

        t_129[k] = pa_x[k] * gg_70[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_x, pb_z, gf_32, gf_56, gg_71, \
                         gg_72, gg_73, gg_75, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_x[k] * gg_71[k];

        t_131[k] = pa_x[k] * gg_72[k];

        t_132[k] = pa_x[k] * gg_73[k];

        t_133[k] = f_7 * gf_56[k]
                   + pa_x[k] * gg_75[k];

        t_134[k] = f_7 * gf_32[k]
                   + pb_z[k] * hf_47[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pa_x, pb_x, gf_58, gf_59, gf_62, \
                         gg_77, gg_78, gg_81, gg_82, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_6 * gf_58[k]
                   + pa_x[k] * gg_77[k];

        t_136[k] = f_6 * gf_59[k]
                   + pa_x[k] * gg_78[k];

        t_137[k] = f_5 * gf_62[k]
                   + pb_x[k] * hf_48[k];

        t_138[k] = pa_x[k] * gg_81[k];

        t_139[k] = pa_x[k] * gg_82[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pb_x, pb_y, gf_35, gg_83, gg_85, \
                         hd0_15, hd1_15, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pa_x[k] * gg_83[k];

        t_141[k] = pa_x[k] * gg_85[k];

        t_142[k] = f_1 * hd0_15[k]
                   - f_2 * hd1_15[k]
                   + pb_x[k] * hf_49[k];

        t_143[k] = f_0 * gf_35[k]
                   + pb_y[k] * hf_49[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pb_x, pb_y, gf_39, hd0_16, hd0_17, \
                         hd1_16, hd1_17, hf_50, hf_51, hf_52, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_3 * hd0_16[k]
                   - f_4 * hd1_16[k]
                   + pb_x[k] * hf_50[k];

        t_145[k] = f_3 * hd0_17[k]
                   - f_4 * hd1_17[k]
                   + pb_x[k] * hf_51[k];

        t_146[k] = pb_x[k] * hf_52[k];

        t_147[k] = pb_x[k] * hf_54[k];

        t_148[k] = f_0 * gf_39[k]
                   + f_1 * hd0_16[k]
                   - f_2 * hd1_16[k]
                   + pb_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pb_y, pb_z, gf_41, hd0_16, hd0_17, \
                         hd1_16, hd1_17, hf_52, hf_53, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_z[k] * hf_52[k];

        t_150[k] = f_3 * hd0_16[k]
                   - f_4 * hd1_16[k]
                   + pb_z[k] * hf_53[k];

        t_151[k] = f_0 * gf_41[k]
                   + pb_y[k] * hf_54[k];

        t_152[k] = f_1 * hd0_17[k]
                   - f_2 * hd1_17[k]
                   + pb_z[k] * hf_54[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pa_z, pb_z, gf_35, gf_36, gg_43, \
                         gg_44, gg_45, gg_48, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * gg_43[k];

        t_154[k] = f_5 * gf_35[k]
                   + pb_z[k] * hf_55[k];

        t_155[k] = pa_z[k] * gg_44[k];

        t_156[k] = f_6 * gf_36[k]
                   + pa_z[k] * gg_45[k];

        t_157[k] = pa_z[k] * gg_48[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pb_y, pb_z, gf_39, gf_40, gf_41, \
                         gf_45, gg_50, gg_52, hf_56, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * gf_39[k]
                   + pb_z[k] * hf_56[k];

        t_159[k] = f_6 * gf_40[k]
                   + pa_z[k] * gg_50[k];

        t_160[k] = f_7 * gf_45[k]
                   + pb_y[k] * hf_57[k];

        t_161[k] = f_7 * gf_41[k]
                   + pa_z[k] * gg_52[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_x, pb_z, gf_42, hd0_18, hd0_19, \
                         hd0_20, hd1_18, hd1_19, hd1_20, hf_58, hf_59, \
                         hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_1 * hd0_18[k]
                   - f_2 * hd1_18[k]
                   + pb_x[k] * hf_58[k];

        t_163[k] = f_6 * gf_42[k]
                   + pb_z[k] * hf_58[k];

        t_164[k] = f_3 * hd0_19[k]
                   - f_4 * hd1_19[k]
                   + pb_x[k] * hf_59[k];

        t_165[k] = f_3 * hd0_20[k]
                   - f_4 * hd1_20[k]
                   + pb_x[k] * hf_60[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pb_x, pb_z, fg0_5, fg1_31, gf_44, \
                         gg_54, hf_61, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_x[k] * hf_61[k];

        t_167[k] = pb_x[k] * hf_63[k];

        t_168[k] = f_8 * fg0_5[k]
                   - f_9 * fg1_31[k]
                   + pa_z[k] * gg_54[k];

        t_169[k] = f_6 * gf_44[k]
                   + pb_z[k] * hf_61[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pa_y, pb_y, fg0_7, fg1_47, gf_50, gf_51, gg_68, \
                         hd0_20, hd1_20, hf_62, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * gf_50[k]
                   + f_3 * hd0_20[k]
                   - f_4 * hd1_20[k]
                   + pb_y[k] * hf_62[k];

        t_171[k] = f_10 * gf_51[k]
                   + pb_y[k] * hf_63[k];

        t_172[k] = f_11 * fg0_7[k]
                   - f_12 * fg1_47[k]
                   + pa_y[k] * gg_68[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pb_x, pb_z, gf_46, hd0_21, hd0_22, \
                         hd0_23, hd1_21, hd1_22, hd1_23, hf_64, hf_65, \
                         hf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_1 * hd0_21[k]
                   - f_2 * hd1_21[k]
                   + pb_x[k] * hf_64[k];

        t_174[k] = f_10 * gf_46[k]
                   + pb_z[k] * hf_64[k];

        t_175[k] = f_3 * hd0_22[k]
                   - f_4 * hd1_22[k]
                   + pb_x[k] * hf_65[k];

        t_176[k] = f_3 * hd0_23[k]
                   - f_4 * hd1_23[k]
                   + pb_x[k] * hf_66[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_z, pb_x, pb_z, fg0_6, fg1_37, gf_49, \
                         gg_64, hf_67, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pb_x[k] * hf_67[k];

        t_178[k] = pb_x[k] * hf_69[k];

        t_179[k] = f_11 * fg0_6[k]
                   - f_12 * fg1_37[k]
                   + pa_z[k] * gg_64[k];

        t_180[k] = f_10 * gf_49[k]
                   + pb_z[k] * hf_67[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_y, pb_y, fg0_8, fg1_58, gf_54, gf_55, \
                         gg_74, gg_75, hd0_23, hd1_23, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * gf_54[k]
                   + f_3 * hd0_23[k]
                   - f_4 * hd1_23[k]
                   + pb_y[k] * hf_68[k];

        t_182[k] = f_6 * gf_55[k]
                   + pb_y[k] * hf_69[k];

        t_183[k] = f_8 * fg0_8[k]
                   - f_9 * fg1_58[k]
                   + pa_y[k] * gg_74[k];

        t_184[k] = pa_y[k] * gg_75[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, pa_y, pb_z, gf_53, gf_57, gf_60, \
                         gg_76, gg_77, gg_78, gg_81, hf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pa_y[k] * gg_76[k];

        t_186[k] = f_6 * gf_57[k]
                   + pa_y[k] * gg_77[k];

        t_187[k] = pa_y[k] * gg_78[k];

        t_188[k] = f_7 * gf_60[k]
                   + pa_y[k] * gg_81[k];

        t_189[k] = f_7 * gf_53[k]
                   + pb_z[k] * hf_70[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_y, pb_x, pb_y, gf_61, gf_62, gg_83, \
                         gg_85, hd0_24, hd1_24, hf_71, hf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_6 * gf_61[k]
                   + pa_y[k] * gg_83[k];

        t_191[k] = f_5 * gf_62[k]
                   + pb_y[k] * hf_71[k];

        t_192[k] = pa_y[k] * gg_85[k];

        t_193[k] = f_1 * hd0_24[k]
                   - f_2 * hd1_24[k]
                   + pb_x[k] * hf_72[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pb_x, pb_z, gf_56, hd0_25, hd0_26, \
                         hd1_25, hd1_26, hf_72, hf_73, hf_74, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_0 * gf_56[k]
                   + pb_z[k] * hf_72[k];

        t_195[k] = f_3 * hd0_25[k]
                   - f_4 * hd1_25[k]
                   + pb_x[k] * hf_73[k];

        t_196[k] = f_3 * hd0_26[k]
                   - f_4 * hd1_26[k]
                   + pb_x[k] * hf_74[k];

        t_197[k] = pb_x[k] * hf_75[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pb_x, pb_y, pb_z, gf_60, hd0_25, \
                         hd0_26, hd1_25, hd1_26, hf_75, hf_76, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_x[k] * hf_77[k];

        t_199[k] = f_1 * hd0_25[k]
                   - f_2 * hd1_25[k]
                   + pb_y[k] * hf_75[k];

        t_200[k] = f_0 * gf_60[k]
                   + pb_z[k] * hf_75[k];

        t_201[k] = f_3 * hd0_26[k]
                   - f_4 * hd1_26[k]
                   + pb_y[k] * hf_76[k];

        t_202[k] = pb_y[k] * hf_77[k];
    }

#pragma omp simd aligned(t_203, pb_z, gf_62, hd0_26, hd1_26, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * gf_62[k]
                   + f_1 * hd0_26[k]
                   - f_2 * hd1_26[k]
                   + pb_z[k] * hf_77[k];
    }
}

auto
compute_prim_hg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_3 = buffer.data(fg0 + 3);
    const auto *fg0_4 = buffer.data(fg0 + 4);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_6 = buffer.data(fg0 + 6);
    const auto *fg0_7 = buffer.data(fg0 + 7);
    const auto *fg0_8 = buffer.data(fg0 + 8);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_1 = buffer.data(fg1 + 1);
    const auto *fg1_2 = buffer.data(fg1 + 2);
    const auto *fg1_3 = buffer.data(fg1 + 3);
    const auto *fg1_4 = buffer.data(fg1 + 4);
    const auto *fg1_5 = buffer.data(fg1 + 5);
    const auto *fg1_6 = buffer.data(fg1 + 6);
    const auto *fg1_7 = buffer.data(fg1 + 7);
    const auto *fg1_8 = buffer.data(fg1 + 8);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_41 = buffer.data(gf + 41);

    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, hd0_1, hd0_2, hd1_1, hd1_2, hf_3, \
                         hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, fg0_0, fg1_0, gf_8, gg_9, hd0_4, \
                         hd1_4, hf_6, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fg0_0[k]
                 - f_6 * fg1_0[k]
                 + pa_y[k] * gg_9[k];

        t_10[k] = pb_z[k] * hf_6[k];

        t_11[k] = f_7 * gf_8[k]
                  + f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_x[k] * hf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, fg0_3, fg1_3, gf_9, gg_16, \
                         hd0_3, hd1_3, hf_7, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hd0_3[k]
                  - f_4 * hd1_3[k]
                  + pb_z[k] * hf_7[k];

        t_13[k] = f_7 * gf_9[k]
                  + pb_x[k] * hf_9[k];

        t_14[k] = f_8 * fg0_3[k]
                  - f_9 * fg1_3[k]
                  + pa_x[k] * gg_16[k];

        t_15[k] = pb_z[k] * hf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, fg0_0, fg1_0, gg_10, hd0_4, hd0_5, \
                         hd1_4, hd1_5, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_z[k] * hf_10[k];

        t_17[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_5[k]
                  + pb_z[k] * hf_11[k];

        t_18[k] = f_5 * fg0_0[k]
                  - f_6 * fg1_0[k]
                  + pa_z[k] * gg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, gf_14, gf_17, hd0_6, hd0_8, \
                         hd1_6, hd1_8, hf_12, hf_13, hf_14, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * hf_12[k];

        t_20[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_y[k] * hf_13[k];

        t_21[k] = f_7 * gf_14[k]
                  + f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_x[k] * hf_14[k];

        t_22[k] = f_7 * gf_17[k]
                  + pb_x[k] * hf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, fg0_4, fg1_4, gg_28, hd0_7, \
                         hd0_8, hd1_7, hd1_8, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_y[k] * hf_15[k];

        t_24[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_16[k];

        t_25[k] = pb_y[k] * hf_17[k];

        t_26[k] = f_8 * fg0_4[k]
                  - f_9 * fg1_4[k]
                  + pa_x[k] * gg_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, fg0_1, fg1_1, gf_18, gg_11, \
                         hd0_10, hd1_10, hf_18, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_8 * fg0_1[k]
                  - f_9 * fg1_1[k]
                  + pa_y[k] * gg_11[k];

        t_28[k] = pb_z[k] * hf_18[k];

        t_29[k] = f_10 * gf_18[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, fg0_5, fg1_5, gf_19, gg_29, \
                         hd0_9, hd1_9, hf_19, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hd0_9[k]
                  - f_4 * hd1_9[k]
                  + pb_z[k] * hf_19[k];

        t_31[k] = f_10 * gf_19[k]
                  + pb_x[k] * hf_21[k];

        t_32[k] = f_5 * fg0_5[k]
                  - f_6 * fg1_5[k]
                  + pa_x[k] * gg_29[k];

        t_33[k] = pb_z[k] * hf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, fg0_2, fg1_2, gg_20, hd0_10, hd0_11, \
                         hd1_10, hd1_11, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_z[k] * hf_22[k];

        t_35[k] = f_1 * hd0_11[k]
                  - f_2 * hd1_11[k]
                  + pb_z[k] * hf_23[k];

        t_36[k] = f_8 * fg0_2[k]
                  - f_9 * fg1_2[k]
                  + pa_z[k] * gg_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, gf_20, gf_21, hd0_12, hd0_14, \
                         hd1_12, hd1_14, hf_24, hf_25, hf_26, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * hf_24[k];

        t_38[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_y[k] * hf_25[k];

        t_39[k] = f_10 * gf_20[k]
                  + f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_x[k] * hf_26[k];

        t_40[k] = f_10 * gf_21[k]
                  + pb_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, fg0_8, fg1_8, gg_30, hd0_13, \
                         hd0_14, hd1_13, hd1_14, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_y[k] * hf_27[k];

        t_42[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_28[k];

        t_43[k] = pb_y[k] * hf_29[k];

        t_44[k] = f_5 * fg0_8[k]
                  - f_6 * fg1_8[k]
                  + pa_x[k] * gg_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, hd0_15, hd0_16, hd0_17, hd1_15, hd1_16, \
                         hd1_17, hf_30, hf_31, hf_32, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_15[k]
                  + pb_x[k] * hf_30[k];

        t_46[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_31[k];

        t_47[k] = f_3 * hd0_17[k]
                  - f_4 * hd1_17[k]
                  + pb_x[k] * hf_32[k];

        t_48[k] = pb_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pb_x, pb_y, pb_z, gf_25, hd0_16, \
                         hd0_17, hd1_16, hd1_17, hf_33, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_x[k] * hf_35[k];

        t_50[k] = f_0 * gf_25[k]
                  + f_1 * hd0_16[k]
                  - f_2 * hd1_16[k]
                  + pb_y[k] * hf_33[k];

        t_51[k] = pb_z[k] * hf_33[k];

        t_52[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_z[k] * hf_34[k];

        t_53[k] = f_1 * hd0_17[k]
                  - f_2 * hd1_17[k]
                  + pb_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, hd0_18, hd0_19, hd0_20, hd1_18, hd1_19, \
                         hd1_20, hf_36, hf_37, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * hd0_18[k]
                  - f_2 * hd1_18[k]
                  + pb_x[k] * hf_36[k];

        t_55[k] = f_3 * hd0_19[k]
                  - f_4 * hd1_19[k]
                  + pb_x[k] * hf_37[k];

        t_56[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_x[k] * hf_38[k];

        t_57[k] = pb_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_x, pb_y, fg0_5, fg1_5, gf_32, gf_33, \
                         gg_40, hd0_20, hd1_20, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_x[k] * hf_41[k];

        t_59[k] = f_5 * fg0_5[k]
                  - f_6 * fg1_5[k]
                  + pa_z[k] * gg_40[k];

        t_60[k] = f_7 * gf_32[k]
                  + f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_y[k] * hf_40[k];

        t_61[k] = f_7 * gf_33[k]
                  + pb_y[k] * hf_41[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pb_x, fg0_7, fg1_7, gg_49, hd0_21, hd0_22, \
                         hd1_21, hd1_22, hf_42, hf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * fg0_7[k]
                  - f_9 * fg1_7[k]
                  + pa_y[k] * gg_49[k];

        t_63[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_21[k]
                  + pb_x[k] * hf_42[k];

        t_64[k] = f_3 * hd0_22[k]
                  - f_4 * hd1_22[k]
                  + pb_x[k] * hf_43[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_x, fg0_6, fg1_6, gg_46, hd0_23, \
                         hd1_23, hf_44, hf_45, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * hd0_23[k]
                  - f_4 * hd1_23[k]
                  + pb_x[k] * hf_44[k];

        t_66[k] = pb_x[k] * hf_45[k];

        t_67[k] = pb_x[k] * hf_47[k];

        t_68[k] = f_8 * fg0_6[k]
                  - f_9 * fg1_6[k]
                  + pa_z[k] * gg_46[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pb_y, fg0_8, fg1_8, gf_34, gf_35, gg_50, \
                         hd0_23, hd1_23, hf_46, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * gf_34[k]
                  + f_3 * hd0_23[k]
                  - f_4 * hd1_23[k]
                  + pb_y[k] * hf_46[k];

        t_70[k] = f_10 * gf_35[k]
                  + pb_y[k] * hf_47[k];

        t_71[k] = f_5 * fg0_8[k]
                  - f_6 * fg1_8[k]
                  + pa_y[k] * gg_50[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, hd0_24, hd0_25, hd0_26, hd1_24, hd1_25, \
                         hd1_26, hf_48, hf_49, hf_50, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * hd0_24[k]
                  - f_2 * hd1_24[k]
                  + pb_x[k] * hf_48[k];

        t_73[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_x[k] * hf_49[k];

        t_74[k] = f_3 * hd0_26[k]
                  - f_4 * hd1_26[k]
                  + pb_x[k] * hf_50[k];

        t_75[k] = pb_x[k] * hf_51[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pb_x, pb_y, pb_z, gf_41, hd0_25, \
                         hd0_26, hd1_25, hd1_26, hf_51, hf_52, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_x[k] * hf_53[k];

        t_77[k] = f_1 * hd0_25[k]
                  - f_2 * hd1_25[k]
                  + pb_y[k] * hf_51[k];

        t_78[k] = f_3 * hd0_26[k]
                  - f_4 * hd1_26[k]
                  + pb_y[k] * hf_52[k];

        t_79[k] = pb_y[k] * hf_53[k];

        t_80[k] = f_0 * gf_41[k]
                  + f_1 * hd0_26[k]
                  - f_2 * hd1_26[k]
                  + pb_z[k] * hf_53[k];
    }
}

auto
compute_prim_hg_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_3 = buffer.data(fg0 + 3);
    const auto *fg0_4 = buffer.data(fg0 + 4);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_6 = buffer.data(fg0 + 6);
    const auto *fg0_7 = buffer.data(fg0 + 7);
    const auto *fg0_8 = buffer.data(fg0 + 8);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_10 = buffer.data(fg1 + 10);
    const auto *fg1_13 = buffer.data(fg1 + 13);
    const auto *fg1_21 = buffer.data(fg1 + 21);
    const auto *fg1_25 = buffer.data(fg1 + 25);
    const auto *fg1_32 = buffer.data(fg1 + 32);
    const auto *fg1_39 = buffer.data(fg1 + 39);
    const auto *fg1_46 = buffer.data(fg1 + 46);
    const auto *fg1_56 = buffer.data(fg1 + 56);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_44 = buffer.data(gf + 44);

    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_78 = buffer.data(gg + 78);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, hd0_1, hd0_2, hd1_1, hd1_2, hf_3, \
                         hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, fg0_0, fg1_0, gf_10, gg_10, hd0_4, \
                         hd1_4, hf_6, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fg0_0[k]
                 - f_6 * fg1_0[k]
                 + pa_y[k] * gg_10[k];

        t_10[k] = pb_z[k] * hf_6[k];

        t_11[k] = f_7 * gf_10[k]
                  + f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_x[k] * hf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, fg0_3, fg1_21, gf_11, \
                         gg_24, hd0_3, hd1_3, hf_7, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hd0_3[k]
                  - f_4 * hd1_3[k]
                  + pb_z[k] * hf_7[k];

        t_13[k] = f_7 * gf_11[k]
                  + pb_x[k] * hf_9[k];

        t_14[k] = f_8 * fg0_3[k]
                  - f_9 * fg1_21[k]
                  + pa_x[k] * gg_24[k];

        t_15[k] = pb_z[k] * hf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, fg0_0, fg1_0, gg_13, hd0_4, hd0_5, \
                         hd1_4, hd1_5, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_z[k] * hf_10[k];

        t_17[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_5[k]
                  + pb_z[k] * hf_11[k];

        t_18[k] = f_5 * fg0_0[k]
                  - f_6 * fg1_0[k]
                  + pa_z[k] * gg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, gf_16, gf_19, hd0_6, hd0_8, \
                         hd1_6, hd1_8, hf_12, hf_13, hf_14, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * hf_12[k];

        t_20[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_y[k] * hf_13[k];

        t_21[k] = f_7 * gf_16[k]
                  + f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_x[k] * hf_14[k];

        t_22[k] = f_7 * gf_19[k]
                  + pb_x[k] * hf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, fg0_4, fg1_25, gg_38, hd0_7, \
                         hd0_8, hd1_7, hd1_8, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_y[k] * hf_15[k];

        t_24[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_16[k];

        t_25[k] = pb_y[k] * hf_17[k];

        t_26[k] = f_8 * fg0_4[k]
                  - f_9 * fg1_25[k]
                  + pa_x[k] * gg_38[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_z, fg0_1, fg1_10, gf_20, gg_19, \
                         hd0_10, hd1_10, hf_18, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_8 * fg0_1[k]
                  - f_9 * fg1_10[k]
                  + pa_y[k] * gg_19[k];

        t_28[k] = pb_z[k] * hf_18[k];

        t_29[k] = f_10 * gf_20[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, fg0_5, fg1_32, gf_21, \
                         gg_42, hd0_9, hd1_9, hf_19, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hd0_9[k]
                  - f_4 * hd1_9[k]
                  + pb_z[k] * hf_19[k];

        t_31[k] = f_10 * gf_21[k]
                  + pb_x[k] * hf_21[k];

        t_32[k] = f_5 * fg0_5[k]
                  - f_6 * fg1_32[k]
                  + pa_x[k] * gg_42[k];

        t_33[k] = pb_z[k] * hf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, fg0_2, fg1_13, gg_30, hd0_10, hd0_11, \
                         hd1_10, hd1_11, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_z[k] * hf_22[k];

        t_35[k] = f_1 * hd0_11[k]
                  - f_2 * hd1_11[k]
                  + pb_z[k] * hf_23[k];

        t_36[k] = f_8 * fg0_2[k]
                  - f_9 * fg1_13[k]
                  + pa_z[k] * gg_30[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, gf_22, gf_23, hd0_12, hd0_14, \
                         hd1_12, hd1_14, hf_24, hf_25, hf_26, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * hf_24[k];

        t_38[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_y[k] * hf_25[k];

        t_39[k] = f_10 * gf_22[k]
                  + f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_x[k] * hf_26[k];

        t_40[k] = f_10 * gf_23[k]
                  + pb_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_y, fg0_8, fg1_56, gg_47, hd0_13, \
                         hd0_14, hd1_13, hd1_14, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_y[k] * hf_27[k];

        t_42[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_28[k];

        t_43[k] = pb_y[k] * hf_29[k];

        t_44[k] = f_5 * fg0_8[k]
                  - f_6 * fg1_56[k]
                  + pa_x[k] * gg_47[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, hd0_15, hd0_16, hd0_17, hd1_15, hd1_16, \
                         hd1_17, hf_30, hf_31, hf_32, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_15[k]
                  + pb_x[k] * hf_30[k];

        t_46[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_31[k];

        t_47[k] = f_3 * hd0_17[k]
                  - f_4 * hd1_17[k]
                  + pb_x[k] * hf_32[k];

        t_48[k] = pb_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pb_x, pb_y, pb_z, gf_27, hd0_16, \
                         hd0_17, hd1_16, hd1_17, hf_33, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_x[k] * hf_35[k];

        t_50[k] = f_0 * gf_27[k]
                  + f_1 * hd0_16[k]
                  - f_2 * hd1_16[k]
                  + pb_y[k] * hf_33[k];

        t_51[k] = pb_z[k] * hf_33[k];

        t_52[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_z[k] * hf_34[k];

        t_53[k] = f_1 * hd0_17[k]
                  - f_2 * hd1_17[k]
                  + pb_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, hd0_18, hd0_19, hd0_20, hd1_18, hd1_19, \
                         hd1_20, hf_36, hf_37, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * hd0_18[k]
                  - f_2 * hd1_18[k]
                  + pb_x[k] * hf_36[k];

        t_55[k] = f_3 * hd0_19[k]
                  - f_4 * hd1_19[k]
                  + pb_x[k] * hf_37[k];

        t_56[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_x[k] * hf_38[k];

        t_57[k] = pb_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_x, pb_y, fg0_5, fg1_32, gf_35, \
                         gf_36, gg_61, hd0_20, hd1_20, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_x[k] * hf_41[k];

        t_59[k] = f_5 * fg0_5[k]
                  - f_6 * fg1_32[k]
                  + pa_z[k] * gg_61[k];

        t_60[k] = f_7 * gf_35[k]
                  + f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_y[k] * hf_40[k];

        t_61[k] = f_7 * gf_36[k]
                  + pb_y[k] * hf_41[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pb_x, fg0_7, fg1_46, gg_72, hd0_21, hd0_22, \
                         hd1_21, hd1_22, hf_42, hf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * fg0_7[k]
                  - f_9 * fg1_46[k]
                  + pa_y[k] * gg_72[k];

        t_63[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_21[k]
                  + pb_x[k] * hf_42[k];

        t_64[k] = f_3 * hd0_22[k]
                  - f_4 * hd1_22[k]
                  + pb_x[k] * hf_43[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_x, fg0_6, fg1_39, gg_69, hd0_23, \
                         hd1_23, hf_44, hf_45, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * hd0_23[k]
                  - f_4 * hd1_23[k]
                  + pb_x[k] * hf_44[k];

        t_66[k] = pb_x[k] * hf_45[k];

        t_67[k] = pb_x[k] * hf_47[k];

        t_68[k] = f_8 * fg0_6[k]
                  - f_9 * fg1_39[k]
                  + pa_z[k] * gg_69[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pb_y, fg0_8, fg1_56, gf_37, gf_38, gg_78, \
                         hd0_23, hd1_23, hf_46, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * gf_37[k]
                  + f_3 * hd0_23[k]
                  - f_4 * hd1_23[k]
                  + pb_y[k] * hf_46[k];

        t_70[k] = f_10 * gf_38[k]
                  + pb_y[k] * hf_47[k];

        t_71[k] = f_5 * fg0_8[k]
                  - f_6 * fg1_56[k]
                  + pa_y[k] * gg_78[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, hd0_24, hd0_25, hd0_26, hd1_24, hd1_25, \
                         hd1_26, hf_48, hf_49, hf_50, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * hd0_24[k]
                  - f_2 * hd1_24[k]
                  + pb_x[k] * hf_48[k];

        t_73[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_x[k] * hf_49[k];

        t_74[k] = f_3 * hd0_26[k]
                  - f_4 * hd1_26[k]
                  + pb_x[k] * hf_50[k];

        t_75[k] = pb_x[k] * hf_51[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pb_x, pb_y, pb_z, gf_44, hd0_25, \
                         hd0_26, hd1_25, hd1_26, hf_51, hf_52, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_x[k] * hf_53[k];

        t_77[k] = f_1 * hd0_25[k]
                  - f_2 * hd1_25[k]
                  + pb_y[k] * hf_51[k];

        t_78[k] = f_3 * hd0_26[k]
                  - f_4 * hd1_26[k]
                  + pb_y[k] * hf_52[k];

        t_79[k] = pb_y[k] * hf_53[k];

        t_80[k] = f_0 * gf_44[k]
                  + f_1 * hd0_26[k]
                  - f_2 * hd1_26[k]
                  + pb_z[k] * hf_53[k];
    }
}

auto
compute_prim_hg_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.5 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_10 = buffer.data(fg0 + 10);
    const auto *fg0_13 = buffer.data(fg0 + 13);
    const auto *fg0_21 = buffer.data(fg0 + 21);
    const auto *fg0_25 = buffer.data(fg0 + 25);
    const auto *fg0_32 = buffer.data(fg0 + 32);
    const auto *fg0_39 = buffer.data(fg0 + 39);
    const auto *fg0_46 = buffer.data(fg0 + 46);
    const auto *fg0_56 = buffer.data(fg0 + 56);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_10 = buffer.data(fg1 + 10);
    const auto *fg1_12 = buffer.data(fg1 + 12);
    const auto *fg1_17 = buffer.data(fg1 + 17);
    const auto *fg1_21 = buffer.data(fg1 + 21);
    const auto *fg1_28 = buffer.data(fg1 + 28);
    const auto *fg1_32 = buffer.data(fg1 + 32);
    const auto *fg1_37 = buffer.data(fg1 + 37);
    const auto *fg1_47 = buffer.data(fg1 + 47);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_47 = buffer.data(gf + 47);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_65 = buffer.data(gg + 65);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, gg_0, hd0_1, hd0_2, hd1_1, \
                         hd1_2, hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];

        t_9[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, fg0_0, fg1_0, gf_3, \
                         gf_5, gg_0, gg_5, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gf_3[k]
                  + pa_y[k] * gg_5[k];

        t_11[k] = pa_y[k] * gg_8[k];

        t_12[k] = pa_z[k] * gg_0[k];

        t_13[k] = pa_z[k] * gg_5[k];

        t_14[k] = f_5 * gf_5[k]
                  + pa_z[k] * gg_8[k];

        t_15[k] = f_6 * fg0_0[k]
                  - f_7 * fg1_0[k]
                  + pa_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_z, gf_11, gf_12, hd0_3, hd0_4, \
                         hd1_3, hd1_4, hf_6, hf_7, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * hf_6[k];

        t_17[k] = f_8 * gf_11[k]
                  + f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_x[k] * hf_8[k];

        t_18[k] = f_3 * hd0_3[k]
                  - f_4 * hd1_3[k]
                  + pb_z[k] * hf_7[k];

        t_19[k] = f_8 * gf_12[k]
                  + pb_x[k] * hf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_z, fg0_21, fg1_17, gg_18, hd0_4, \
                         hd0_5, hd1_4, hd1_5, hf_9, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_9 * fg0_21[k]
                  - f_10 * fg1_17[k]
                  + pa_x[k] * gg_18[k];

        t_21[k] = pb_z[k] * hf_9[k];

        t_22[k] = f_3 * hd0_4[k]
                  - f_4 * hd1_4[k]
                  + pb_z[k] * hf_10[k];

        t_23[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_5[k]
                  + pb_z[k] * hf_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pa_z, pb_y, fg0_0, fg1_0, gg_10, gg_11, \
                         gg_12, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_z[k] * gg_10[k];

        t_25[k] = pa_y[k] * gg_12[k];

        t_26[k] = f_6 * fg0_0[k]
                  - f_7 * fg1_0[k]
                  + pa_z[k] * gg_11[k];

        t_27[k] = pb_y[k] * hf_12[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, pb_y, gf_17, gf_20, hd0_6, hd0_8, hd1_6, \
                         hd1_8, hf_13, hf_14, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_y[k] * hf_13[k];

        t_29[k] = f_8 * gf_17[k]
                  + f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_x[k] * hf_14[k];

        t_30[k] = f_8 * gf_20[k]
                  + pb_x[k] * hf_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_y, fg0_25, fg1_21, gg_30, hd0_7, \
                         hd0_8, hd1_7, hd1_8, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_y[k] * hf_15[k];

        t_32[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_16[k];

        t_33[k] = pb_y[k] * hf_17[k];

        t_34[k] = f_9 * fg0_25[k]
                  - f_10 * fg1_21[k]
                  + pa_x[k] * gg_30[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_x, pb_z, fg0_10, fg1_10, gf_21, gg_13, \
                         hd0_10, hd1_10, hf_18, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_9 * fg0_10[k]
                  - f_10 * fg1_10[k]
                  + pa_y[k] * gg_13[k];

        t_36[k] = pb_z[k] * hf_18[k];

        t_37[k] = f_11 * gf_21[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_x, pb_z, fg0_32, fg1_28, gf_22, \
                         gg_32, hd0_9, hd1_9, hf_19, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * hd0_9[k]
                  - f_4 * hd1_9[k]
                  + pb_z[k] * hf_19[k];

        t_39[k] = f_11 * gf_22[k]
                  + pb_x[k] * hf_21[k];

        t_40[k] = f_6 * fg0_32[k]
                  - f_7 * fg1_28[k]
                  + pa_x[k] * gg_32[k];

        t_41[k] = pb_z[k] * hf_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_z, gg_13, gg_18, hd0_10, hd0_11, \
                         hd1_10, hd1_11, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_z[k] * hf_22[k];

        t_43[k] = f_1 * hd0_11[k]
                  - f_2 * hd1_11[k]
                  + pb_z[k] * hf_23[k];

        t_44[k] = pa_z[k] * gg_13[k];

        t_45[k] = pa_z[k] * gg_18[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, fg0_13, fg1_12, gf_14, gf_18, \
                         gg_21, gg_22, gg_27, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_5 * gf_14[k]
                  + pa_z[k] * gg_21[k];

        t_47[k] = f_5 * gf_18[k]
                  + pa_y[k] * gg_27[k];

        t_48[k] = pa_y[k] * gg_30[k];

        t_49[k] = f_9 * fg0_13[k]
                  - f_10 * fg1_12[k]
                  + pa_z[k] * gg_22[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pb_y, gf_23, gf_24, hd0_12, hd0_14, \
                         hd1_12, hd1_14, hf_24, hf_25, hf_26, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_y[k] * hf_24[k];

        t_51[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_y[k] * hf_25[k];

        t_52[k] = f_11 * gf_23[k]
                  + f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_x[k] * hf_26[k];

        t_53[k] = f_11 * gf_24[k]
                  + pb_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_y, fg0_56, fg1_47, gg_34, hd0_13, \
                         hd0_14, hd1_13, hd1_14, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_y[k] * hf_27[k];

        t_55[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_28[k];

        t_56[k] = pb_y[k] * hf_29[k];

        t_57[k] = f_6 * fg0_56[k]
                  - f_7 * fg1_47[k]
                  + pa_x[k] * gg_34[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pa_z, gf_25, gf_33, gf_42, gg_31, \
                         gg_35, gg_40, gg_46, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_5 * gf_25[k]
                  + pa_x[k] * gg_35[k];

        t_59[k] = pa_x[k] * gg_40[k];

        t_60[k] = pa_z[k] * gg_31[k];

        t_61[k] = f_5 * gf_33[k]
                  + pa_x[k] * gg_46[k];

        t_62[k] = f_5 * gf_42[k]
                  + pa_x[k] * gg_57[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_x, pb_x, gg_65, hd0_15, hd0_16, hd0_17, \
                         hd1_15, hd1_16, hd1_17, hf_30, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_x[k] * gg_65[k];

        t_64[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_15[k]
                  + pb_x[k] * hf_30[k];

        t_65[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_31[k];

        t_66[k] = f_3 * hd0_17[k]
                  - f_4 * hd1_17[k]
                  + pb_x[k] * hf_32[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, gf_28, hd0_16, \
                         hd1_16, hf_33, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * hf_33[k];

        t_68[k] = pb_x[k] * hf_35[k];

        t_69[k] = f_0 * gf_28[k]
                  + f_1 * hd0_16[k]
                  - f_2 * hd1_16[k]
                  + pb_y[k] * hf_33[k];

        t_70[k] = pb_z[k] * hf_33[k];

        t_71[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_z[k] * hf_34[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_z, gf_30, gg_35, gg_40, gg_43, \
                         hd0_17, hd1_17, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * hd0_17[k]
                  - f_2 * hd1_17[k]
                  + pb_z[k] * hf_35[k];

        t_73[k] = pa_z[k] * gg_35[k];

        t_74[k] = pa_z[k] * gg_40[k];

        t_75[k] = f_5 * gf_30[k]
                  + pa_z[k] * gg_43[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, hd0_18, hd0_19, hd0_20, hd1_18, hd1_19, \
                         hd1_20, hf_36, hf_37, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_1 * hd0_18[k]
                  - f_2 * hd1_18[k]
                  + pb_x[k] * hf_36[k];

        t_77[k] = f_3 * hd0_19[k]
                  - f_4 * hd1_19[k]
                  + pb_x[k] * hf_37[k];

        t_78[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_x[k] * hf_38[k];

        t_79[k] = pb_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_z, pb_x, pb_y, fg0_32, fg1_28, gf_37, \
                         gf_38, gg_44, hd0_20, hd1_20, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * hf_41[k];

        t_81[k] = f_6 * fg0_32[k]
                  - f_7 * fg1_28[k]
                  + pa_z[k] * gg_44[k];

        t_82[k] = f_8 * gf_37[k]
                  + f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_y[k] * hf_40[k];

        t_83[k] = f_8 * gf_38[k]
                  + pb_y[k] * hf_41[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_x, fg0_46, fg1_37, gg_54, hd0_21, hd0_22, \
                         hd1_21, hd1_22, hf_42, hf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_9 * fg0_46[k]
                  - f_10 * fg1_37[k]
                  + pa_y[k] * gg_54[k];

        t_85[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_21[k]
                  + pb_x[k] * hf_42[k];

        t_86[k] = f_3 * hd0_22[k]
                  - f_4 * hd1_22[k]
                  + pb_x[k] * hf_43[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_z, pb_x, fg0_39, fg1_32, gg_51, hd0_23, \
                         hd1_23, hf_44, hf_45, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * hd0_23[k]
                  - f_4 * hd1_23[k]
                  + pb_x[k] * hf_44[k];

        t_88[k] = pb_x[k] * hf_45[k];

        t_89[k] = pb_x[k] * hf_47[k];

        t_90[k] = f_9 * fg0_39[k]
                  - f_10 * fg1_32[k]
                  + pa_z[k] * gg_51[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_y, pb_y, fg0_56, fg1_47, gf_40, gf_41, gg_56, \
                         hd0_23, hd1_23, hf_46, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_11 * gf_40[k]
                  + f_3 * hd0_23[k]
                  - f_4 * hd1_23[k]
                  + pb_y[k] * hf_46[k];

        t_92[k] = f_11 * gf_41[k]
                  + pb_y[k] * hf_47[k];

        t_93[k] = f_6 * fg0_56[k]
                  - f_7 * fg1_47[k]
                  + pa_y[k] * gg_56[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pa_y, pb_x, gf_45, gg_62, gg_65, hd0_24, \
                         hd0_25, hd1_24, hd1_25, hf_48, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * gf_45[k]
                  + pa_y[k] * gg_62[k];

        t_95[k] = pa_y[k] * gg_65[k];

        t_96[k] = f_1 * hd0_24[k]
                  - f_2 * hd1_24[k]
                  + pb_x[k] * hf_48[k];

        t_97[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_x[k] * hf_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, pb_x, pb_y, hd0_25, hd0_26, \
                         hd1_25, hd1_26, hf_50, hf_51, hf_52, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_3 * hd0_26[k]
                  - f_4 * hd1_26[k]
                  + pb_x[k] * hf_50[k];

        t_99[k] = pb_x[k] * hf_51[k];

        t_100[k] = pb_x[k] * hf_53[k];

        t_101[k] = f_1 * hd0_25[k]
                   - f_2 * hd1_25[k]
                   + pb_y[k] * hf_51[k];

        t_102[k] = f_3 * hd0_26[k]
                   - f_4 * hd1_26[k]
                   + pb_y[k] * hf_52[k];

        t_103[k] = pb_y[k] * hf_53[k];
    }

#pragma omp simd aligned(t_104, pb_z, gf_47, hd0_26, hd1_26, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_0 * gf_47[k]
                   + f_1 * hd0_26[k]
                   - f_2 * hd1_26[k]
                   + pb_z[k] * hf_53[k];
    }
}

auto
compute_prim_hg_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);

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

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_3 = buffer.data(fg0 + 3);
    const auto *fg0_4 = buffer.data(fg0 + 4);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_6 = buffer.data(fg0 + 6);
    const auto *fg0_7 = buffer.data(fg0 + 7);
    const auto *fg0_8 = buffer.data(fg0 + 8);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_1 = buffer.data(fg1 + 1);
    const auto *fg1_2 = buffer.data(fg1 + 2);
    const auto *fg1_3 = buffer.data(fg1 + 3);
    const auto *fg1_4 = buffer.data(fg1 + 4);
    const auto *fg1_5 = buffer.data(fg1 + 5);
    const auto *fg1_6 = buffer.data(fg1 + 6);
    const auto *fg1_7 = buffer.data(fg1 + 7);
    const auto *fg1_8 = buffer.data(fg1 + 8);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_31 = buffer.data(hd0 + 31);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_34 = buffer.data(hd1 + 34);
    const auto *hd1_35 = buffer.data(hd1 + 35);
    const auto *hd1_36 = buffer.data(hd1 + 36);
    const auto *hd1_40 = buffer.data(hd1 + 40);
    const auto *hd1_41 = buffer.data(hd1 + 41);
    const auto *hd1_42 = buffer.data(hd1 + 42);
    const auto *hd1_43 = buffer.data(hd1 + 43);
    const auto *hd1_44 = buffer.data(hd1 + 44);
    const auto *hd1_45 = buffer.data(hd1 + 45);
    const auto *hd1_48 = buffer.data(hd1 + 48);
    const auto *hd1_49 = buffer.data(hd1 + 49);
    const auto *hd1_50 = buffer.data(hd1 + 50);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gf_0, gf_3, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2, hf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_2[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];

        t_3[k] = f_0 * gf_3[k]
                 + pb_x[k] * hf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, gf_5, hd0_1, hd0_2, hd1_1, \
                         hd1_2, hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * gf_5[k]
                 + pb_x[k] * hf_5[k];

        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pb_x, pb_y, gf_0, gf_1, gf_7, gg_0, gg_1, \
                         hf_6, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_y[k] * gg_0[k];

        t_9[k] = f_5 * gf_0[k]
                 + pb_y[k] * hf_6[k];

        t_10[k] = f_6 * gf_1[k]
                  + pa_y[k] * gg_1[k];

        t_11[k] = f_7 * gf_7[k]
                  + pb_x[k] * hf_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, pb_z, gf_0, gf_2, gf_3, gg_0, \
                         gg_2, gg_3, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gf_3[k]
                  + pa_y[k] * gg_3[k];

        t_13[k] = pa_z[k] * gg_0[k];

        t_14[k] = f_5 * gf_0[k]
                  + pb_z[k] * hf_8[k];

        t_15[k] = f_6 * gf_2[k]
                  + pa_z[k] * gg_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_x, fg0_0, fg1_0, gf_4, gf_5, \
                         gf_10, gg_4, gg_5, gg_6, hf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gf_10[k]
                  + pb_x[k] * hf_10[k];

        t_17[k] = f_6 * gf_4[k]
                  + pa_z[k] * gg_4[k];

        t_18[k] = f_7 * gf_5[k]
                  + pa_z[k] * gg_5[k];

        t_19[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, pb_z, gf_6, gf_12, hd0_5, hd0_6, hd1_9, \
                         hd1_10, hf_11, hf_12, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * gf_6[k]
                  + pb_y[k] * hf_11[k];

        t_21[k] = f_10 * gf_12[k]
                  + f_3 * hd0_6[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_13[k];

        t_22[k] = f_3 * hd0_5[k]
                  - f_4 * hd1_9[k]
                  + pb_z[k] * hf_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_x, pb_z, fg0_3, fg1_3, gf_13, gg_10, \
                         hd0_6, hd1_10, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_10 * gf_13[k]
                  + pb_x[k] * hf_14[k];

        t_24[k] = f_11 * fg0_3[k]
                  - f_12 * fg1_3[k]
                  + pa_x[k] * gg_10[k];

        t_25[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_10[k]
                  + pb_z[k] * hf_15[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_z, fg0_0, fg1_0, gf_8, gg_7, hd0_7, \
                         hd1_11, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_11[k]
                  + pb_z[k] * hf_16[k];

        t_27[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_z[k] * gg_7[k];

        t_28[k] = f_6 * gf_8[k]
                  + pb_z[k] * hf_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, gf_17, gf_19, hd0_8, hd0_10, hd1_14, \
                         hd1_16, hf_18, hf_20, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_18[k];

        t_30[k] = f_10 * gf_17[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_20[k];

        t_31[k] = f_10 * gf_19[k]
                  + pb_x[k] * hf_23[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pb_y, fg0_4, fg1_4, gg_13, hd0_9, hd0_10, \
                         hd1_15, hd1_16, hf_21, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * hd0_9[k]
                  - f_2 * hd1_15[k]
                  + pb_y[k] * hf_21[k];

        t_33[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_16[k]
                  + pb_y[k] * hf_22[k];

        t_34[k] = f_11 * fg0_4[k]
                  - f_12 * fg1_4[k]
                  + pa_x[k] * gg_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_x, pb_y, fg0_1, fg1_1, gf_11, gf_21, gg_8, \
                         hd0_12, hd1_18, hf_24, hf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_11 * fg0_1[k]
                  - f_12 * fg1_1[k]
                  + pa_y[k] * gg_8[k];

        t_36[k] = f_10 * gf_11[k]
                  + pb_y[k] * hf_24[k];

        t_37[k] = f_6 * gf_21[k]
                  + f_3 * hd0_12[k]
                  - f_4 * hd1_18[k]
                  + pb_x[k] * hf_26[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_x, pb_z, fg0_5, fg1_5, gf_22, gg_14, \
                         hd0_11, hd1_17, hf_25, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * hd0_11[k]
                  - f_4 * hd1_17[k]
                  + pb_z[k] * hf_25[k];

        t_39[k] = f_6 * gf_22[k]
                  + pb_x[k] * hf_27[k];

        t_40[k] = f_8 * fg0_5[k]
                  - f_9 * fg1_5[k]
                  + pa_x[k] * gg_14[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pa_z, pb_z, gg_9, gg_11, hd0_12, \
                         hd0_13, hd1_18, hd1_19, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_18[k]
                  + pb_z[k] * hf_28[k];

        t_42[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_19[k]
                  + pb_z[k] * hf_29[k];

        t_43[k] = pa_z[k] * gg_9[k];

        t_44[k] = pa_y[k] * gg_11[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pa_z, pb_z, fg0_2, fg1_2, gf_15, gg_11, \
                         gg_12, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * gg_12[k];

        t_46[k] = f_11 * fg0_2[k]
                  - f_12 * fg1_2[k]
                  + pa_z[k] * gg_11[k];

        t_47[k] = f_10 * gf_15[k]
                  + pb_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, gf_25, gf_26, hd0_14, hd0_16, hd1_25, \
                         hd1_27, hf_32, hf_34, hf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_25[k]
                  + pb_y[k] * hf_32[k];

        t_49[k] = f_6 * gf_25[k]
                  + f_3 * hd0_16[k]
                  - f_4 * hd1_27[k]
                  + pb_x[k] * hf_34[k];

        t_50[k] = f_6 * gf_26[k]
                  + pb_x[k] * hf_37[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_x, pb_y, fg0_8, fg1_8, gg_15, hd0_15, hd0_16, \
                         hd1_26, hd1_27, hf_35, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_26[k]
                  + pb_y[k] * hf_35[k];

        t_52[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_27[k]
                  + pb_y[k] * hf_36[k];

        t_53[k] = f_8 * fg0_8[k]
                  - f_9 * fg1_8[k]
                  + pa_x[k] * gg_15[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_x, pb_y, gf_20, gf_27, gf_29, gf_30, \
                         gg_16, gg_17, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * gf_27[k]
                  + pa_x[k] * gg_16[k];

        t_55[k] = f_7 * gf_20[k]
                  + pb_y[k] * hf_38[k];

        t_56[k] = f_6 * gf_29[k]
                  + pa_x[k] * gg_17[k];

        t_57[k] = f_5 * gf_30[k]
                  + pb_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, pa_x, pb_z, gf_23, gf_47, gg_19, \
                         gg_23, gg_24, gg_25, gg_27, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_x[k] * gg_19[k];

        t_59[k] = pa_x[k] * gg_23[k];

        t_60[k] = pa_x[k] * gg_24[k];

        t_61[k] = pa_x[k] * gg_25[k];

        t_62[k] = f_7 * gf_47[k]
                  + pa_x[k] * gg_27[k];

        t_63[k] = f_7 * gf_23[k]
                  + pb_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, gf_50, gf_54, gg_29, gg_32, \
                         hd0_19, hd1_34, hf_44, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_6 * gf_50[k]
                  + pa_x[k] * gg_29[k];

        t_65[k] = f_5 * gf_54[k]
                  + pb_x[k] * hf_44[k];

        t_66[k] = pa_x[k] * gg_32[k];

        t_67[k] = f_1 * hd0_19[k]
                  - f_2 * hd1_34[k]
                  + pb_x[k] * hf_45[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, gf_27, gf_30, hd0_20, hd0_21, \
                         hd1_35, hd1_36, hf_45, hf_46, hf_47, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * gf_27[k]
                  + pb_y[k] * hf_45[k];

        t_69[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_35[k]
                  + pb_x[k] * hf_46[k];

        t_70[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_36[k]
                  + pb_x[k] * hf_47[k];

        t_71[k] = f_0 * gf_30[k]
                  + f_1 * hd0_20[k]
                  - f_2 * hd1_35[k]
                  + pb_y[k] * hf_48[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_y, pb_z, gf_28, gf_33, gg_18, \
                         hd0_20, hd0_21, hd1_35, hd1_36, hf_49, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_35[k]
                  + pb_z[k] * hf_49[k];

        t_73[k] = f_0 * gf_33[k]
                  + pb_y[k] * hf_51[k];

        t_74[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_36[k]
                  + pb_z[k] * hf_51[k];

        t_75[k] = f_6 * gf_28[k]
                  + pa_z[k] * gg_18[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, gf_30, gf_31, gf_37, gg_19, \
                         gg_20, hf_52, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_z[k] * gg_19[k];

        t_77[k] = f_5 * gf_30[k]
                  + pb_z[k] * hf_52[k];

        t_78[k] = f_6 * gf_31[k]
                  + pa_z[k] * gg_20[k];

        t_79[k] = f_7 * gf_37[k]
                  + pb_y[k] * hf_55[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_z, pb_x, gf_33, gg_21, hd0_23, hd0_24, hd1_40, \
                         hd1_41, hf_56, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_7 * gf_33[k]
                  + pa_z[k] * gg_21[k];

        t_81[k] = f_1 * hd0_23[k]
                  - f_2 * hd1_40[k]
                  + pb_x[k] * hf_56[k];

        t_82[k] = f_3 * hd0_24[k]
                  - f_4 * hd1_41[k]
                  + pb_x[k] * hf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_z, pb_x, pb_z, fg0_5, fg1_5, gf_34, gg_22, \
                         hd0_25, hd1_42, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_42[k]
                  + pb_x[k] * hf_58[k];

        t_84[k] = f_8 * fg0_5[k]
                  - f_9 * fg1_5[k]
                  + pa_z[k] * gg_22[k];

        t_85[k] = f_6 * gf_34[k]
                  + pb_z[k] * hf_59[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_y, pb_y, fg0_7, fg1_7, gf_41, gf_42, gg_25, \
                         hd0_25, hd1_42, hf_61, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_10 * gf_41[k]
                  + f_3 * hd0_25[k]
                  - f_4 * hd1_42[k]
                  + pb_y[k] * hf_61[k];

        t_87[k] = f_10 * gf_42[k]
                  + pb_y[k] * hf_62[k];

        t_88[k] = f_11 * fg0_7[k]
                  - f_12 * fg1_7[k]
                  + pa_y[k] * gg_25[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pb_x, hd0_26, hd0_27, hd0_28, hd1_43, hd1_44, \
                         hd1_45, hf_63, hf_64, hf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_1 * hd0_26[k]
                  - f_2 * hd1_43[k]
                  + pb_x[k] * hf_63[k];

        t_90[k] = f_3 * hd0_27[k]
                  - f_4 * hd1_44[k]
                  + pb_x[k] * hf_64[k];

        t_91[k] = f_3 * hd0_28[k]
                  - f_4 * hd1_45[k]
                  + pb_x[k] * hf_65[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_z, pb_y, pb_z, fg0_6, fg1_6, gf_39, gf_45, \
                         gg_23, hd0_28, hd1_45, hf_66, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_11 * fg0_6[k]
                  - f_12 * fg1_6[k]
                  + pa_z[k] * gg_23[k];

        t_93[k] = f_10 * gf_39[k]
                  + pb_z[k] * hf_66[k];

        t_94[k] = f_6 * gf_45[k]
                  + f_3 * hd0_28[k]
                  - f_4 * hd1_45[k]
                  + pb_y[k] * hf_68[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_y, pb_y, fg0_8, fg1_8, gf_46, gf_48, \
                         gf_51, gg_26, gg_28, gg_30, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_6 * gf_46[k]
                  + pb_y[k] * hf_69[k];

        t_96[k] = f_8 * fg0_8[k]
                  - f_9 * fg1_8[k]
                  + pa_y[k] * gg_26[k];

        t_97[k] = f_6 * gf_48[k]
                  + pa_y[k] * gg_28[k];

        t_98[k] = f_7 * gf_51[k]
                  + pa_y[k] * gg_30[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_y, pb_y, pb_z, gf_43, gf_53, gf_54, \
                         gg_31, gg_32, hf_70, hf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_7 * gf_43[k]
                  + pb_z[k] * hf_70[k];

        t_100[k] = f_6 * gf_53[k]
                   + pa_y[k] * gg_31[k];

        t_101[k] = f_5 * gf_54[k]
                   + pb_y[k] * hf_73[k];

        t_102[k] = pa_y[k] * gg_32[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pb_x, pb_z, gf_47, hd0_30, hd0_31, \
                         hd0_32, hd1_48, hd1_49, hd1_50, hf_74, hf_76, \
                         hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * hd0_30[k]
                   - f_2 * hd1_48[k]
                   + pb_x[k] * hf_74[k];

        t_104[k] = f_0 * gf_47[k]
                   + pb_z[k] * hf_74[k];

        t_105[k] = f_3 * hd0_31[k]
                   - f_4 * hd1_49[k]
                   + pb_x[k] * hf_76[k];

        t_106[k] = f_3 * hd0_32[k]
                   - f_4 * hd1_50[k]
                   + pb_x[k] * hf_77[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_y, pb_z, gf_51, gf_54, hd0_31, hd0_32, \
                         hd1_49, hd1_50, hf_78, hf_80, hf_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_1 * hd0_31[k]
                   - f_2 * hd1_49[k]
                   + pb_y[k] * hf_78[k];

        t_108[k] = f_0 * gf_51[k]
                   + pb_z[k] * hf_78[k];

        t_109[k] = f_3 * hd0_32[k]
                   - f_4 * hd1_50[k]
                   + pb_y[k] * hf_80[k];

        t_110[k] = f_0 * gf_54[k]
                   + f_1 * hd0_32[k]
                   - f_2 * hd1_50[k]
                   + pb_z[k] * hf_81[k];
    }
}

auto
compute_prim_hg_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);

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

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_8 = buffer.data(fg0 + 8);
    const auto *fg0_9 = buffer.data(fg0 + 9);
    const auto *fg0_10 = buffer.data(fg0 + 10);
    const auto *fg0_13 = buffer.data(fg0 + 13);
    const auto *fg0_14 = buffer.data(fg0 + 14);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_1 = buffer.data(fg1 + 1);
    const auto *fg1_2 = buffer.data(fg1 + 2);
    const auto *fg1_5 = buffer.data(fg1 + 5);
    const auto *fg1_8 = buffer.data(fg1 + 8);
    const auto *fg1_9 = buffer.data(fg1 + 9);
    const auto *fg1_10 = buffer.data(fg1 + 10);
    const auto *fg1_13 = buffer.data(fg1 + 13);
    const auto *fg1_14 = buffer.data(fg1 + 14);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
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
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
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
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_46 = buffer.data(gg + 46);
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
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_81 = buffer.data(gg + 81);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_29 = buffer.data(hd0 + 29);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_31 = buffer.data(hd1 + 31);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, hd0_1, hd0_2, hd1_1, hd1_2, \
                         hf_3, hf_5, hf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = pb_z[k] * hf_3[k];

        t_7[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_5[k];

        t_8[k] = pb_y[k] * hf_6[k];

        t_9[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, gf_1, gf_3, gf_5, gg_0, gg_3, \
                         gg_4, gg_5, gg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gg_0[k];

        t_11[k] = f_5 * gf_1[k]
                  + pa_y[k] * gg_3[k];

        t_12[k] = pa_y[k] * gg_4[k];

        t_13[k] = f_6 * gf_3[k]
                  + pa_y[k] * gg_5[k];

        t_14[k] = f_5 * gf_5[k]
                  + pa_y[k] * gg_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pa_z, pb_y, pb_z, gf_0, gf_6, \
                         gg_0, gg_3, gg_8, hf_9, hf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * gf_6[k]
                  + pb_y[k] * hf_9[k];

        t_16[k] = pa_y[k] * gg_8[k];

        t_17[k] = pa_z[k] * gg_0[k];

        t_18[k] = f_7 * gf_0[k]
                  + pb_z[k] * hf_10[k];

        t_19[k] = pa_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, gf_2, gf_3, gf_4, \
                         gg_4, gg_5, gg_7, hf_11, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * gf_2[k]
                  + pa_z[k] * gg_4[k];

        t_21[k] = pa_z[k] * gg_5[k];

        t_22[k] = f_7 * gf_3[k]
                  + pb_z[k] * hf_11[k];

        t_23[k] = f_5 * gf_4[k]
                  + pa_z[k] * gg_7[k];

        t_24[k] = pb_y[k] * hf_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_z, fg0_0, fg1_0, gf_6, gg_8, gg_9, \
                         hf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * gf_6[k]
                  + pa_z[k] * gg_8[k];

        t_26[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_y[k] * gg_9[k];

        t_27[k] = pb_z[k] * hf_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, pb_z, gf_16, gf_17, hd0_3, hd0_4, hd1_5, \
                         hd1_6, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_10 * gf_16[k]
                  + f_3 * hd0_4[k]
                  - f_4 * hd1_6[k]
                  + pb_x[k] * hf_16[k];

        t_29[k] = f_3 * hd0_3[k]
                  - f_4 * hd1_5[k]
                  + pb_z[k] * hf_15[k];

        t_30[k] = f_10 * gf_17[k]
                  + pb_x[k] * hf_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_y, pb_z, fg0_5, fg1_5, gf_9, gg_21, \
                         hd0_4, hd1_6, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_11 * fg0_5[k]
                  - f_12 * fg1_5[k]
                  + pa_x[k] * gg_21[k];

        t_32[k] = pb_z[k] * hf_17[k];

        t_33[k] = f_3 * hd0_4[k]
                  - f_4 * hd1_6[k]
                  + pb_z[k] * hf_18[k];

        t_34[k] = f_5 * gf_9[k]
                  + pb_y[k] * hf_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, pb_z, gg_10, gg_11, gg_13, \
                         gg_14, hd0_5, hd1_7, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_7[k]
                  + pb_z[k] * hf_19[k];

        t_36[k] = pa_y[k] * gg_13[k];

        t_37[k] = pa_z[k] * gg_10[k];

        t_38[k] = pa_y[k] * gg_14[k];

        t_39[k] = pa_z[k] * gg_11[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_y, pb_z, gf_8, gf_12, gf_13, gg_15, \
                         gg_16, hf_20, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * gf_8[k]
                  + pb_z[k] * hf_20[k];

        t_41[k] = f_5 * gf_12[k]
                  + pa_y[k] * gg_15[k];

        t_42[k] = f_7 * gf_13[k]
                  + pb_y[k] * hf_21[k];

        t_43[k] = pa_y[k] * gg_16[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_y, pb_z, fg0_0, fg1_0, gf_10, gg_12, \
                         hd0_6, hd1_8, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_z[k] * gg_12[k];

        t_45[k] = pb_y[k] * hf_22[k];

        t_46[k] = f_5 * gf_10[k]
                  + pb_z[k] * hf_22[k];

        t_47[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, gf_24, gf_27, hd0_7, hd0_8, hd1_9, \
                         hd1_10, hf_24, hf_25, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * gf_24[k]
                  + f_3 * hd0_8[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_24[k];

        t_49[k] = f_10 * gf_27[k]
                  + pb_x[k] * hf_27[k];

        t_50[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_9[k]
                  + pb_y[k] * hf_25[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_x, pb_y, pb_z, fg0_8, fg1_8, gf_11, gg_31, \
                         hd0_8, hd1_10, hf_25, hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * gf_11[k]
                  + pb_z[k] * hf_25[k];

        t_52[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_10[k]
                  + pb_y[k] * hf_26[k];

        t_53[k] = pb_y[k] * hf_27[k];

        t_54[k] = f_11 * fg0_8[k]
                  - f_12 * fg1_8[k]
                  + pa_x[k] * gg_31[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pb_x, pb_z, fg0_1, fg1_1, gf_29, gg_17, \
                         hd0_10, hd1_12, hf_28, hf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_11 * fg0_1[k]
                  - f_12 * fg1_1[k]
                  + pa_y[k] * gg_17[k];

        t_56[k] = pb_z[k] * hf_28[k];

        t_57[k] = f_5 * gf_29[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_12[k]
                  + pb_x[k] * hf_30[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pb_x, pb_z, fg0_9, fg1_9, gf_30, gg_35, \
                         hd0_9, hd1_11, hf_29, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * hd0_9[k]
                  - f_4 * hd1_11[k]
                  + pb_z[k] * hf_29[k];

        t_59[k] = f_5 * gf_30[k]
                  + pb_x[k] * hf_31[k];

        t_60[k] = f_8 * fg0_9[k]
                  - f_9 * fg1_9[k]
                  + pa_x[k] * gg_35[k];

        t_61[k] = pb_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pb_y, pb_z, gf_19, gg_17, hd0_10, \
                         hd0_11, hd1_12, hd1_13, hf_32, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_12[k]
                  + pb_z[k] * hf_32[k];

        t_63[k] = f_10 * gf_19[k]
                  + pb_y[k] * hf_33[k];

        t_64[k] = f_1 * hd0_11[k]
                  - f_2 * hd1_13[k]
                  + pb_z[k] * hf_33[k];

        t_65[k] = pa_z[k] * gg_17[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_z, pb_z, gf_14, gf_15, gf_17, gg_18, \
                         gg_19, gg_21, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_7 * gf_14[k]
                  + pb_z[k] * hf_34[k];

        t_67[k] = pa_z[k] * gg_18[k];

        t_68[k] = f_5 * gf_15[k]
                  + pa_z[k] * gg_19[k];

        t_69[k] = pa_z[k] * gg_21[k];

        t_70[k] = f_7 * gf_17[k]
                  + pb_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pa_z, pb_y, gf_18, gf_19, gf_21, \
                         gg_22, gg_23, gg_24, gg_25, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * gf_18[k]
                  + pa_z[k] * gg_22[k];

        t_72[k] = f_5 * gf_21[k]
                  + pb_y[k] * hf_36[k];

        t_73[k] = f_6 * gf_19[k]
                  + pa_z[k] * gg_23[k];

        t_74[k] = pa_y[k] * gg_24[k];

        t_75[k] = pa_y[k] * gg_25[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pa_y, pb_z, gf_20, gf_23, gf_25, gf_26, \
                         gg_26, gg_27, gg_29, gg_30, hf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * gf_23[k]
                  + pa_y[k] * gg_26[k];

        t_77[k] = pa_y[k] * gg_27[k];

        t_78[k] = f_6 * gf_25[k]
                  + pa_y[k] * gg_29[k];

        t_79[k] = f_5 * gf_20[k]
                  + pb_z[k] * hf_37[k];

        t_80[k] = f_5 * gf_26[k]
                  + pa_y[k] * gg_30[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, fg0_2, fg1_2, gf_27, gg_24, \
                         gg_31, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_7 * gf_27[k]
                  + pb_y[k] * hf_38[k];

        t_82[k] = pa_y[k] * gg_31[k];

        t_83[k] = f_11 * fg0_2[k]
                  - f_12 * fg1_2[k]
                  + pa_z[k] * gg_24[k];

        t_84[k] = pb_y[k] * hf_39[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pb_x, pb_y, pb_z, gf_22, gf_33, hd0_12, hd0_14, \
                         hd1_14, hd1_16, hf_39, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_10 * gf_22[k]
                  + pb_z[k] * hf_39[k];

        t_86[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_40[k];

        t_87[k] = f_5 * gf_33[k]
                  + f_3 * hd0_14[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_41[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_x, pb_y, pb_z, gf_25, gf_34, hd0_13, \
                         hd0_14, hd1_15, hd1_16, hf_42, hf_43, hf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_5 * gf_34[k]
                  + pb_x[k] * hf_44[k];

        t_89[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_15[k]
                  + pb_y[k] * hf_42[k];

        t_90[k] = f_10 * gf_25[k]
                  + pb_z[k] * hf_42[k];

        t_91[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_16[k]
                  + pb_y[k] * hf_43[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_y, fg0_14, fg1_14, gf_35, gf_37, \
                         gg_40, gg_41, gg_43, hf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_y[k] * hf_44[k];

        t_93[k] = f_8 * fg0_14[k]
                  - f_9 * fg1_14[k]
                  + pa_x[k] * gg_40[k];

        t_94[k] = f_6 * gf_35[k]
                  + pa_x[k] * gg_41[k];

        t_95[k] = f_5 * gf_37[k]
                  + pa_x[k] * gg_43[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pa_x, pb_x, gf_38, gf_39, \
                         gg_44, gg_46, gg_48, gg_49, gg_50, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * gf_38[k]
                  + pa_x[k] * gg_44[k];

        t_97[k] = f_7 * gf_39[k]
                  + pb_x[k] * hf_47[k];

        t_98[k] = pa_x[k] * gg_46[k];

        t_99[k] = pa_x[k] * gg_48[k];

        t_100[k] = pa_x[k] * gg_49[k];

        t_101[k] = pa_x[k] * gg_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_x, pa_z, pb_z, gf_28, gf_43, \
                         gg_32, gg_33, gg_51, gg_53, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pa_z[k] * gg_32[k];

        t_103[k] = f_7 * gf_28[k]
                   + pb_z[k] * hf_48[k];

        t_104[k] = pa_z[k] * gg_33[k];

        t_105[k] = f_5 * gf_43[k]
                   + pa_x[k] * gg_51[k];

        t_106[k] = pa_x[k] * gg_53[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, pa_x, pb_z, gf_31, gf_46, gg_54, \
                         gg_55, gg_56, gg_57, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * gg_54[k];

        t_108[k] = pa_x[k] * gg_55[k];

        t_109[k] = pa_x[k] * gg_56[k];

        t_110[k] = f_6 * gf_46[k]
                   + pa_x[k] * gg_57[k];

        t_111[k] = f_5 * gf_31[k]
                   + pb_z[k] * hf_49[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, t_117, pa_x, gf_47, gf_48, gg_58, \
                         gg_59, gg_60, gg_61, gg_62, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_5 * gf_47[k]
                   + pa_x[k] * gg_58[k];

        t_113[k] = f_5 * gf_48[k]
                   + pa_x[k] * gg_59[k];

        t_114[k] = pa_x[k] * gg_60[k];

        t_115[k] = pa_x[k] * gg_61[k];

        t_116[k] = pa_x[k] * gg_62[k];

        t_117[k] = pa_x[k] * gg_63[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, t_123, pa_x, pa_y, gf_52, gg_36, \
                         gg_37, gg_38, gg_64, gg_65, gg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_x[k] * gg_64[k];

        t_119[k] = pa_y[k] * gg_36[k];

        t_120[k] = pa_y[k] * gg_37[k];

        t_121[k] = f_5 * gf_52[k]
                   + pa_x[k] * gg_65[k];

        t_122[k] = pa_y[k] * gg_38[k];

        t_123[k] = pa_x[k] * gg_66[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_x, pb_z, gf_32, gf_56, gg_67, \
                         gg_68, gg_69, gg_71, hf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_x[k] * gg_67[k];

        t_125[k] = pa_x[k] * gg_68[k];

        t_126[k] = pa_x[k] * gg_69[k];

        t_127[k] = f_6 * gf_56[k]
                   + pa_x[k] * gg_71[k];

        t_128[k] = f_6 * gf_32[k]
                   + pb_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pa_x, pb_x, gf_58, gf_59, gf_62, \
                         gg_74, gg_75, gg_77, gg_78, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_5 * gf_58[k]
                   + pa_x[k] * gg_74[k];

        t_130[k] = f_5 * gf_59[k]
                   + pa_x[k] * gg_75[k];

        t_131[k] = f_7 * gf_62[k]
                   + pb_x[k] * hf_52[k];

        t_132[k] = pa_x[k] * gg_77[k];

        t_133[k] = pa_x[k] * gg_78[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_x, pb_z, gg_79, gg_81, \
                         hd0_17, hd0_18, hd1_19, hd1_20, hf_53, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_x[k] * gg_79[k];

        t_135[k] = pa_x[k] * gg_81[k];

        t_136[k] = f_1 * hd0_17[k]
                   - f_2 * hd1_19[k]
                   + pb_x[k] * hf_53[k];

        t_137[k] = pb_z[k] * hf_53[k];

        t_138[k] = f_3 * hd0_18[k]
                   - f_4 * hd1_20[k]
                   + pb_x[k] * hf_55[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, pb_z, gf_39, hd0_18, \
                         hd0_19, hd1_20, hd1_21, hf_56, hf_57, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * hd0_19[k]
                   - f_4 * hd1_21[k]
                   + pb_x[k] * hf_56[k];

        t_140[k] = pb_x[k] * hf_57[k];

        t_141[k] = pb_x[k] * hf_59[k];

        t_142[k] = f_0 * gf_39[k]
                   + f_1 * hd0_18[k]
                   - f_2 * hd1_20[k]
                   + pb_y[k] * hf_57[k];

        t_143[k] = pb_z[k] * hf_57[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_z, pb_y, pb_z, gf_41, gg_41, hd0_18, \
                         hd0_19, hd1_20, hd1_21, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_3 * hd0_18[k]
                   - f_4 * hd1_20[k]
                   + pb_z[k] * hf_58[k];

        t_145[k] = f_0 * gf_41[k]
                   + pb_y[k] * hf_59[k];

        t_146[k] = f_1 * hd0_19[k]
                   - f_2 * hd1_21[k]
                   + pb_z[k] * hf_59[k];

        t_147[k] = pa_z[k] * gg_41[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, pa_z, pb_x, pb_z, gf_35, gf_36, \
                         gg_43, gg_44, gg_46, hf_60, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_7 * gf_35[k]
                   + pb_z[k] * hf_60[k];

        t_149[k] = pa_z[k] * gg_43[k];

        t_150[k] = f_5 * gf_36[k]
                   + pa_z[k] * gg_44[k];

        t_151[k] = pb_x[k] * hf_63[k];

        t_152[k] = pa_z[k] * gg_46[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_z, pb_y, pb_z, gf_39, gf_40, gf_41, \
                         gf_45, gg_48, gg_50, hf_62, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_7 * gf_39[k]
                   + pb_z[k] * hf_62[k];

        t_154[k] = f_5 * gf_40[k]
                   + pa_z[k] * gg_48[k];

        t_155[k] = f_6 * gf_45[k]
                   + pb_y[k] * hf_63[k];

        t_156[k] = f_6 * gf_41[k]
                   + pa_z[k] * gg_50[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_x, pb_z, gf_42, hd0_20, hd0_21, \
                         hd0_22, hd1_23, hd1_24, hd1_25, hf_64, hf_65, \
                         hf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_1 * hd0_20[k]
                   - f_2 * hd1_23[k]
                   + pb_x[k] * hf_64[k];

        t_158[k] = f_5 * gf_42[k]
                   + pb_z[k] * hf_64[k];

        t_159[k] = f_3 * hd0_21[k]
                   - f_4 * hd1_24[k]
                   + pb_x[k] * hf_65[k];

        t_160[k] = f_3 * hd0_22[k]
                   - f_4 * hd1_25[k]
                   + pb_x[k] * hf_66[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pa_z, pb_x, pb_z, fg0_9, fg1_9, gf_44, \
                         gg_52, hf_67, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_x[k] * hf_67[k];

        t_162[k] = pb_x[k] * hf_69[k];

        t_163[k] = f_8 * fg0_9[k]
                   - f_9 * fg1_9[k]
                   + pa_z[k] * gg_52[k];

        t_164[k] = f_5 * gf_44[k]
                   + pb_z[k] * hf_67[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_y, pb_y, fg0_13, fg1_13, gf_50, gf_51, gg_64, \
                         hd0_22, hd1_25, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_10 * gf_50[k]
                   + f_3 * hd0_22[k]
                   - f_4 * hd1_25[k]
                   + pb_y[k] * hf_68[k];

        t_166[k] = f_10 * gf_51[k]
                   + pb_y[k] * hf_69[k];

        t_167[k] = f_11 * fg0_13[k]
                   - f_12 * fg1_13[k]
                   + pa_y[k] * gg_64[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pb_x, pb_z, gf_46, hd0_23, hd0_24, \
                         hd0_25, hd1_26, hd1_27, hd1_28, hf_70, hf_71, \
                         hf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_1 * hd0_23[k]
                   - f_2 * hd1_26[k]
                   + pb_x[k] * hf_70[k];

        t_169[k] = f_10 * gf_46[k]
                   + pb_z[k] * hf_70[k];

        t_170[k] = f_3 * hd0_24[k]
                   - f_4 * hd1_27[k]
                   + pb_x[k] * hf_71[k];

        t_171[k] = f_3 * hd0_25[k]
                   - f_4 * hd1_28[k]
                   + pb_x[k] * hf_72[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_z, pb_x, pb_z, fg0_10, fg1_10, gf_49, \
                         gg_60, hf_73, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_x[k] * hf_73[k];

        t_173[k] = pb_x[k] * hf_75[k];

        t_174[k] = f_11 * fg0_10[k]
                   - f_12 * fg1_10[k]
                   + pa_z[k] * gg_60[k];

        t_175[k] = f_10 * gf_49[k]
                   + pb_z[k] * hf_73[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, fg0_14, fg1_14, gf_54, gf_55, \
                         gg_70, gg_71, hd0_25, hd1_28, hf_74, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_5 * gf_54[k]
                   + f_3 * hd0_25[k]
                   - f_4 * hd1_28[k]
                   + pb_y[k] * hf_74[k];

        t_177[k] = f_5 * gf_55[k]
                   + pb_y[k] * hf_75[k];

        t_178[k] = f_8 * fg0_14[k]
                   - f_9 * fg1_14[k]
                   + pa_y[k] * gg_70[k];

        t_179[k] = pa_y[k] * gg_71[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pa_y, pb_x, gf_57, gf_60, gg_73, \
                         gg_74, gg_75, gg_77, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * gg_73[k];

        t_181[k] = f_5 * gf_57[k]
                   + pa_y[k] * gg_74[k];

        t_182[k] = pa_y[k] * gg_75[k];

        t_183[k] = pb_x[k] * hf_77[k];

        t_184[k] = f_6 * gf_60[k]
                   + pa_y[k] * gg_77[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pb_y, pb_z, gf_53, gf_61, gf_62, \
                         gg_79, gg_81, hf_77, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_6 * gf_53[k]
                   + pb_z[k] * hf_77[k];

        t_186[k] = f_5 * gf_61[k]
                   + pa_y[k] * gg_79[k];

        t_187[k] = f_7 * gf_62[k]
                   + pb_y[k] * hf_79[k];

        t_188[k] = pa_y[k] * gg_81[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, pb_y, pb_z, gf_56, hd0_27, hd0_28, \
                         hd1_30, hd1_31, hf_80, hf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_1 * hd0_27[k]
                   - f_2 * hd1_30[k]
                   + pb_x[k] * hf_80[k];

        t_190[k] = pb_y[k] * hf_80[k];

        t_191[k] = f_0 * gf_56[k]
                   + pb_z[k] * hf_80[k];

        t_192[k] = f_3 * hd0_28[k]
                   - f_4 * hd1_31[k]
                   + pb_x[k] * hf_82[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pb_x, pb_y, pb_z, gf_60, hd0_28, \
                         hd0_29, hd1_31, hd1_32, hf_83, hf_84, hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * hd0_29[k]
                   - f_4 * hd1_32[k]
                   + pb_x[k] * hf_83[k];

        t_194[k] = pb_x[k] * hf_84[k];

        t_195[k] = pb_x[k] * hf_86[k];

        t_196[k] = f_1 * hd0_28[k]
                   - f_2 * hd1_31[k]
                   + pb_y[k] * hf_84[k];

        t_197[k] = f_0 * gf_60[k]
                   + pb_z[k] * hf_84[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pb_z, gf_62, hd0_29, hd1_32, hf_85, \
                         hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * hd0_29[k]
                   - f_4 * hd1_32[k]
                   + pb_y[k] * hf_85[k];

        t_199[k] = pb_y[k] * hf_86[k];

        t_200[k] = f_0 * gf_62[k]
                   + f_1 * hd0_29[k]
                   - f_2 * hd1_32[k]
                   + pb_z[k] * hf_86[k];
    }
}

auto
compute_prim_hg_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_8 = buffer.data(fg0 + 8);
    const auto *fg0_9 = buffer.data(fg0 + 9);
    const auto *fg0_10 = buffer.data(fg0 + 10);
    const auto *fg0_13 = buffer.data(fg0 + 13);
    const auto *fg0_14 = buffer.data(fg0 + 14);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_6 = buffer.data(fg1 + 6);
    const auto *fg1_7 = buffer.data(fg1 + 7);
    const auto *fg1_10 = buffer.data(fg1 + 10);
    const auto *fg1_13 = buffer.data(fg1 + 13);
    const auto *fg1_17 = buffer.data(fg1 + 17);
    const auto *fg1_20 = buffer.data(fg1 + 20);
    const auto *fg1_23 = buffer.data(fg1 + 23);
    const auto *fg1_29 = buffer.data(fg1 + 29);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_44 = buffer.data(gg + 44);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_31 = buffer.data(hd0 + 31);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_33 = buffer.data(hd1 + 33);
    const auto *hd1_34 = buffer.data(hd1 + 34);
    const auto *hd1_35 = buffer.data(hd1 + 35);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, gg_0, hd0_1, hd0_2, hd1_1, \
                         hd1_2, hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];

        t_9[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_z, gf_0, gf_1, gf_3, gg_0, \
                         gg_3, gg_5, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gf_1[k]
                  + pa_y[k] * gg_3[k];

        t_11[k] = f_6 * gf_3[k]
                  + pa_y[k] * gg_5[k];

        t_12[k] = pa_z[k] * gg_0[k];

        t_13[k] = f_7 * gf_0[k]
                  + pb_z[k] * hf_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, fg0_0, fg1_0, gf_2, gf_4, gf_6, \
                         gg_4, gg_6, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * gf_2[k]
                  + pa_z[k] * gg_4[k];

        t_15[k] = f_5 * gf_4[k]
                  + pa_z[k] * gg_6[k];

        t_16[k] = f_6 * gf_6[k]
                  + pa_z[k] * gg_8[k];

        t_17[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pb_z, gf_13, gf_14, hd0_5, hd0_6, \
                         hd1_6, hd1_7, hf_10, hf_11, hf_12, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_z[k] * hf_10[k];

        t_19[k] = f_10 * gf_13[k]
                  + f_3 * hd0_6[k]
                  - f_4 * hd1_7[k]
                  + pb_x[k] * hf_12[k];

        t_20[k] = f_3 * hd0_5[k]
                  - f_4 * hd1_6[k]
                  + pb_z[k] * hf_11[k];

        t_21[k] = f_10 * gf_14[k]
                  + pb_x[k] * hf_13[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_z, fg0_5, fg1_10, gg_14, hd0_6, \
                         hd0_7, hd1_7, hd1_8, hf_13, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_11 * fg0_5[k]
                  - f_12 * fg1_10[k]
                  + pa_x[k] * gg_14[k];

        t_23[k] = pb_z[k] * hf_13[k];

        t_24[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_7[k]
                  + pb_z[k] * hf_14[k];

        t_25[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_8[k]
                  + pb_z[k] * hf_15[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_y, pb_z, fg0_0, fg1_0, gf_9, gg_10, \
                         hd0_8, hd1_9, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_z[k] * gg_10[k];

        t_27[k] = pb_y[k] * hf_16[k];

        t_28[k] = f_5 * gf_9[k]
                  + pb_z[k] * hf_16[k];

        t_29[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_9[k]
                  + pb_y[k] * hf_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_x, pb_y, gf_19, gf_22, hd0_9, hd0_10, \
                         hd1_10, hd1_11, hf_18, hf_19, hf_20, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_10 * gf_19[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_11[k]
                  + pb_x[k] * hf_18[k];

        t_31[k] = f_10 * gf_22[k]
                  + pb_x[k] * hf_21[k];

        t_32[k] = f_1 * hd0_9[k]
                  - f_2 * hd1_10[k]
                  + pb_y[k] * hf_19[k];

        t_33[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_11[k]
                  + pb_y[k] * hf_20[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pa_y, pb_y, pb_z, fg0_1, fg0_8, fg1_6, \
                         fg1_13, gg_11, gg_18, hf_21, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * hf_21[k];

        t_35[k] = f_11 * fg0_8[k]
                  - f_12 * fg1_13[k]
                  + pa_x[k] * gg_18[k];

        t_36[k] = f_11 * fg0_1[k]
                  - f_12 * fg1_6[k]
                  + pa_y[k] * gg_11[k];

        t_37[k] = pb_z[k] * hf_22[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, pb_z, gf_24, gf_25, hd0_11, hd0_12, hd1_12, \
                         hd1_13, hf_23, hf_24, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_5 * gf_24[k]
                  + f_3 * hd0_12[k]
                  - f_4 * hd1_13[k]
                  + pb_x[k] * hf_24[k];

        t_39[k] = f_3 * hd0_11[k]
                  - f_4 * hd1_12[k]
                  + pb_z[k] * hf_23[k];

        t_40[k] = f_5 * gf_25[k]
                  + pb_x[k] * hf_25[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_z, fg0_9, fg1_17, gg_19, hd0_12, \
                         hd0_13, hd1_13, hd1_14, hf_25, hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_8 * fg0_9[k]
                  - f_9 * fg1_17[k]
                  + pa_x[k] * gg_19[k];

        t_42[k] = pb_z[k] * hf_25[k];

        t_43[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_13[k]
                  + pb_z[k] * hf_26[k];

        t_44[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_14[k]
                  + pb_z[k] * hf_27[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pa_z, pb_y, fg0_2, fg1_7, gg_12, \
                         gg_15, gg_16, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * gg_12[k];

        t_46[k] = pa_y[k] * gg_15[k];

        t_47[k] = pa_y[k] * gg_16[k];

        t_48[k] = f_11 * fg0_2[k]
                  - f_12 * fg1_7[k]
                  + pa_z[k] * gg_15[k];

        t_49[k] = pb_y[k] * hf_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, pb_y, pb_z, gf_17, gf_27, hd0_14, hd0_16, \
                         hd1_15, hd1_17, hf_29, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_10 * gf_17[k]
                  + pb_z[k] * hf_29[k];

        t_51[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_15[k]
                  + pb_y[k] * hf_30[k];

        t_52[k] = f_5 * gf_27[k]
                  + f_3 * hd0_16[k]
                  - f_4 * hd1_17[k]
                  + pb_x[k] * hf_31[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, gf_28, hd0_15, hd0_16, hd1_16, \
                         hd1_17, hf_32, hf_33, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_5 * gf_28[k]
                  + pb_x[k] * hf_34[k];

        t_54[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_16[k]
                  + pb_y[k] * hf_32[k];

        t_55[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_17[k]
                  + pb_y[k] * hf_33[k];

        t_56[k] = pb_y[k] * hf_34[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_x, fg0_14, fg1_29, gf_29, gf_31, \
                         gg_20, gg_21, gg_22, gg_26, gg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * fg0_14[k]
                  - f_9 * fg1_29[k]
                  + pa_x[k] * gg_20[k];

        t_58[k] = f_6 * gf_29[k]
                  + pa_x[k] * gg_21[k];

        t_59[k] = f_5 * gf_31[k]
                  + pa_x[k] * gg_22[k];

        t_60[k] = pa_x[k] * gg_26[k];

        t_61[k] = pa_x[k] * gg_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_x, pb_z, gf_26, gf_47, gf_50, gg_32, \
                         gg_34, gg_36, gg_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_x[k] * gg_32[k];

        t_63[k] = pa_x[k] * gg_34[k];

        t_64[k] = f_6 * gf_47[k]
                  + pa_x[k] * gg_36[k];

        t_65[k] = f_6 * gf_26[k]
                  + pb_z[k] * hf_39[k];

        t_66[k] = f_5 * gf_50[k]
                  + pa_x[k] * gg_38[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_x, pb_x, gg_44, hd0_19, hd0_20, hd0_21, \
                         hd1_20, hd1_21, hd1_22, hf_41, hf_42, hf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pa_x[k] * gg_44[k];

        t_68[k] = f_1 * hd0_19[k]
                  - f_2 * hd1_20[k]
                  + pb_x[k] * hf_41[k];

        t_69[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_21[k]
                  + pb_x[k] * hf_42[k];

        t_70[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_22[k]
                  + pb_x[k] * hf_43[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, pb_x, pb_y, pb_z, gf_33, gf_35, \
                         hd0_20, hd1_21, hf_44, hf_45, hf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_x[k] * hf_44[k];

        t_72[k] = pb_x[k] * hf_46[k];

        t_73[k] = f_0 * gf_33[k]
                  + f_1 * hd0_20[k]
                  - f_2 * hd1_21[k]
                  + pb_y[k] * hf_44[k];

        t_74[k] = pb_z[k] * hf_44[k];

        t_75[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_21[k]
                  + pb_z[k] * hf_45[k];

        t_76[k] = f_0 * gf_35[k]
                  + pb_y[k] * hf_46[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pb_z, gf_30, gf_33, gg_23, gg_26, \
                         hd0_21, hd1_22, hf_46, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_22[k]
                  + pb_z[k] * hf_46[k];

        t_78[k] = f_5 * gf_30[k]
                  + pa_z[k] * gg_23[k];

        t_79[k] = pa_z[k] * gg_26[k];

        t_80[k] = f_7 * gf_33[k]
                  + pb_z[k] * hf_47[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_z, pb_x, pb_y, gf_34, gf_35, gf_37, gg_28, \
                         gg_29, hd0_23, hd1_25, hf_48, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * gf_34[k]
                  + pa_z[k] * gg_28[k];

        t_82[k] = f_6 * gf_37[k]
                  + pb_y[k] * hf_48[k];

        t_83[k] = f_6 * gf_35[k]
                  + pa_z[k] * gg_29[k];

        t_84[k] = f_1 * hd0_23[k]
                  - f_2 * hd1_25[k]
                  + pb_x[k] * hf_49[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, hd0_24, hd0_25, hd1_26, hd1_27, hf_50, \
                         hf_51, hf_52, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_3 * hd0_24[k]
                  - f_4 * hd1_26[k]
                  + pb_x[k] * hf_50[k];

        t_86[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_27[k]
                  + pb_x[k] * hf_51[k];

        t_87[k] = pb_x[k] * hf_52[k];

        t_88[k] = pb_x[k] * hf_54[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_z, pb_y, pb_z, fg0_9, fg1_17, gf_36, gf_42, \
                         gg_30, hd0_25, hd1_27, hf_52, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_8 * fg0_9[k]
                  - f_9 * fg1_17[k]
                  + pa_z[k] * gg_30[k];

        t_90[k] = f_5 * gf_36[k]
                  + pb_z[k] * hf_52[k];

        t_91[k] = f_10 * gf_42[k]
                  + f_3 * hd0_25[k]
                  - f_4 * hd1_27[k]
                  + pb_y[k] * hf_53[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_y, pb_x, pb_y, fg0_13, fg1_23, gf_43, gg_34, \
                         hd0_26, hd1_28, hf_54, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_10 * gf_43[k]
                  + pb_y[k] * hf_54[k];

        t_93[k] = f_11 * fg0_13[k]
                  - f_12 * fg1_23[k]
                  + pa_y[k] * gg_34[k];

        t_94[k] = f_1 * hd0_26[k]
                  - f_2 * hd1_28[k]
                  + pb_x[k] * hf_55[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_x, hd0_27, hd0_28, hd1_29, hd1_30, hf_56, \
                         hf_57, hf_58, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * hd0_27[k]
                  - f_4 * hd1_29[k]
                  + pb_x[k] * hf_56[k];

        t_96[k] = f_3 * hd0_28[k]
                  - f_4 * hd1_30[k]
                  + pb_x[k] * hf_57[k];

        t_97[k] = pb_x[k] * hf_58[k];

        t_98[k] = pb_x[k] * hf_60[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pb_y, pb_z, fg0_10, fg1_20, gf_41, gf_45, \
                         gg_31, hd0_28, hd1_30, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * fg0_10[k]
                  - f_12 * fg1_20[k]
                  + pa_z[k] * gg_31[k];

        t_100[k] = f_10 * gf_41[k]
                   + pb_z[k] * hf_58[k];

        t_101[k] = f_5 * gf_45[k]
                   + f_3 * hd0_28[k]
                   - f_4 * hd1_30[k]
                   + pb_y[k] * hf_59[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_y, pb_y, fg0_14, fg1_29, gf_46, gf_48, \
                         gf_51, gg_35, gg_37, gg_41, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_5 * gf_46[k]
                   + pb_y[k] * hf_60[k];

        t_103[k] = f_8 * fg0_14[k]
                   - f_9 * fg1_29[k]
                   + pa_y[k] * gg_35[k];

        t_104[k] = f_5 * gf_48[k]
                   + pa_y[k] * gg_37[k];

        t_105[k] = f_6 * gf_51[k]
                   + pa_y[k] * gg_41[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_y, pb_y, pb_z, gf_44, gf_52, gf_53, \
                         gg_42, gg_44, hf_61, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_6 * gf_44[k]
                   + pb_z[k] * hf_61[k];

        t_107[k] = f_5 * gf_52[k]
                   + pa_y[k] * gg_42[k];

        t_108[k] = f_7 * gf_53[k]
                   + pb_y[k] * hf_62[k];

        t_109[k] = pa_y[k] * gg_44[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_x, pb_z, gf_47, hd0_30, hd0_31, \
                         hd0_32, hd1_33, hd1_34, hd1_35, hf_63, hf_64, \
                         hf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_1 * hd0_30[k]
                   - f_2 * hd1_33[k]
                   + pb_x[k] * hf_63[k];

        t_111[k] = f_0 * gf_47[k]
                   + pb_z[k] * hf_63[k];

        t_112[k] = f_3 * hd0_31[k]
                   - f_4 * hd1_34[k]
                   + pb_x[k] * hf_64[k];

        t_113[k] = f_3 * hd0_32[k]
                   - f_4 * hd1_35[k]
                   + pb_x[k] * hf_65[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pb_x, pb_y, pb_z, gf_51, hd0_31, \
                         hd0_32, hd1_34, hd1_35, hf_66, hf_67, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_x[k] * hf_66[k];

        t_115[k] = pb_x[k] * hf_68[k];

        t_116[k] = f_1 * hd0_31[k]
                   - f_2 * hd1_34[k]
                   + pb_y[k] * hf_66[k];

        t_117[k] = f_0 * gf_51[k]
                   + pb_z[k] * hf_66[k];

        t_118[k] = f_3 * hd0_32[k]
                   - f_4 * hd1_35[k]
                   + pb_y[k] * hf_67[k];
    }

#pragma omp simd aligned(t_119, t_120, pb_y, pb_z, gf_53, hd0_32, hd1_35, \
                         hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_y[k] * hf_68[k];

        t_120[k] = f_0 * gf_53[k]
                   + f_1 * hd0_32[k]
                   - f_2 * hd1_35[k]
                   + pb_z[k] * hf_68[k];
    }
}

auto
compute_prim_hg_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_3 = buffer.data(fg0 + 3);
    const auto *fg0_4 = buffer.data(fg0 + 4);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_6 = buffer.data(fg0 + 6);
    const auto *fg0_7 = buffer.data(fg0 + 7);
    const auto *fg0_8 = buffer.data(fg0 + 8);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_1 = buffer.data(fg1 + 1);
    const auto *fg1_2 = buffer.data(fg1 + 2);
    const auto *fg1_5 = buffer.data(fg1 + 5);
    const auto *fg1_8 = buffer.data(fg1 + 8);
    const auto *fg1_9 = buffer.data(fg1 + 9);
    const auto *fg1_10 = buffer.data(fg1 + 10);
    const auto *fg1_13 = buffer.data(fg1 + 13);
    const auto *fg1_14 = buffer.data(fg1 + 14);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_44 = buffer.data(gf + 44);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_86 = buffer.data(gg + 86);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_31 = buffer.data(hd0 + 31);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_31 = buffer.data(hd1 + 31);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, gg_0, hd0_1, hd0_2, hd1_1, \
                         hd1_2, hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];

        t_9[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pa_z, pb_z, fg0_0, fg1_0, gg_0, gg_10, \
                         hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * gg_0[k];

        t_11[k] = f_5 * fg0_0[k]
                  - f_6 * fg1_0[k]
                  + pa_y[k] * gg_10[k];

        t_12[k] = pb_z[k] * hf_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_z, gf_10, gf_11, hd0_5, hd0_6, hd1_5, \
                         hd1_6, hf_9, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * gf_10[k]
                  + f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_x[k] * hf_10[k];

        t_14[k] = f_3 * hd0_5[k]
                  - f_4 * hd1_5[k]
                  + pb_z[k] * hf_9[k];

        t_15[k] = f_7 * gf_11[k]
                  + pb_x[k] * hf_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_z, fg0_3, fg1_5, gg_23, hd0_6, \
                         hd0_7, hd1_6, hd1_7, hf_11, hf_12, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * fg0_3[k]
                  - f_9 * fg1_5[k]
                  + pa_x[k] * gg_23[k];

        t_17[k] = pb_z[k] * hf_11[k];

        t_18[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_z[k] * hf_12[k];

        t_19[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, fg0_0, fg1_0, gg_13, hd0_8, hd1_8, \
                         hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * fg0_0[k]
                  - f_6 * fg1_0[k]
                  + pa_z[k] * gg_13[k];

        t_21[k] = pb_y[k] * hf_14[k];

        t_22[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, gf_16, gf_19, hd0_9, hd0_10, \
                         hd1_9, hd1_10, hf_16, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * gf_16[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_16[k];

        t_24[k] = f_7 * gf_19[k]
                  + pb_x[k] * hf_19[k];

        t_25[k] = f_1 * hd0_9[k]
                  - f_2 * hd1_9[k]
                  + pb_y[k] * hf_17[k];

        t_26[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pb_y, pb_z, fg0_1, fg0_4, fg1_1, \
                         fg1_8, gg_18, gg_37, hf_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * hf_19[k];

        t_28[k] = f_8 * fg0_4[k]
                  - f_9 * fg1_8[k]
                  + pa_x[k] * gg_37[k];

        t_29[k] = f_8 * fg0_1[k]
                  - f_9 * fg1_1[k]
                  + pa_y[k] * gg_18[k];

        t_30[k] = pb_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, pb_z, gf_20, gf_21, hd0_11, hd0_12, hd1_11, \
                         hd1_12, hf_21, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * gf_20[k]
                  + f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_x[k] * hf_22[k];

        t_32[k] = f_3 * hd0_11[k]
                  - f_4 * hd1_11[k]
                  + pb_z[k] * hf_21[k];

        t_33[k] = f_10 * gf_21[k]
                  + pb_x[k] * hf_23[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_z, fg0_5, fg1_9, gg_41, hd0_12, \
                         hd0_13, hd1_12, hd1_13, hf_23, hf_24, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * fg0_5[k]
                  - f_6 * fg1_9[k]
                  + pa_x[k] * gg_41[k];

        t_35[k] = pb_z[k] * hf_23[k];

        t_36[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_z[k] * hf_24[k];

        t_37[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_z[k] * hf_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_z, pb_y, fg0_2, fg1_2, gg_29, hd0_14, hd1_14, \
                         hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * fg0_2[k]
                  - f_9 * fg1_2[k]
                  + pa_z[k] * gg_29[k];

        t_39[k] = pb_y[k] * hf_26[k];

        t_40[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_27[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pb_x, pb_y, gf_22, gf_23, hd0_15, hd0_16, \
                         hd1_15, hd1_16, hf_28, hf_29, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * gf_22[k]
                  + f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_28[k];

        t_42[k] = f_10 * gf_23[k]
                  + pb_x[k] * hf_31[k];

        t_43[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_15[k]
                  + pb_y[k] * hf_29[k];

        t_44[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_y[k] * hf_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_x, pb_y, fg0_8, fg1_14, gg_46, \
                         gg_53, gg_68, gg_86, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * hf_31[k];

        t_46[k] = f_5 * fg0_8[k]
                  - f_6 * fg1_14[k]
                  + pa_x[k] * gg_46[k];

        t_47[k] = pa_x[k] * gg_53[k];

        t_48[k] = pa_x[k] * gg_68[k];

        t_49[k] = pa_x[k] * gg_86[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, hd0_19, hd0_20, hd0_21, hd1_19, hd1_20, \
                         hd1_21, hf_36, hf_37, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_1 * hd0_19[k]
                  - f_2 * hd1_19[k]
                  + pb_x[k] * hf_36[k];

        t_51[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_x[k] * hf_37[k];

        t_52[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_21[k]
                  + pb_x[k] * hf_38[k];

        t_53[k] = pb_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pb_x, pb_y, pb_z, gf_27, hd0_20, \
                         hd0_21, hd1_20, hd1_21, hf_39, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_x[k] * hf_41[k];

        t_55[k] = f_0 * gf_27[k]
                  + f_1 * hd0_20[k]
                  - f_2 * hd1_20[k]
                  + pb_y[k] * hf_39[k];

        t_56[k] = pb_z[k] * hf_39[k];

        t_57[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_z[k] * hf_40[k];

        t_58[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_21[k]
                  + pb_z[k] * hf_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, gg_53, hd0_23, hd0_24, hd0_25, \
                         hd1_23, hd1_24, hd1_25, hf_43, hf_44, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_z[k] * gg_53[k];

        t_60[k] = f_1 * hd0_23[k]
                  - f_2 * hd1_23[k]
                  + pb_x[k] * hf_43[k];

        t_61[k] = f_3 * hd0_24[k]
                  - f_4 * hd1_24[k]
                  + pb_x[k] * hf_44[k];

        t_62[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_x[k] * hf_45[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_x, pb_y, fg0_5, fg1_9, gf_35, gg_60, \
                         hd0_25, hd1_25, hf_46, hf_47, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * hf_46[k];

        t_64[k] = pb_x[k] * hf_48[k];

        t_65[k] = f_5 * fg0_5[k]
                  - f_6 * fg1_9[k]
                  + pa_z[k] * gg_60[k];

        t_66[k] = f_7 * gf_35[k]
                  + f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_y[k] * hf_47[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pb_x, pb_y, fg0_7, fg1_13, gf_36, gg_70, \
                         hd0_26, hd1_26, hf_48, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_7 * gf_36[k]
                  + pb_y[k] * hf_48[k];

        t_68[k] = f_8 * fg0_7[k]
                  - f_9 * fg1_13[k]
                  + pa_y[k] * gg_70[k];

        t_69[k] = f_1 * hd0_26[k]
                  - f_2 * hd1_26[k]
                  + pb_x[k] * hf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, hd0_27, hd0_28, hd1_27, hd1_28, hf_50, \
                         hf_51, hf_52, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * hd0_27[k]
                  - f_4 * hd1_27[k]
                  + pb_x[k] * hf_50[k];

        t_71[k] = f_3 * hd0_28[k]
                  - f_4 * hd1_28[k]
                  + pb_x[k] * hf_51[k];

        t_72[k] = pb_x[k] * hf_52[k];

        t_73[k] = pb_x[k] * hf_54[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_z, pb_y, fg0_6, fg1_10, gf_37, gf_38, gg_67, \
                         hd0_28, hd1_28, hf_53, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_8 * fg0_6[k]
                  - f_9 * fg1_10[k]
                  + pa_z[k] * gg_67[k];

        t_75[k] = f_10 * gf_37[k]
                  + f_3 * hd0_28[k]
                  - f_4 * hd1_28[k]
                  + pb_y[k] * hf_53[k];

        t_76[k] = f_10 * gf_38[k]
                  + pb_y[k] * hf_54[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_x, fg0_8, fg1_14, gg_76, gg_86, \
                         hd0_30, hd0_31, hd1_30, hd1_31, hf_57, hf_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * fg0_8[k]
                  - f_6 * fg1_14[k]
                  + pa_y[k] * gg_76[k];

        t_78[k] = pa_y[k] * gg_86[k];

        t_79[k] = f_1 * hd0_30[k]
                  - f_2 * hd1_30[k]
                  + pb_x[k] * hf_57[k];

        t_80[k] = f_3 * hd0_31[k]
                  - f_4 * hd1_31[k]
                  + pb_x[k] * hf_58[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, pb_x, pb_y, hd0_31, hd0_32, \
                         hd1_31, hd1_32, hf_59, hf_60, hf_61, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pb_x[k] * hf_59[k];

        t_82[k] = pb_x[k] * hf_60[k];

        t_83[k] = pb_x[k] * hf_62[k];

        t_84[k] = f_1 * hd0_31[k]
                  - f_2 * hd1_31[k]
                  + pb_y[k] * hf_60[k];

        t_85[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pb_y[k] * hf_61[k];

        t_86[k] = pb_y[k] * hf_62[k];
    }

#pragma omp simd aligned(t_87, pb_z, gf_44, hd0_32, hd1_32, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * gf_44[k]
                  + f_1 * hd0_32[k]
                  - f_2 * hd1_32[k]
                  + pb_z[k] * hf_62[k];
    }
}

auto
compute_prim_hg_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fg0, const size_t fg1,
                                     const size_t gf, const size_t gg, const size_t hd0,
                                     const size_t hd1, const size_t hf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 1.5 / p;
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_1 = buffer.data(fg0 + 1);
    const auto *fg0_2 = buffer.data(fg0 + 2);
    const auto *fg0_5 = buffer.data(fg0 + 5);
    const auto *fg0_8 = buffer.data(fg0 + 8);
    const auto *fg0_9 = buffer.data(fg0 + 9);
    const auto *fg0_10 = buffer.data(fg0 + 10);
    const auto *fg0_13 = buffer.data(fg0 + 13);
    const auto *fg0_14 = buffer.data(fg0 + 14);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_9 = buffer.data(fg1 + 9);
    const auto *fg1_11 = buffer.data(fg1 + 11);
    const auto *fg1_16 = buffer.data(fg1 + 16);
    const auto *fg1_20 = buffer.data(fg1 + 20);
    const auto *fg1_26 = buffer.data(fg1 + 26);
    const auto *fg1_30 = buffer.data(fg1 + 30);
    const auto *fg1_35 = buffer.data(fg1 + 35);
    const auto *fg1_44 = buffer.data(fg1 + 44);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_74 = buffer.data(gg + 74);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_31 = buffer.data(hd0 + 31);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_31 = buffer.data(hd1 + 31);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, hd0_1, hd0_2, hd1_1, hd1_2, \
                         hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = pb_z[k] * hf_3[k];

        t_7[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_8[k] = pb_y[k] * hf_5[k];

        t_9[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, gf_2, gf_3, gg_0, \
                         gg_4, gg_5, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gg_0[k];

        t_11[k] = f_5 * gf_3[k]
                  + pa_y[k] * gg_5[k];

        t_12[k] = pa_y[k] * gg_9[k];

        t_13[k] = pa_z[k] * gg_0[k];

        t_14[k] = f_6 * gf_2[k]
                  + pa_z[k] * gg_4[k];

        t_15[k] = pa_z[k] * gg_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_y, fg0_0, fg1_0, gf_4, gf_6, \
                         gg_7, gg_9, gg_10, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * gf_4[k]
                  + pa_z[k] * gg_7[k];

        t_17[k] = pb_y[k] * hf_8[k];

        t_18[k] = f_5 * gf_6[k]
                  + pa_z[k] * gg_9[k];

        t_19[k] = f_7 * fg0_0[k]
                  - f_8 * fg1_0[k]
                  + pa_y[k] * gg_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_x, pb_z, gf_12, gf_13, hd0_5, hd0_6, \
                         hd1_5, hd1_6, hf_9, hf_10, hf_11, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * hf_9[k];

        t_21[k] = f_9 * gf_12[k]
                  + f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_x[k] * hf_11[k];

        t_22[k] = f_3 * hd0_5[k]
                  - f_4 * hd1_5[k]
                  + pb_z[k] * hf_10[k];

        t_23[k] = f_9 * gf_13[k]
                  + pb_x[k] * hf_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_z, fg0_5, fg1_16, gg_19, hd0_6, \
                         hd0_7, hd1_6, hd1_7, hf_12, hf_13, hf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_10 * fg0_5[k]
                  - f_11 * fg1_16[k]
                  + pa_x[k] * gg_19[k];

        t_25[k] = pb_z[k] * hf_12[k];

        t_26[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_z[k] * hf_13[k];

        t_27[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_z[k] * hf_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, fg0_0, fg1_0, gg_11, gg_12, \
                         gg_13, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_z[k] * gg_11[k];

        t_29[k] = pa_y[k] * gg_13[k];

        t_30[k] = f_7 * fg0_0[k]
                  - f_8 * fg1_0[k]
                  + pa_z[k] * gg_12[k];

        t_31[k] = pb_y[k] * hf_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, gf_18, gf_21, hd0_8, hd0_10, hd1_8, \
                         hd1_10, hf_16, hf_17, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_16[k];

        t_33[k] = f_9 * gf_18[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_17[k];

        t_34[k] = f_9 * gf_21[k]
                  + pb_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, fg0_8, fg1_20, gg_31, hd0_9, \
                         hd0_10, hd1_9, hd1_10, hf_18, hf_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * hd0_9[k]
                  - f_2 * hd1_9[k]
                  + pb_y[k] * hf_18[k];

        t_36[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_y[k] * hf_19[k];

        t_37[k] = pb_y[k] * hf_20[k];

        t_38[k] = f_10 * fg0_8[k]
                  - f_11 * fg1_20[k]
                  + pa_x[k] * gg_31[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_z, fg0_1, fg1_9, gf_22, gg_14, \
                         hd0_12, hd1_12, hf_21, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * fg0_1[k]
                  - f_11 * fg1_9[k]
                  + pa_y[k] * gg_14[k];

        t_40[k] = pb_z[k] * hf_21[k];

        t_41[k] = f_6 * gf_22[k]
                  + f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_x[k] * hf_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_z, fg0_9, fg1_26, gf_23, \
                         gg_35, hd0_11, hd1_11, hf_22, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hd0_11[k]
                  - f_4 * hd1_11[k]
                  + pb_z[k] * hf_22[k];

        t_43[k] = f_6 * gf_23[k]
                  + pb_x[k] * hf_24[k];

        t_44[k] = f_7 * fg0_9[k]
                  - f_8 * fg1_26[k]
                  + pa_x[k] * gg_35[k];

        t_45[k] = pb_z[k] * hf_24[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_z, gg_14, gg_19, hd0_12, hd0_13, \
                         hd1_12, hd1_13, hf_25, hf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_z[k] * hf_25[k];

        t_47[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_z[k] * hf_26[k];

        t_48[k] = pa_z[k] * gg_14[k];

        t_49[k] = pa_z[k] * gg_19[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pa_z, fg0_2, fg1_11, gf_15, gf_19, \
                         gg_22, gg_23, gg_28, gg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_5 * gf_15[k]
                  + pa_z[k] * gg_22[k];

        t_51[k] = f_5 * gf_19[k]
                  + pa_y[k] * gg_28[k];

        t_52[k] = pa_y[k] * gg_31[k];

        t_53[k] = f_10 * fg0_2[k]
                  - f_11 * fg1_11[k]
                  + pa_z[k] * gg_23[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, pb_y, gf_24, gf_25, hd0_14, hd0_16, \
                         hd1_14, hd1_16, hf_27, hf_28, hf_29, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * hf_27[k];

        t_55[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_28[k];

        t_56[k] = f_6 * gf_24[k]
                  + f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_29[k];

        t_57[k] = f_6 * gf_25[k]
                  + pb_x[k] * hf_32[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pb_y, fg0_14, fg1_44, gg_39, hd0_15, \
                         hd0_16, hd1_15, hd1_16, hf_30, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_15[k]
                  + pb_y[k] * hf_30[k];

        t_59[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_y[k] * hf_31[k];

        t_60[k] = pb_y[k] * hf_32[k];

        t_61[k] = f_7 * fg0_14[k]
                  - f_8 * fg1_44[k]
                  + pa_x[k] * gg_39[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_x, pa_z, pb_x, gf_26, gf_28, gf_30, \
                         gg_32, gg_40, gg_42, gg_46, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_5 * gf_26[k]
                  + pa_x[k] * gg_40[k];

        t_63[k] = f_6 * gf_28[k]
                  + pa_x[k] * gg_42[k];

        t_64[k] = f_12 * gf_30[k]
                  + pb_x[k] * hf_34[k];

        t_65[k] = pa_x[k] * gg_46[k];

        t_66[k] = pa_z[k] * gg_32[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pa_x, pb_x, gf_35, gf_44, gf_47, gf_50, \
                         gg_52, gg_58, gg_65, gg_68, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * gf_35[k]
                  + pa_x[k] * gg_52[k];

        t_68[k] = pa_x[k] * gg_58[k];

        t_69[k] = f_5 * gf_44[k]
                  + pa_x[k] * gg_65[k];

        t_70[k] = f_6 * gf_47[k]
                  + pa_x[k] * gg_68[k];

        t_71[k] = f_12 * gf_50[k]
                  + pb_x[k] * hf_36[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_x, pb_x, pb_z, gg_74, hd0_19, hd0_20, \
                         hd1_19, hd1_20, hf_37, hf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_x[k] * gg_74[k];

        t_73[k] = f_1 * hd0_19[k]
                  - f_2 * hd1_19[k]
                  + pb_x[k] * hf_37[k];

        t_74[k] = pb_z[k] * hf_37[k];

        t_75[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_x[k] * hf_38[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pb_x, pb_y, pb_z, gf_30, hd0_20, \
                         hd0_21, hd1_20, hd1_21, hf_39, hf_40, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_21[k]
                  + pb_x[k] * hf_39[k];

        t_77[k] = pb_x[k] * hf_40[k];

        t_78[k] = pb_x[k] * hf_42[k];

        t_79[k] = f_0 * gf_30[k]
                  + f_1 * hd0_20[k]
                  - f_2 * hd1_20[k]
                  + pb_y[k] * hf_40[k];

        t_80[k] = pb_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_z, pb_z, gf_27, gg_40, gg_43, hd0_20, \
                         hd0_21, hd1_20, hd1_21, hf_41, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_z[k] * hf_41[k];

        t_82[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_21[k]
                  + pb_z[k] * hf_42[k];

        t_83[k] = pa_z[k] * gg_40[k];

        t_84[k] = f_6 * gf_27[k]
                  + pa_z[k] * gg_43[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_x, gf_31, gf_32, gg_46, gg_48, \
                         gg_49, hd0_23, hd1_23, hf_44, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_x[k] * hf_44[k];

        t_86[k] = pa_z[k] * gg_46[k];

        t_87[k] = f_6 * gf_31[k]
                  + pa_z[k] * gg_48[k];

        t_88[k] = f_5 * gf_32[k]
                  + pa_z[k] * gg_49[k];

        t_89[k] = f_1 * hd0_23[k]
                  - f_2 * hd1_23[k]
                  + pb_x[k] * hf_45[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, hd0_24, hd0_25, hd1_24, hd1_25, hf_46, \
                         hf_47, hf_48, hf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * hd0_24[k]
                  - f_4 * hd1_24[k]
                  + pb_x[k] * hf_46[k];

        t_91[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_x[k] * hf_47[k];

        t_92[k] = pb_x[k] * hf_48[k];

        t_93[k] = pb_x[k] * hf_50[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_z, pb_y, fg0_9, fg1_26, gf_39, gf_40, gg_50, \
                         hd0_25, hd1_25, hf_49, hf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_7 * fg0_9[k]
                  - f_8 * fg1_26[k]
                  + pa_z[k] * gg_50[k];

        t_95[k] = f_9 * gf_39[k]
                  + f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_y[k] * hf_49[k];

        t_96[k] = f_9 * gf_40[k]
                  + pb_y[k] * hf_50[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pa_y, pb_x, fg0_13, fg1_35, gg_60, hd0_26, hd0_27, \
                         hd1_26, hd1_27, hf_51, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_10 * fg0_13[k]
                  - f_11 * fg1_35[k]
                  + pa_y[k] * gg_60[k];

        t_98[k] = f_1 * hd0_26[k]
                  - f_2 * hd1_26[k]
                  + pb_x[k] * hf_51[k];

        t_99[k] = f_3 * hd0_27[k]
                  - f_4 * hd1_27[k]
                  + pb_x[k] * hf_52[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pb_x, fg0_10, fg1_30, gg_57, \
                         hd0_28, hd1_28, hf_53, hf_54, hf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * hd0_28[k]
                   - f_4 * hd1_28[k]
                   + pb_x[k] * hf_53[k];

        t_101[k] = pb_x[k] * hf_54[k];

        t_102[k] = pb_x[k] * hf_56[k];

        t_103[k] = f_10 * fg0_10[k]
                   - f_11 * fg1_30[k]
                   + pa_z[k] * gg_57[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pb_y, fg0_14, fg1_44, gf_42, gf_43, gg_64, \
                         hd0_28, hd1_28, hf_55, hf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_6 * gf_42[k]
                   + f_3 * hd0_28[k]
                   - f_4 * hd1_28[k]
                   + pb_y[k] * hf_55[k];

        t_105[k] = f_6 * gf_43[k]
                   + pb_y[k] * hf_56[k];

        t_106[k] = f_7 * fg0_14[k]
                   - f_8 * fg1_44[k]
                   + pa_y[k] * gg_64[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pb_x, gf_45, gf_48, gf_49, gg_67, \
                         gg_71, gg_72, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * gf_45[k]
                   + pa_y[k] * gg_67[k];

        t_108[k] = pb_x[k] * hf_57[k];

        t_109[k] = f_5 * gf_48[k]
                   + pa_y[k] * gg_71[k];

        t_110[k] = f_6 * gf_49[k]
                   + pa_y[k] * gg_72[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_y, pb_x, pb_y, gf_50, gg_74, hd0_30, \
                         hd1_30, hf_59, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * gf_50[k]
                   + pb_y[k] * hf_59[k];

        t_112[k] = pa_y[k] * gg_74[k];

        t_113[k] = f_1 * hd0_30[k]
                   - f_2 * hd1_30[k]
                   + pb_x[k] * hf_60[k];

        t_114[k] = pb_y[k] * hf_60[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, hd0_31, hd0_32, \
                         hd1_31, hd1_32, hf_61, hf_62, hf_63, hf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * hd0_31[k]
                   - f_4 * hd1_31[k]
                   + pb_x[k] * hf_61[k];

        t_116[k] = f_3 * hd0_32[k]
                   - f_4 * hd1_32[k]
                   + pb_x[k] * hf_62[k];

        t_117[k] = pb_x[k] * hf_63[k];

        t_118[k] = pb_x[k] * hf_65[k];

        t_119[k] = f_1 * hd0_31[k]
                   - f_2 * hd1_31[k]
                   + pb_y[k] * hf_63[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pb_y, pb_z, gf_50, hd0_32, hd1_32, hf_64, \
                         hf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_3 * hd0_32[k]
                   - f_4 * hd1_32[k]
                   + pb_y[k] * hf_64[k];

        t_121[k] = pb_y[k] * hf_65[k];

        t_122[k] = f_0 * gf_50[k]
                   + f_1 * hd0_32[k]
                   - f_2 * hd1_32[k]
                   + pb_z[k] * hf_65[k];
    }
}

auto
compute_prim_hg_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fg0, const size_t fg1,
                                      const size_t gf, const size_t gg, const size_t hd0,
                                      const size_t hd1, const size_t hf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.5 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg0_0 = buffer.data(fg0 + 0);
    const auto *fg0_9 = buffer.data(fg0 + 9);
    const auto *fg0_11 = buffer.data(fg0 + 11);
    const auto *fg0_16 = buffer.data(fg0 + 16);
    const auto *fg0_20 = buffer.data(fg0 + 20);
    const auto *fg0_26 = buffer.data(fg0 + 26);
    const auto *fg0_30 = buffer.data(fg0 + 30);
    const auto *fg0_35 = buffer.data(fg0 + 35);
    const auto *fg0_44 = buffer.data(fg0 + 44);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_9 = buffer.data(fg1 + 9);
    const auto *fg1_10 = buffer.data(fg1 + 10);
    const auto *fg1_13 = buffer.data(fg1 + 13);
    const auto *fg1_16 = buffer.data(fg1 + 16);
    const auto *fg1_22 = buffer.data(fg1 + 22);
    const auto *fg1_26 = buffer.data(fg1 + 26);
    const auto *fg1_29 = buffer.data(fg1 + 29);
    const auto *fg1_38 = buffer.data(fg1 + 38);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_44 = buffer.data(gf + 44);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_44 = buffer.data(gg + 44);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_31 = buffer.data(hd0 + 31);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_31 = buffer.data(hd1 + 31);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gf_0, hd0_0, hd1_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd0_0[k]
                 - f_2 * hd1_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, gg_0, hd0_1, hd0_2, hd1_1, \
                         hd1_2, hf_3, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hd0_1[k]
                 - f_2 * hd1_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_3 * hd0_2[k]
                 - f_4 * hd1_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd0_2[k]
                 - f_2 * hd1_2[k]
                 + pb_z[k] * hf_5[k];

        t_9[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, fg0_0, fg1_0, gf_3, gf_5, gg_0, \
                         gg_5, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gf_3[k]
                  + pa_y[k] * gg_5[k];

        t_11[k] = pa_z[k] * gg_0[k];

        t_12[k] = f_5 * gf_5[k]
                  + pa_z[k] * gg_8[k];

        t_13[k] = f_6 * fg0_0[k]
                  - f_7 * fg1_0[k]
                  + pa_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_z, gf_10, gf_11, hd0_5, hd0_6, \
                         hd1_5, hd1_6, hf_8, hf_9, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_z[k] * hf_8[k];

        t_15[k] = f_8 * gf_10[k]
                  + f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_x[k] * hf_10[k];

        t_16[k] = f_3 * hd0_5[k]
                  - f_4 * hd1_5[k]
                  + pb_z[k] * hf_9[k];

        t_17[k] = f_8 * gf_11[k]
                  + pb_x[k] * hf_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_z, fg0_16, fg1_13, gg_14, hd0_6, \
                         hd0_7, hd1_6, hd1_7, hf_11, hf_12, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * fg0_16[k]
                  - f_10 * fg1_13[k]
                  + pa_x[k] * gg_14[k];

        t_19[k] = pb_z[k] * hf_11[k];

        t_20[k] = f_3 * hd0_6[k]
                  - f_4 * hd1_6[k]
                  + pb_z[k] * hf_12[k];

        t_21[k] = f_1 * hd0_7[k]
                  - f_2 * hd1_7[k]
                  + pb_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_z, pb_y, fg0_0, fg1_0, gg_10, hd0_8, hd1_8, \
                         hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_6 * fg0_0[k]
                  - f_7 * fg1_0[k]
                  + pa_z[k] * gg_10[k];

        t_23[k] = pb_y[k] * hf_14[k];

        t_24[k] = f_3 * hd0_8[k]
                  - f_4 * hd1_8[k]
                  + pb_y[k] * hf_15[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_x, pb_y, gf_16, gf_19, hd0_9, hd0_10, \
                         hd1_9, hd1_10, hf_16, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * gf_16[k]
                  + f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_x[k] * hf_16[k];

        t_26[k] = f_8 * gf_19[k]
                  + pb_x[k] * hf_19[k];

        t_27[k] = f_1 * hd0_9[k]
                  - f_2 * hd1_9[k]
                  + pb_y[k] * hf_17[k];

        t_28[k] = f_3 * hd0_10[k]
                  - f_4 * hd1_10[k]
                  + pb_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_y, pb_y, pb_z, fg0_9, fg0_20, fg1_9, \
                         fg1_16, gg_11, gg_18, hf_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * hf_19[k];

        t_30[k] = f_9 * fg0_20[k]
                  - f_10 * fg1_16[k]
                  + pa_x[k] * gg_18[k];

        t_31[k] = f_9 * fg0_9[k]
                  - f_10 * fg1_9[k]
                  + pa_y[k] * gg_11[k];

        t_32[k] = pb_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_z, gf_20, gf_21, hd0_11, hd0_12, hd1_11, \
                         hd1_12, hf_21, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_11 * gf_20[k]
                  + f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_x[k] * hf_22[k];

        t_34[k] = f_3 * hd0_11[k]
                  - f_4 * hd1_11[k]
                  + pb_z[k] * hf_21[k];

        t_35[k] = f_11 * gf_21[k]
                  + pb_x[k] * hf_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_z, fg0_26, fg1_22, gg_19, hd0_12, \
                         hd0_13, hd1_12, hd1_13, hf_23, hf_24, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * fg0_26[k]
                  - f_7 * fg1_22[k]
                  + pa_x[k] * gg_19[k];

        t_37[k] = pb_z[k] * hf_23[k];

        t_38[k] = f_3 * hd0_12[k]
                  - f_4 * hd1_12[k]
                  + pb_z[k] * hf_24[k];

        t_39[k] = f_1 * hd0_13[k]
                  - f_2 * hd1_13[k]
                  + pb_z[k] * hf_25[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pb_y, fg0_11, fg1_10, gg_15, \
                         hd0_14, hd1_14, hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * gg_15[k];

        t_41[k] = f_9 * fg0_11[k]
                  - f_10 * fg1_10[k]
                  + pa_z[k] * gg_15[k];

        t_42[k] = pb_y[k] * hf_26[k];

        t_43[k] = f_3 * hd0_14[k]
                  - f_4 * hd1_14[k]
                  + pb_y[k] * hf_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, gf_22, gf_23, hd0_15, hd0_16, \
                         hd1_15, hd1_16, hf_28, hf_29, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * gf_22[k]
                  + f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_x[k] * hf_28[k];

        t_45[k] = f_11 * gf_23[k]
                  + pb_x[k] * hf_31[k];

        t_46[k] = f_1 * hd0_15[k]
                  - f_2 * hd1_15[k]
                  + pb_y[k] * hf_29[k];

        t_47[k] = f_3 * hd0_16[k]
                  - f_4 * hd1_16[k]
                  + pb_y[k] * hf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pb_y, fg0_44, fg1_38, gf_24, \
                         gg_20, gg_21, gg_26, gg_31, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * hf_31[k];

        t_49[k] = f_6 * fg0_44[k]
                  - f_7 * fg1_38[k]
                  + pa_x[k] * gg_20[k];

        t_50[k] = f_5 * gf_24[k]
                  + pa_x[k] * gg_21[k];

        t_51[k] = pa_x[k] * gg_26[k];

        t_52[k] = pa_x[k] * gg_31[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_x, pb_x, gf_39, gg_32, gg_34, gg_36, \
                         gg_44, hd0_19, hd1_19, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * gg_32[k];

        t_54[k] = pa_x[k] * gg_34[k];

        t_55[k] = f_5 * gf_39[k]
                  + pa_x[k] * gg_36[k];

        t_56[k] = pa_x[k] * gg_44[k];

        t_57[k] = f_1 * hd0_19[k]
                  - f_2 * hd1_19[k]
                  + pb_x[k] * hf_34[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, gf_27, hd0_20, hd0_21, \
                         hd1_20, hd1_21, hf_35, hf_36, hf_37, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_x[k] * hf_35[k];

        t_59[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_21[k]
                  + pb_x[k] * hf_36[k];

        t_60[k] = pb_x[k] * hf_37[k];

        t_61[k] = pb_x[k] * hf_39[k];

        t_62[k] = f_0 * gf_27[k]
                  + f_1 * hd0_20[k]
                  - f_2 * hd1_20[k]
                  + pb_y[k] * hf_37[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_z, gg_26, hd0_20, hd0_21, hd1_20, \
                         hd1_21, hf_37, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_z[k] * hf_37[k];

        t_64[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pb_z[k] * hf_38[k];

        t_65[k] = f_1 * hd0_21[k]
                  - f_2 * hd1_21[k]
                  + pb_z[k] * hf_39[k];

        t_66[k] = pa_z[k] * gg_26[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_z, pb_x, gf_29, gg_29, hd0_23, hd0_24, hd1_23, \
                         hd1_24, hf_41, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * gf_29[k]
                  + pa_z[k] * gg_29[k];

        t_68[k] = f_1 * hd0_23[k]
                  - f_2 * hd1_23[k]
                  + pb_x[k] * hf_41[k];

        t_69[k] = f_3 * hd0_24[k]
                  - f_4 * hd1_24[k]
                  + pb_x[k] * hf_42[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pb_x, fg0_26, fg1_22, gg_30, hd0_25, \
                         hd1_25, hf_43, hf_44, hf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_x[k] * hf_43[k];

        t_71[k] = pb_x[k] * hf_44[k];

        t_72[k] = pb_x[k] * hf_46[k];

        t_73[k] = f_6 * fg0_26[k]
                  - f_7 * fg1_22[k]
                  + pa_z[k] * gg_30[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pb_y, fg0_35, fg1_29, gf_35, gf_36, gg_34, \
                         hd0_25, hd1_25, hf_45, hf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_8 * gf_35[k]
                  + f_3 * hd0_25[k]
                  - f_4 * hd1_25[k]
                  + pb_y[k] * hf_45[k];

        t_75[k] = f_8 * gf_36[k]
                  + pb_y[k] * hf_46[k];

        t_76[k] = f_9 * fg0_35[k]
                  - f_10 * fg1_29[k]
                  + pa_y[k] * gg_34[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pb_x, hd0_26, hd0_27, hd0_28, hd1_26, hd1_27, \
                         hd1_28, hf_47, hf_48, hf_49, hf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * hd0_26[k]
                  - f_2 * hd1_26[k]
                  + pb_x[k] * hf_47[k];

        t_78[k] = f_3 * hd0_27[k]
                  - f_4 * hd1_27[k]
                  + pb_x[k] * hf_48[k];

        t_79[k] = f_3 * hd0_28[k]
                  - f_4 * hd1_28[k]
                  + pb_x[k] * hf_49[k];

        t_80[k] = pb_x[k] * hf_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_z, pb_x, pb_y, fg0_30, fg1_26, gf_37, \
                         gf_38, gg_31, hd0_28, hd1_28, hf_51, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pb_x[k] * hf_52[k];

        t_82[k] = f_9 * fg0_30[k]
                  - f_10 * fg1_26[k]
                  + pa_z[k] * gg_31[k];

        t_83[k] = f_11 * gf_37[k]
                  + f_3 * hd0_28[k]
                  - f_4 * hd1_28[k]
                  + pb_y[k] * hf_51[k];

        t_84[k] = f_11 * gf_38[k]
                  + pb_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pb_x, fg0_44, fg1_38, gf_42, gg_35, \
                         gg_41, gg_44, hd0_30, hd1_30, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_6 * fg0_44[k]
                  - f_7 * fg1_38[k]
                  + pa_y[k] * gg_35[k];

        t_86[k] = f_5 * gf_42[k]
                  + pa_y[k] * gg_41[k];

        t_87[k] = pa_y[k] * gg_44[k];

        t_88[k] = f_1 * hd0_30[k]
                  - f_2 * hd1_30[k]
                  + pb_x[k] * hf_54[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pb_x, pb_y, hd0_31, hd0_32, hd1_31, \
                         hd1_32, hf_55, hf_56, hf_57, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_3 * hd0_31[k]
                  - f_4 * hd1_31[k]
                  + pb_x[k] * hf_55[k];

        t_90[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pb_x[k] * hf_56[k];

        t_91[k] = pb_x[k] * hf_57[k];

        t_92[k] = pb_x[k] * hf_59[k];

        t_93[k] = f_1 * hd0_31[k]
                  - f_2 * hd1_31[k]
                  + pb_y[k] * hf_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_y, pb_z, gf_44, hd0_32, hd1_32, hf_58, \
                         hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pb_y[k] * hf_58[k];

        t_95[k] = pb_y[k] * hf_59[k];

        t_96[k] = f_0 * gf_44[k]
                  + f_1 * hd0_32[k]
                  - f_2 * hd1_32[k]
                  + pb_z[k] * hf_59[k];
    }
}

}  // namespace simdt2ceri
