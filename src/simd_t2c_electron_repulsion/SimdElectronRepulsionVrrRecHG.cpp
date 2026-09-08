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
    const auto *fg0_15 = buffer.data(fg0 + 15);
    const auto *fg0_30 = buffer.data(fg0 + 30);
    const auto *fg0_55 = buffer.data(fg0 + 55);
    const auto *fg0_89 = buffer.data(fg0 + 89);
    const auto *fg0_100 = buffer.data(fg0 + 100);
    const auto *fg0_115 = buffer.data(fg0 + 115);
    const auto *fg0_134 = buffer.data(fg0 + 134);
    const auto *fg0_149 = buffer.data(fg0 + 149);

    const auto *fg1_0 = buffer.data(fg1 + 0);
    const auto *fg1_15 = buffer.data(fg1 + 15);
    const auto *fg1_30 = buffer.data(fg1 + 30);
    const auto *fg1_55 = buffer.data(fg1 + 55);
    const auto *fg1_89 = buffer.data(fg1 + 89);
    const auto *fg1_100 = buffer.data(fg1 + 100);
    const auto *fg1_115 = buffer.data(fg1 + 115);
    const auto *fg1_134 = buffer.data(fg1 + 134);
    const auto *fg1_149 = buffer.data(fg1 + 149);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_33 = buffer.data(hd0 + 33);
    const auto *hd0_35 = buffer.data(hd0 + 35);
    const auto *hd0_36 = buffer.data(hd0 + 36);
    const auto *hd0_39 = buffer.data(hd0 + 39);
    const auto *hd0_41 = buffer.data(hd0 + 41);
    const auto *hd0_54 = buffer.data(hd0 + 54);
    const auto *hd0_57 = buffer.data(hd0 + 57);
    const auto *hd0_59 = buffer.data(hd0 + 59);
    const auto *hd0_90 = buffer.data(hd0 + 90);
    const auto *hd0_93 = buffer.data(hd0 + 93);
    const auto *hd0_95 = buffer.data(hd0 + 95);
    const auto *hd0_102 = buffer.data(hd0 + 102);
    const auto *hd0_105 = buffer.data(hd0 + 105);
    const auto *hd0_107 = buffer.data(hd0 + 107);
    const auto *hd0_108 = buffer.data(hd0 + 108);
    const auto *hd0_111 = buffer.data(hd0 + 111);
    const auto *hd0_113 = buffer.data(hd0 + 113);
    const auto *hd0_120 = buffer.data(hd0 + 120);
    const auto *hd0_123 = buffer.data(hd0 + 123);
    const auto *hd0_125 = buffer.data(hd0 + 125);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_33 = buffer.data(hd1 + 33);
    const auto *hd1_35 = buffer.data(hd1 + 35);
    const auto *hd1_36 = buffer.data(hd1 + 36);
    const auto *hd1_39 = buffer.data(hd1 + 39);
    const auto *hd1_41 = buffer.data(hd1 + 41);
    const auto *hd1_54 = buffer.data(hd1 + 54);
    const auto *hd1_57 = buffer.data(hd1 + 57);
    const auto *hd1_59 = buffer.data(hd1 + 59);
    const auto *hd1_90 = buffer.data(hd1 + 90);
    const auto *hd1_93 = buffer.data(hd1 + 93);
    const auto *hd1_95 = buffer.data(hd1 + 95);
    const auto *hd1_102 = buffer.data(hd1 + 102);
    const auto *hd1_105 = buffer.data(hd1 + 105);
    const auto *hd1_107 = buffer.data(hd1 + 107);
    const auto *hd1_108 = buffer.data(hd1 + 108);
    const auto *hd1_111 = buffer.data(hd1 + 111);
    const auto *hd1_113 = buffer.data(hd1 + 113);
    const auto *hd1_120 = buffer.data(hd1 + 120);
    const auto *hd1_123 = buffer.data(hd1 + 123);
    const auto *hd1_125 = buffer.data(hd1 + 125);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

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

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, gf_6, gf_9, hd0_3, hd1_3, \
                         hf_3, hf_5, hf_6, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gf_6[k]
                 + pb_x[k] * hf_6[k];

        t_7[k] = pb_z[k] * hf_3[k];

        t_8[k] = pb_y[k] * hf_5[k];

        t_9[k] = f_0 * gf_9[k]
                 + pb_x[k] * hf_9[k];

        t_10[k] = f_1 * hd0_3[k]
                  - f_2 * hd1_3[k]
                  + pb_y[k] * hf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, gg_0, hd0_5, hd1_5, \
                         hf_6, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * hf_6[k];

        t_12[k] = f_3 * hd0_5[k]
                  - f_4 * hd1_5[k]
                  + pb_y[k] * hf_8[k];

        t_13[k] = pb_y[k] * hf_9[k];

        t_14[k] = f_1 * hd0_5[k]
                  - f_2 * hd1_5[k]
                  + pb_z[k] * hf_9[k];

        t_15[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, gf_0, gf_1, gg_3, \
                         gg_5, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * gf_0[k]
                  + pb_y[k] * hf_10[k];

        t_17[k] = pb_z[k] * hf_10[k];

        t_18[k] = f_6 * gf_1[k]
                  + pa_y[k] * gg_3[k];

        t_19[k] = pb_z[k] * hf_11[k];

        t_20[k] = pa_y[k] * gg_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, gf_6, gf_16, gf_18, \
                         gg_9, gg_10, hf_13, hf_16, hf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * gf_16[k]
                  + pb_x[k] * hf_16[k];

        t_22[k] = pb_z[k] * hf_13[k];

        t_23[k] = f_7 * gf_18[k]
                  + pb_x[k] * hf_18[k];

        t_24[k] = pa_y[k] * gg_9[k];

        t_25[k] = f_7 * gf_6[k]
                  + pa_y[k] * gg_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, gf_8, gf_9, \
                         gg_0, gg_12, gg_14, hf_16, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * hf_16[k];

        t_27[k] = f_6 * gf_8[k]
                  + pa_y[k] * gg_12[k];

        t_28[k] = f_5 * gf_9[k]
                  + pb_y[k] * hf_19[k];

        t_29[k] = pa_y[k] * gg_14[k];

        t_30[k] = pa_z[k] * gg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, gf_0, gf_2, \
                         gg_3, gg_5, gg_6, hf_20, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * hf_20[k];

        t_32[k] = f_5 * gf_0[k]
                  + pb_z[k] * hf_20[k];

        t_33[k] = pa_z[k] * gg_3[k];

        t_34[k] = pb_y[k] * hf_22[k];

        t_35[k] = f_6 * gf_2[k]
                  + pa_z[k] * gg_5[k];

        t_36[k] = pa_z[k] * gg_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, gf_27, gf_29, gg_10, hf_25, \
                         hf_27, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * gf_27[k]
                  + pb_x[k] * hf_27[k];

        t_38[k] = pb_y[k] * hf_25[k];

        t_39[k] = f_7 * gf_29[k]
                  + pb_x[k] * hf_29[k];

        t_40[k] = pa_z[k] * gg_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, gf_6, gf_7, gf_9, gg_12, \
                         gg_14, hf_26, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * gf_6[k]
                  + pb_z[k] * hf_26[k];

        t_42[k] = f_6 * gf_7[k]
                  + pa_z[k] * gg_12[k];

        t_43[k] = pb_y[k] * hf_29[k];

        t_44[k] = f_7 * gf_9[k]
                  + pa_z[k] * gg_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, fg0_0, fg1_0, gf_10, gg_15, \
                         hf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_y[k] * gg_15[k];

        t_46[k] = f_6 * gf_10[k]
                  + pb_y[k] * hf_30[k];

        t_47[k] = pb_z[k] * hf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, gf_33, gf_36, hd0_18, hd0_21, \
                         hd1_18, hd1_21, hf_31, hf_32, hf_33, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * gf_33[k]
                  + f_3 * hd0_21[k]
                  - f_4 * hd1_21[k]
                  + pb_x[k] * hf_33[k];

        t_49[k] = pb_z[k] * hf_31[k];

        t_50[k] = f_3 * hd0_18[k]
                  - f_4 * hd1_18[k]
                  + pb_z[k] * hf_32[k];

        t_51[k] = f_10 * gf_36[k]
                  + pb_x[k] * hf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, fg0_55, fg1_55, gf_38, \
                         gf_39, gg_55, hf_33, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * hf_33[k];

        t_53[k] = f_10 * gf_38[k]
                  + pb_x[k] * hf_38[k];

        t_54[k] = f_10 * gf_39[k]
                  + pb_x[k] * hf_39[k];

        t_55[k] = f_11 * fg0_55[k]
                  - f_12 * fg1_55[k]
                  + pa_x[k] * gg_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, gf_19, hd0_21, hd0_23, hd1_21, \
                         hd1_23, hf_36, hf_37, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * hf_36[k];

        t_57[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_21[k]
                  + pb_z[k] * hf_37[k];

        t_58[k] = f_6 * gf_19[k]
                  + pb_y[k] * hf_39[k];

        t_59[k] = f_1 * hd0_23[k]
                  - f_2 * hd1_23[k]
                  + pb_z[k] * hf_39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, gf_22, gg_16, \
                         gg_18, gg_30, gg_32, gg_35, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * gg_30[k];

        t_61[k] = pa_z[k] * gg_16[k];

        t_62[k] = pa_y[k] * gg_32[k];

        t_63[k] = pa_z[k] * gg_18[k];

        t_64[k] = f_5 * gf_22[k]
                  + pb_y[k] * hf_42[k];

        t_65[k] = pa_y[k] * gg_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, gf_47, gf_48, gg_21, \
                         gg_25, gg_39, hf_47, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * gg_21[k];

        t_67[k] = f_10 * gf_47[k]
                  + pb_x[k] * hf_47[k];

        t_68[k] = f_10 * gf_48[k]
                  + pb_x[k] * hf_48[k];

        t_69[k] = pa_y[k] * gg_39[k];

        t_70[k] = pa_z[k] * gg_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, gf_16, gf_28, gf_29, gg_42, \
                         gg_44, hf_46, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * gf_16[k]
                  + pb_z[k] * hf_46[k];

        t_72[k] = f_6 * gf_28[k]
                  + pa_y[k] * gg_42[k];

        t_73[k] = f_5 * gf_29[k]
                  + pb_y[k] * hf_49[k];

        t_74[k] = pa_y[k] * gg_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, fg0_0, fg1_0, gf_20, gg_30, \
                         hd0_30, hd1_30, hf_50, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * fg0_0[k]
                  - f_9 * fg1_0[k]
                  + pa_z[k] * gg_30[k];

        t_76[k] = pb_y[k] * hf_50[k];

        t_77[k] = f_6 * gf_20[k]
                  + pb_z[k] * hf_50[k];

        t_78[k] = f_3 * hd0_30[k]
                  - f_4 * hd1_30[k]
                  + pb_y[k] * hf_51[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, gf_55, gf_56, gf_57, \
                         hd0_35, hd1_35, hf_52, hf_55, hf_56, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * hf_52[k];

        t_80[k] = f_10 * gf_55[k]
                  + f_3 * hd0_35[k]
                  - f_4 * hd1_35[k]
                  + pb_x[k] * hf_55[k];

        t_81[k] = f_10 * gf_56[k]
                  + pb_x[k] * hf_56[k];

        t_82[k] = f_10 * gf_57[k]
                  + pb_x[k] * hf_57[k];

        t_83[k] = pb_y[k] * hf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, gf_26, gf_59, hd0_33, \
                         hd0_35, hd1_33, hd1_35, hf_56, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_10 * gf_59[k]
                  + pb_x[k] * hf_59[k];

        t_85[k] = f_1 * hd0_33[k]
                  - f_2 * hd1_33[k]
                  + pb_y[k] * hf_56[k];

        t_86[k] = f_6 * gf_26[k]
                  + pb_z[k] * hf_56[k];

        t_87[k] = f_3 * hd0_35[k]
                  - f_4 * hd1_35[k]
                  + pb_y[k] * hf_58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, fg0_15, fg0_89, fg1_15, \
                         fg1_89, gf_30, gg_45, gg_89, hf_59, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * hf_59[k];

        t_89[k] = f_11 * fg0_89[k]
                  - f_12 * fg1_89[k]
                  + pa_x[k] * gg_89[k];

        t_90[k] = f_11 * fg0_15[k]
                  - f_12 * fg1_15[k]
                  + pa_y[k] * gg_45[k];

        t_91[k] = f_10 * gf_30[k]
                  + pb_y[k] * hf_60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, gf_63, hd0_36, hd0_39, hd1_36, \
                         hd1_39, hf_60, hf_61, hf_62, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * hf_60[k];

        t_93[k] = f_6 * gf_63[k]
                  + f_3 * hd0_39[k]
                  - f_4 * hd1_39[k]
                  + pb_x[k] * hf_63[k];

        t_94[k] = pb_z[k] * hf_61[k];

        t_95[k] = f_3 * hd0_36[k]
                  - f_4 * hd1_36[k]
                  + pb_z[k] * hf_62[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, gf_66, gf_68, gf_69, hf_63, \
                         hf_66, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_6 * gf_66[k]
                  + pb_x[k] * hf_66[k];

        t_97[k] = pb_z[k] * hf_63[k];

        t_98[k] = f_6 * gf_68[k]
                  + pb_x[k] * hf_68[k];

        t_99[k] = f_6 * gf_69[k]
                  + pb_x[k] * hf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, fg0_100, fg1_100, \
                         gf_39, gg_100, hd0_39, hd1_39, hf_66, hf_67, \
                         hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_8 * fg0_100[k]
                   - f_9 * fg1_100[k]
                   + pa_x[k] * gg_100[k];

        t_101[k] = pb_z[k] * hf_66[k];

        t_102[k] = f_3 * hd0_39[k]
                   - f_4 * hd1_39[k]
                   + pb_z[k] * hf_67[k];

        t_103[k] = f_10 * gf_39[k]
                   + pb_y[k] * hf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, gf_30, gg_45, gg_46, \
                         gg_48, hd0_41, hd1_41, hf_69, hf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * hd0_41[k]
                   - f_2 * hd1_41[k]
                   + pb_z[k] * hf_69[k];

        t_105[k] = pa_z[k] * gg_45[k];

        t_106[k] = pa_z[k] * gg_46[k];

        t_107[k] = f_5 * gf_30[k]
                   + pb_z[k] * hf_70[k];

        t_108[k] = pa_z[k] * gg_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, gf_32, gf_42, gf_77, \
                         gg_50, gg_51, hf_72, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * gf_42[k]
                   + pb_y[k] * hf_72[k];

        t_110[k] = f_6 * gf_32[k]
                   + pa_z[k] * gg_50[k];

        t_111[k] = pa_z[k] * gg_51[k];

        t_112[k] = f_6 * gf_77[k]
                   + pb_x[k] * hf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, gf_36, gf_78, gf_79, \
                         gg_55, hf_76, hf_78, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_6 * gf_78[k]
                   + pb_x[k] * hf_78[k];

        t_114[k] = f_6 * gf_79[k]
                   + pb_x[k] * hf_79[k];

        t_115[k] = pa_z[k] * gg_55[k];

        t_116[k] = f_5 * gf_36[k]
                   + pb_z[k] * hf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, gf_37, gf_39, gf_49, \
                         gg_57, gg_59, gg_75, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * gf_37[k]
                   + pa_z[k] * gg_57[k];

        t_118[k] = f_6 * gf_49[k]
                   + pb_y[k] * hf_79[k];

        t_119[k] = f_7 * gf_39[k]
                   + pa_z[k] * gg_59[k];

        t_120[k] = pa_y[k] * gg_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, gf_50, gf_51, gf_52, \
                         gg_77, gg_78, gg_80, hf_80, hf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * gf_50[k]
                   + pb_y[k] * hf_80[k];

        t_122[k] = pa_y[k] * gg_77[k];

        t_123[k] = f_6 * gf_51[k]
                   + pa_y[k] * gg_78[k];

        t_124[k] = f_5 * gf_52[k]
                   + pb_y[k] * hf_82[k];

        t_125[k] = pa_y[k] * gg_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, gf_56, gf_86, gf_87, \
                         gf_88, gg_84, gg_85, hf_86, hf_87, hf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_6 * gf_86[k]
                   + pb_x[k] * hf_86[k];

        t_127[k] = f_6 * gf_87[k]
                   + pb_x[k] * hf_87[k];

        t_128[k] = f_6 * gf_88[k]
                   + pb_x[k] * hf_88[k];

        t_129[k] = pa_y[k] * gg_84[k];

        t_130[k] = f_7 * gf_56[k]
                   + pa_y[k] * gg_85[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, gf_46, gf_58, gf_59, \
                         gg_87, gg_89, hf_86, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * gf_46[k]
                   + pb_z[k] * hf_86[k];

        t_132[k] = f_6 * gf_58[k]
                   + pa_y[k] * gg_87[k];

        t_133[k] = f_5 * gf_59[k]
                   + pb_y[k] * hf_89[k];

        t_134[k] = pa_y[k] * gg_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, fg0_30, fg1_30, gf_50, \
                         gg_75, hd0_54, hd1_54, hf_90, hf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * fg0_30[k]
                   - f_12 * fg1_30[k]
                   + pa_z[k] * gg_75[k];

        t_136[k] = pb_y[k] * hf_90[k];

        t_137[k] = f_10 * gf_50[k]
                   + pb_z[k] * hf_90[k];

        t_138[k] = f_3 * hd0_54[k]
                   - f_4 * hd1_54[k]
                   + pb_y[k] * hf_91[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, gf_95, gf_96, gf_97, \
                         hd0_59, hd1_59, hf_92, hf_95, hf_96, hf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * hf_92[k];

        t_140[k] = f_6 * gf_95[k]
                   + f_3 * hd0_59[k]
                   - f_4 * hd1_59[k]
                   + pb_x[k] * hf_95[k];

        t_141[k] = f_6 * gf_96[k]
                   + pb_x[k] * hf_96[k];

        t_142[k] = f_6 * gf_97[k]
                   + pb_x[k] * hf_97[k];

        t_143[k] = pb_y[k] * hf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, gf_56, gf_99, hd0_57, \
                         hd0_59, hd1_57, hd1_59, hf_96, hf_98, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_6 * gf_99[k]
                   + pb_x[k] * hf_99[k];

        t_145[k] = f_1 * hd0_57[k]
                   - f_2 * hd1_57[k]
                   + pb_y[k] * hf_96[k];

        t_146[k] = f_10 * gf_56[k]
                   + pb_z[k] * hf_96[k];

        t_147[k] = f_3 * hd0_59[k]
                   - f_4 * hd1_59[k]
                   + pb_y[k] * hf_98[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, pa_x, pb_y, pb_z, fg0_149, \
                         fg1_149, gf_60, gf_100, gg_149, gg_150, hf_99, \
                         hf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * hf_99[k];

        t_149[k] = f_8 * fg0_149[k]
                   - f_9 * fg1_149[k]
                   + pa_x[k] * gg_149[k];

        t_150[k] = f_7 * gf_100[k]
                   + pa_x[k] * gg_150[k];

        t_151[k] = f_7 * gf_60[k]
                   + pb_y[k] * hf_100[k];

        t_152[k] = pb_z[k] * hf_100[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, gf_103, gf_105, \
                         gf_106, gg_153, gg_155, hf_101, hf_103, \
                         hf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_6 * gf_103[k]
                   + pa_x[k] * gg_153[k];

        t_154[k] = pb_z[k] * hf_101[k];

        t_155[k] = f_6 * gf_105[k]
                   + pa_x[k] * gg_155[k];

        t_156[k] = f_5 * gf_106[k]
                   + pb_x[k] * hf_106[k];

        t_157[k] = pb_z[k] * hf_103[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_x, pb_x, pb_z, gf_108, gf_109, \
                         gg_160, gg_162, hf_106, hf_108, hf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * gf_108[k]
                   + pb_x[k] * hf_108[k];

        t_159[k] = f_5 * gf_109[k]
                   + pb_x[k] * hf_109[k];

        t_160[k] = pa_x[k] * gg_160[k];

        t_161[k] = pb_z[k] * hf_106[k];

        t_162[k] = pa_x[k] * gg_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pa_z, pb_z, gf_60, \
                         gg_90, gg_91, gg_93, gg_163, gg_164, hf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * gg_163[k];

        t_164[k] = pa_x[k] * gg_164[k];

        t_165[k] = pa_z[k] * gg_90[k];

        t_166[k] = pa_z[k] * gg_91[k];

        t_167[k] = f_5 * gf_60[k]
                   + pb_z[k] * hf_110[k];

        t_168[k] = pa_z[k] * gg_93[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_x, pa_z, pb_x, pb_y, gf_72, gf_115, \
                         gf_117, gg_96, gg_170, hf_112, hf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_10 * gf_72[k]
                   + pb_y[k] * hf_112[k];

        t_170[k] = f_6 * gf_115[k]
                   + pa_x[k] * gg_170[k];

        t_171[k] = pa_z[k] * gg_96[k];

        t_172[k] = f_5 * gf_117[k]
                   + pb_x[k] * hf_117[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, pa_x, pb_x, gf_118, gf_119, \
                         gg_175, gg_176, gg_177, gg_178, hf_118, \
                         hf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_5 * gf_118[k]
                   + pb_x[k] * hf_118[k];

        t_174[k] = f_5 * gf_119[k]
                   + pb_x[k] * hf_119[k];

        t_175[k] = pa_x[k] * gg_175[k];

        t_176[k] = pa_x[k] * gg_176[k];

        t_177[k] = pa_x[k] * gg_177[k];

        t_178[k] = pa_x[k] * gg_178[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pa_x, pb_y, pb_z, gf_70, gf_80, \
                         gf_120, gf_123, gg_179, gg_180, gg_183, \
                         hf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_x[k] * gg_179[k];

        t_180[k] = f_7 * gf_120[k]
                   + pa_x[k] * gg_180[k];

        t_181[k] = f_6 * gf_80[k]
                   + pb_y[k] * hf_120[k];

        t_182[k] = f_6 * gf_70[k]
                   + pb_z[k] * hf_120[k];

        t_183[k] = f_6 * gf_123[k]
                   + pa_x[k] * gg_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pb_x, pb_y, gf_82, gf_125, gf_126, \
                         gf_127, gg_185, hf_122, hf_126, hf_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_6 * gf_82[k]
                   + pb_y[k] * hf_122[k];

        t_185[k] = f_6 * gf_125[k]
                   + pa_x[k] * gg_185[k];

        t_186[k] = f_5 * gf_126[k]
                   + pb_x[k] * hf_126[k];

        t_187[k] = f_5 * gf_127[k]
                   + pb_x[k] * hf_127[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, t_193, pa_x, pb_x, gf_128, gf_129, \
                         gg_190, gg_191, gg_192, gg_193, hf_128, \
                         hf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_5 * gf_128[k]
                   + pb_x[k] * hf_128[k];

        t_189[k] = f_5 * gf_129[k]
                   + pb_x[k] * hf_129[k];

        t_190[k] = pa_x[k] * gg_190[k];

        t_191[k] = pa_x[k] * gg_191[k];

        t_192[k] = pa_x[k] * gg_192[k];

        t_193[k] = pa_x[k] * gg_193[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, pa_x, pa_y, pb_y, gf_90, gf_133, \
                         gg_135, gg_137, gg_194, gg_198, hf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * gg_194[k];

        t_195[k] = pa_y[k] * gg_135[k];

        t_196[k] = f_5 * gf_90[k]
                   + pb_y[k] * hf_130[k];

        t_197[k] = pa_y[k] * gg_137[k];

        t_198[k] = f_6 * gf_133[k]
                   + pa_x[k] * gg_198[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pb_x, pb_y, gf_92, gf_136, gf_137, \
                         gg_140, hf_132, hf_136, hf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_5 * gf_92[k]
                   + pb_y[k] * hf_132[k];

        t_200[k] = pa_y[k] * gg_140[k];

        t_201[k] = f_5 * gf_136[k]
                   + pb_x[k] * hf_136[k];

        t_202[k] = f_5 * gf_137[k]
                   + pb_x[k] * hf_137[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, pa_x, pa_y, pb_x, gf_138, \
                         gg_144, gg_205, gg_206, gg_207, gg_208, \
                         hf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_5 * gf_138[k]
                   + pb_x[k] * hf_138[k];

        t_204[k] = pa_y[k] * gg_144[k];

        t_205[k] = pa_x[k] * gg_205[k];

        t_206[k] = pa_x[k] * gg_206[k];

        t_207[k] = pa_x[k] * gg_207[k];

        t_208[k] = pa_x[k] * gg_208[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pa_x, pb_y, pb_z, gf_90, gf_140, \
                         gf_143, gg_209, gg_210, gg_213, hf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_x[k] * gg_209[k];

        t_210[k] = f_7 * gf_140[k]
                   + pa_x[k] * gg_210[k];

        t_211[k] = pb_y[k] * hf_140[k];

        t_212[k] = f_7 * gf_90[k]
                   + pb_z[k] * hf_140[k];

        t_213[k] = f_6 * gf_143[k]
                   + pa_x[k] * gg_213[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pa_x, pb_x, pb_y, gf_145, gf_146, \
                         gf_147, gg_215, hf_142, hf_145, hf_146, \
                         hf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * hf_142[k];

        t_215[k] = f_6 * gf_145[k]
                   + pa_x[k] * gg_215[k];

        t_216[k] = f_5 * gf_146[k]
                   + pb_x[k] * hf_146[k];

        t_217[k] = f_5 * gf_147[k]
                   + pb_x[k] * hf_147[k];

        t_218[k] = pb_y[k] * hf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, pa_x, pb_x, pb_y, gf_149, \
                         gg_220, gg_221, gg_222, gg_224, hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_5 * gf_149[k]
                   + pb_x[k] * hf_149[k];

        t_220[k] = pa_x[k] * gg_220[k];

        t_221[k] = pa_x[k] * gg_221[k];

        t_222[k] = pa_x[k] * gg_222[k];

        t_223[k] = pb_y[k] * hf_149[k];

        t_224[k] = pa_x[k] * gg_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, pb_x, pb_y, pb_z, gf_100, hd0_90, \
                         hd0_93, hd1_90, hd1_93, hf_150, hf_151, \
                         hf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_1 * hd0_90[k]
                   - f_2 * hd1_90[k]
                   + pb_x[k] * hf_150[k];

        t_226[k] = f_0 * gf_100[k]
                   + pb_y[k] * hf_150[k];

        t_227[k] = pb_z[k] * hf_150[k];

        t_228[k] = f_3 * hd0_93[k]
                   - f_4 * hd1_93[k]
                   + pb_x[k] * hf_153[k];

        t_229[k] = pb_z[k] * hf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, pb_x, hd0_95, hd1_95, hf_155, \
                         hf_156, hf_157, hf_158, hf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_3 * hd0_95[k]
                   - f_4 * hd1_95[k]
                   + pb_x[k] * hf_155[k];

        t_231[k] = pb_x[k] * hf_156[k];

        t_232[k] = pb_x[k] * hf_157[k];

        t_233[k] = pb_x[k] * hf_158[k];

        t_234[k] = pb_x[k] * hf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pb_y, pb_z, gf_106, gf_109, \
                         hd0_93, hd0_95, hd1_93, hd1_95, hf_156, hf_157, \
                         hf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_0 * gf_106[k]
                   + f_1 * hd0_93[k]
                   - f_2 * hd1_93[k]
                   + pb_y[k] * hf_156[k];

        t_236[k] = pb_z[k] * hf_156[k];

        t_237[k] = f_3 * hd0_93[k]
                   - f_4 * hd1_93[k]
                   + pb_z[k] * hf_157[k];

        t_238[k] = f_0 * gf_109[k]
                   + pb_y[k] * hf_159[k];

        t_239[k] = f_1 * hd0_95[k]
                   - f_2 * hd1_95[k]
                   + pb_z[k] * hf_159[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pa_z, pb_y, pb_z, gf_100, gf_112, \
                         gg_150, gg_151, gg_153, hf_160, hf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_z[k] * gg_150[k];

        t_241[k] = pa_z[k] * gg_151[k];

        t_242[k] = f_5 * gf_100[k]
                   + pb_z[k] * hf_160[k];

        t_243[k] = pa_z[k] * gg_153[k];

        t_244[k] = f_7 * gf_112[k]
                   + pb_y[k] * hf_162[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, t_250, pa_z, pb_x, gf_102, gg_155, \
                         gg_160, hf_166, hf_167, hf_168, hf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_6 * gf_102[k]
                   + pa_z[k] * gg_155[k];

        t_246[k] = pb_x[k] * hf_166[k];

        t_247[k] = pb_x[k] * hf_167[k];

        t_248[k] = pb_x[k] * hf_168[k];

        t_249[k] = pb_x[k] * hf_169[k];

        t_250[k] = pa_z[k] * gg_160[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_z, pb_y, pb_z, gf_106, gf_107, gf_109, \
                         gf_119, gg_162, gg_164, hf_166, hf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_5 * gf_106[k]
                   + pb_z[k] * hf_166[k];

        t_252[k] = f_6 * gf_107[k]
                   + pa_z[k] * gg_162[k];

        t_253[k] = f_7 * gf_119[k]
                   + pb_y[k] * hf_169[k];

        t_254[k] = f_7 * gf_109[k]
                   + pa_z[k] * gg_164[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pb_x, pb_y, pb_z, gf_110, gf_120, \
                         hd0_102, hd0_105, hd1_102, hd1_105, hf_170, \
                         hf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_1 * hd0_102[k]
                   - f_2 * hd1_102[k]
                   + pb_x[k] * hf_170[k];

        t_256[k] = f_10 * gf_120[k]
                   + pb_y[k] * hf_170[k];

        t_257[k] = f_6 * gf_110[k]
                   + pb_z[k] * hf_170[k];

        t_258[k] = f_3 * hd0_105[k]
                   - f_4 * hd1_105[k]
                   + pb_x[k] * hf_173[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, pb_x, pb_y, gf_122, hd0_107, \
                         hd1_107, hf_172, hf_175, hf_176, hf_177, \
                         hf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_10 * gf_122[k]
                   + pb_y[k] * hf_172[k];

        t_260[k] = f_3 * hd0_107[k]
                   - f_4 * hd1_107[k]
                   + pb_x[k] * hf_175[k];

        t_261[k] = pb_x[k] * hf_176[k];

        t_262[k] = pb_x[k] * hf_177[k];

        t_263[k] = pb_x[k] * hf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_z, pb_x, pb_z, fg0_100, fg1_100, gf_116, \
                         gg_175, hf_176, hf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = pb_x[k] * hf_179[k];

        t_265[k] = f_8 * fg0_100[k]
                   - f_9 * fg1_100[k]
                   + pa_z[k] * gg_175[k];

        t_266[k] = f_6 * gf_116[k]
                   + pb_z[k] * hf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_y, pb_y, fg0_134, fg1_134, gf_128, gf_129, \
                         gg_194, hd0_107, hd1_107, hf_178, hf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_10 * gf_128[k]
                   + f_3 * hd0_107[k]
                   - f_4 * hd1_107[k]
                   + pb_y[k] * hf_178[k];

        t_268[k] = f_10 * gf_129[k]
                   + pb_y[k] * hf_179[k];

        t_269[k] = f_11 * fg0_134[k]
                   - f_12 * fg1_134[k]
                   + pa_y[k] * gg_194[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_x, pb_y, pb_z, gf_120, gf_130, \
                         hd0_108, hd0_111, hd1_108, hd1_111, hf_180, \
                         hf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * hd0_108[k]
                   - f_2 * hd1_108[k]
                   + pb_x[k] * hf_180[k];

        t_271[k] = f_6 * gf_130[k]
                   + pb_y[k] * hf_180[k];

        t_272[k] = f_10 * gf_120[k]
                   + pb_z[k] * hf_180[k];

        t_273[k] = f_3 * hd0_111[k]
                   - f_4 * hd1_111[k]
                   + pb_x[k] * hf_183[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_x, pb_y, gf_132, hd0_113, \
                         hd1_113, hf_182, hf_185, hf_186, hf_187, \
                         hf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_6 * gf_132[k]
                   + pb_y[k] * hf_182[k];

        t_275[k] = f_3 * hd0_113[k]
                   - f_4 * hd1_113[k]
                   + pb_x[k] * hf_185[k];

        t_276[k] = pb_x[k] * hf_186[k];

        t_277[k] = pb_x[k] * hf_187[k];

        t_278[k] = pb_x[k] * hf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pa_z, pb_x, pb_z, fg0_115, fg1_115, gf_126, \
                         gg_190, hf_186, hf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_x[k] * hf_189[k];

        t_280[k] = f_11 * fg0_115[k]
                   - f_12 * fg1_115[k]
                   + pa_z[k] * gg_190[k];

        t_281[k] = f_10 * gf_126[k]
                   + pb_z[k] * hf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_y, pb_y, fg0_149, fg1_149, gf_138, \
                         gf_139, gg_209, gg_210, hd0_113, hd1_113, hf_188, \
                         hf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_6 * gf_138[k]
                   + f_3 * hd0_113[k]
                   - f_4 * hd1_113[k]
                   + pb_y[k] * hf_188[k];

        t_283[k] = f_6 * gf_139[k]
                   + pb_y[k] * hf_189[k];

        t_284[k] = f_8 * fg0_149[k]
                   - f_9 * fg1_149[k]
                   + pa_y[k] * gg_209[k];

        t_285[k] = pa_y[k] * gg_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, pa_y, pb_y, gf_140, gf_141, \
                         gf_142, gg_212, gg_213, gg_215, hf_190, \
                         hf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_5 * gf_140[k]
                   + pb_y[k] * hf_190[k];

        t_287[k] = pa_y[k] * gg_212[k];

        t_288[k] = f_6 * gf_141[k]
                   + pa_y[k] * gg_213[k];

        t_289[k] = f_5 * gf_142[k]
                   + pb_y[k] * hf_192[k];

        t_290[k] = pa_y[k] * gg_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, t_295, t_296, pa_y, pb_x, pb_z, gf_136, \
                         gf_146, gg_220, hf_196, hf_197, hf_198, \
                         hf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = pb_x[k] * hf_196[k];

        t_292[k] = pb_x[k] * hf_197[k];

        t_293[k] = pb_x[k] * hf_198[k];

        t_294[k] = pb_x[k] * hf_199[k];

        t_295[k] = f_7 * gf_146[k]
                   + pa_y[k] * gg_220[k];

        t_296[k] = f_7 * gf_136[k]
                   + pb_z[k] * hf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, pa_y, pb_x, pb_y, gf_148, gf_149, \
                         gg_222, gg_224, hd0_120, hd1_120, hf_199, \
                         hf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_6 * gf_148[k]
                   + pa_y[k] * gg_222[k];

        t_298[k] = f_5 * gf_149[k]
                   + pb_y[k] * hf_199[k];

        t_299[k] = pa_y[k] * gg_224[k];

        t_300[k] = f_1 * hd0_120[k]
                   - f_2 * hd1_120[k]
                   + pb_x[k] * hf_200[k];

        t_301[k] = pb_y[k] * hf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pb_x, pb_y, pb_z, gf_140, hd0_123, \
                         hd0_125, hd1_123, hd1_125, hf_200, hf_202, hf_203, \
                         hf_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * gf_140[k]
                   + pb_z[k] * hf_200[k];

        t_303[k] = f_3 * hd0_123[k]
                   - f_4 * hd1_123[k]
                   + pb_x[k] * hf_203[k];

        t_304[k] = pb_y[k] * hf_202[k];

        t_305[k] = f_3 * hd0_125[k]
                   - f_4 * hd1_125[k]
                   + pb_x[k] * hf_205[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, pb_x, pb_y, pb_z, gf_146, \
                         hd0_123, hd1_123, hf_206, hf_207, hf_208, \
                         hf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = pb_x[k] * hf_206[k];

        t_307[k] = pb_x[k] * hf_207[k];

        t_308[k] = pb_x[k] * hf_208[k];

        t_309[k] = pb_x[k] * hf_209[k];

        t_310[k] = f_1 * hd0_123[k]
                   - f_2 * hd1_123[k]
                   + pb_y[k] * hf_206[k];

        t_311[k] = f_0 * gf_146[k]
                   + pb_z[k] * hf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_y, pb_z, gf_149, hd0_125, hd1_125, hf_208, \
                         hf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_3 * hd0_125[k]
                   - f_4 * hd1_125[k]
                   + pb_y[k] * hf_208[k];

        t_313[k] = pb_y[k] * hf_209[k];

        t_314[k] = f_0 * gf_149[k]
                   + f_1 * hd0_125[k]
                   - f_2 * hd1_125[k]
                   + pb_z[k] * hf_209[k];
    }
}

}  // namespace simdt2ceri
