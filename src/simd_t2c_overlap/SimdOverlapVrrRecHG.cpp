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


#include "SimdOverlapVrrRecHG.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_hg_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fg, const size_t gf, const size_t gg,
                          const size_t hd, const size_t hf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_26 = buffer.data(fg + 26);

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
    const auto *gf_53 = buffer.data(gf + 53);
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
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);

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
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);

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
    const auto *hf_141 = buffer.data(hf + 141);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_143 = buffer.data(hf + 143);
    const auto *hf_144 = buffer.data(hf + 144);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_148 = buffer.data(hf + 148);
    const auto *hf_149 = buffer.data(hf + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gf_0, hd_0, hf_0, \
                         hf_1, hf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_2 * hd_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = pb_y[k] * hf_2[k];

        t_5[k] = f_2 * hd_0[k]
                 + pb_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, gf_3, gf_4, hd_1, \
                         hf_3, hf_4, hf_5, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gf_3[k]
                 + pb_x[k] * hf_5[k];

        t_7[k] = pb_z[k] * hf_3[k];

        t_8[k] = pb_y[k] * hf_4[k];

        t_9[k] = f_0 * gf_4[k]
                 + pb_x[k] * hf_7[k];

        t_10[k] = f_1 * hd_1[k]
                  + pb_y[k] * hf_5[k];

        t_11[k] = pb_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_y, pb_y, pb_z, gf_0, gg_0, \
                         hd_2, hf_6, hf_7, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * hd_2[k]
                  + pb_y[k] * hf_6[k];

        t_13[k] = pb_y[k] * hf_7[k];

        t_14[k] = f_1 * hd_2[k]
                  + pb_z[k] * hf_7[k];

        t_15[k] = pa_y[k] * gg_0[k];

        t_16[k] = f_2 * gf_0[k]
                  + pb_y[k] * hf_8[k];

        t_17[k] = pb_z[k] * hf_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pb_x, pb_z, gf_1, gf_6, gg_1, \
                         gg_2, hf_9, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gf_1[k]
                  + pa_y[k] * gg_1[k];

        t_19[k] = pb_z[k] * hf_9[k];

        t_20[k] = pa_y[k] * gg_2[k];

        t_21[k] = f_4 * gf_6[k]
                  + pb_x[k] * hf_11[k];

        t_22[k] = pb_z[k] * hf_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, pb_x, pb_z, fg_4, gf_7, gg_4, \
                         gg_11, hf_11, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gf_7[k]
                  + pb_x[k] * hf_13[k];

        t_24[k] = pa_y[k] * gg_4[k];

        t_25[k] = f_1 * fg_4[k]
                  + pa_x[k] * gg_11[k];

        t_26[k] = pb_z[k] * hf_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, pb_z, gf_4, gg_0, \
                         gg_6, hd_4, hf_12, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * hd_4[k]
                  + pb_z[k] * hf_12[k];

        t_28[k] = f_2 * gf_4[k]
                  + pb_y[k] * hf_14[k];

        t_29[k] = pa_y[k] * gg_6[k];

        t_30[k] = pa_z[k] * gg_0[k];

        t_31[k] = pb_y[k] * hf_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, gf_0, gf_2, gg_1, \
                         gg_2, gg_3, hf_15, hf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * gf_0[k]
                  + pb_z[k] * hf_15[k];

        t_33[k] = pa_z[k] * gg_1[k];

        t_34[k] = pb_y[k] * hf_16[k];

        t_35[k] = f_3 * gf_2[k]
                  + pa_z[k] * gg_2[k];

        t_36[k] = pa_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf_11, gf_12, gg_5, \
                         hd_7, hf_17, hf_18, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_4 * gf_11[k]
                  + pb_x[k] * hf_18[k];

        t_38[k] = pb_y[k] * hf_17[k];

        t_39[k] = f_4 * gf_12[k]
                  + pb_x[k] * hf_20[k];

        t_40[k] = pa_z[k] * gg_5[k];

        t_41[k] = f_3 * hd_7[k]
                  + pb_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, fg_0, fg_7, gg_7, gg_16, \
                         hd_8, hf_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * hd_8[k]
                  + pb_y[k] * hf_19[k];

        t_43[k] = pb_y[k] * hf_20[k];

        t_44[k] = f_1 * fg_7[k]
                  + pa_x[k] * gg_16[k];

        t_45[k] = f_2 * fg_0[k]
                  + pa_y[k] * gg_7[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, gf_5, gf_14, hd_9, \
                         hd_10, hf_21, hf_22, hf_23, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * gf_5[k]
                  + pb_y[k] * hf_21[k];

        t_47[k] = pb_z[k] * hf_21[k];

        t_48[k] = f_1 * gf_14[k]
                  + f_2 * hd_10[k]
                  + pb_x[k] * hf_24[k];

        t_49[k] = pb_z[k] * hf_22[k];

        t_50[k] = f_2 * hd_9[k]
                  + pb_z[k] * hf_23[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_z, gf_15, gf_16, gf_17, hf_24, \
                         hf_25, hf_27, hf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * gf_15[k]
                  + pb_x[k] * hf_25[k];

        t_52[k] = pb_z[k] * hf_24[k];

        t_53[k] = f_1 * gf_16[k]
                  + pb_x[k] * hf_27[k];

        t_54[k] = f_1 * gf_17[k]
                  + pb_x[k] * hf_28[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, fg_8, gf_8, gg_21, \
                         hd_10, hd_11, hf_25, hf_26, hf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_3 * fg_8[k]
                  + pa_x[k] * gg_21[k];

        t_56[k] = pb_z[k] * hf_25[k];

        t_57[k] = f_2 * hd_10[k]
                  + pb_z[k] * hf_26[k];

        t_58[k] = f_3 * gf_8[k]
                  + pb_y[k] * hf_28[k];

        t_59[k] = f_1 * hd_11[k]
                  + pb_z[k] * hf_28[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, gf_10, gg_8, \
                         gg_9, gg_12, gg_13, gg_14, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * gg_12[k];

        t_61[k] = pa_z[k] * gg_8[k];

        t_62[k] = pa_y[k] * gg_13[k];

        t_63[k] = pa_z[k] * gg_9[k];

        t_64[k] = f_2 * gf_10[k]
                  + pb_y[k] * hf_29[k];

        t_65[k] = pa_y[k] * gg_14[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, gf_20, gf_21, gg_10, \
                         gg_11, gg_15, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * gg_10[k];

        t_67[k] = f_1 * gf_20[k]
                  + pb_x[k] * hf_31[k];

        t_68[k] = f_1 * gf_21[k]
                  + pb_x[k] * hf_32[k];

        t_69[k] = pa_y[k] * gg_15[k];

        t_70[k] = pa_z[k] * gg_11[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pa_y, pb_y, pb_z, fg_9, gf_6, gf_12, \
                         gg_16, gg_23, hf_30, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_2 * gf_6[k]
                  + pb_z[k] * hf_30[k];

        t_72[k] = f_3 * fg_9[k]
                  + pa_x[k] * gg_23[k];

        t_73[k] = f_2 * gf_12[k]
                  + pb_y[k] * hf_33[k];

        t_74[k] = pa_y[k] * gg_16[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, fg_0, gf_9, gg_12, \
                         hd_14, hf_34, hf_35, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * fg_0[k]
                  + pa_z[k] * gg_12[k];

        t_76[k] = pb_y[k] * hf_34[k];

        t_77[k] = f_3 * gf_9[k]
                  + pb_z[k] * hf_34[k];

        t_78[k] = f_2 * hd_14[k]
                  + pb_y[k] * hf_35[k];

        t_79[k] = pb_y[k] * hf_36[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pb_x, pb_y, gf_26, gf_27, gf_28, gf_29, \
                         hd_17, hf_37, hf_38, hf_39, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * gf_26[k]
                  + f_2 * hd_17[k]
                  + pb_x[k] * hf_37[k];

        t_81[k] = f_1 * gf_27[k]
                  + pb_x[k] * hf_38[k];

        t_82[k] = f_1 * gf_28[k]
                  + pb_x[k] * hf_39[k];

        t_83[k] = pb_y[k] * hf_37[k];

        t_84[k] = f_1 * gf_29[k]
                  + pb_x[k] * hf_41[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_x, pb_y, fg_10, gg_29, hd_15, hd_16, \
                         hd_17, hf_38, hf_39, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * hd_15[k]
                  + pb_y[k] * hf_38[k];

        t_86[k] = f_3 * hd_16[k]
                  + pb_y[k] * hf_39[k];

        t_87[k] = f_2 * hd_17[k]
                  + pb_y[k] * hf_40[k];

        t_88[k] = pb_y[k] * hf_41[k];

        t_89[k] = f_3 * fg_10[k]
                  + pa_x[k] * gg_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_y, pb_x, pb_y, pb_z, fg_3, gf_13, gf_31, \
                         gg_17, hd_19, hf_42, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * fg_3[k]
                  + pa_y[k] * gg_17[k];

        t_91[k] = f_1 * gf_13[k]
                  + pb_y[k] * hf_42[k];

        t_92[k] = pb_z[k] * hf_42[k];

        t_93[k] = f_3 * gf_31[k]
                  + f_2 * hd_19[k]
                  + pb_x[k] * hf_45[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_z, gf_32, gf_33, hd_18, hf_43, \
                         hf_44, hf_45, hf_46, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pb_z[k] * hf_43[k];

        t_95[k] = f_2 * hd_18[k]
                  + pb_z[k] * hf_44[k];

        t_96[k] = f_3 * gf_32[k]
                  + pb_x[k] * hf_46[k];

        t_97[k] = pb_z[k] * hf_45[k];

        t_98[k] = f_3 * gf_33[k]
                  + pb_x[k] * hf_48[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pb_x, pb_z, fg_13, gf_34, gg_34, \
                         hd_19, hf_46, hf_47, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_3 * gf_34[k]
                  + pb_x[k] * hf_49[k];

        t_100[k] = f_2 * fg_13[k]
                   + pa_x[k] * gg_34[k];

        t_101[k] = pb_z[k] * hf_46[k];

        t_102[k] = f_2 * hd_19[k]
                   + pb_z[k] * hf_47[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, gf_13, gf_17, \
                         gg_17, gg_18, hd_20, hf_49, hf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * gf_17[k]
                   + pb_y[k] * hf_49[k];

        t_104[k] = f_1 * hd_20[k]
                   + pb_z[k] * hf_49[k];

        t_105[k] = pa_z[k] * gg_17[k];

        t_106[k] = pa_z[k] * gg_18[k];

        t_107[k] = f_2 * gf_13[k]
                   + pb_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pa_z, pb_y, fg_6, gf_18, gg_19, \
                         gg_20, gg_22, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * gg_19[k];

        t_109[k] = f_3 * gf_18[k]
                   + pb_y[k] * hf_51[k];

        t_110[k] = f_2 * fg_6[k]
                   + pa_y[k] * gg_22[k];

        t_111[k] = pa_z[k] * gg_20[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_z, pb_x, gf_37, gf_38, gf_39, gg_21, \
                         hf_53, hf_54, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_3 * gf_37[k]
                   + pb_x[k] * hf_53[k];

        t_113[k] = f_3 * gf_38[k]
                   + pb_x[k] * hf_54[k];

        t_114[k] = f_3 * gf_39[k]
                   + pb_x[k] * hf_55[k];

        t_115[k] = pa_z[k] * gg_21[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_x, pb_y, pb_z, fg_16, fg_17, gf_15, \
                         gf_22, gg_35, gg_36, hf_52, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_2 * gf_15[k]
                   + pb_z[k] * hf_52[k];

        t_117[k] = f_2 * fg_16[k]
                   + pa_x[k] * gg_35[k];

        t_118[k] = f_3 * gf_22[k]
                   + pb_y[k] * hf_55[k];

        t_119[k] = f_2 * fg_17[k]
                   + pa_x[k] * gg_36[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pa_y, pb_y, gf_23, gf_24, gf_25, \
                         gg_24, gg_25, gg_26, hf_56, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * gg_24[k];

        t_121[k] = f_2 * gf_23[k]
                   + pb_y[k] * hf_56[k];

        t_122[k] = pa_y[k] * gg_25[k];

        t_123[k] = f_3 * gf_24[k]
                   + pa_y[k] * gg_26[k];

        t_124[k] = f_2 * gf_25[k]
                   + pb_y[k] * hf_57[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_y, pb_x, gf_42, gf_43, gf_44, \
                         gg_27, gg_28, hf_58, hf_59, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_y[k] * gg_27[k];

        t_126[k] = f_3 * gf_42[k]
                   + pb_x[k] * hf_58[k];

        t_127[k] = f_3 * gf_43[k]
                   + pb_x[k] * hf_59[k];

        t_128[k] = f_3 * gf_44[k]
                   + pb_x[k] * hf_60[k];

        t_129[k] = pa_y[k] * gg_28[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_x, pb_y, pb_z, fg_18, fg_19, gf_19, \
                         gf_29, gg_37, gg_38, hf_58, hf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * fg_18[k]
                   + pa_x[k] * gg_37[k];

        t_131[k] = f_3 * gf_19[k]
                   + pb_z[k] * hf_58[k];

        t_132[k] = f_2 * fg_19[k]
                   + pa_x[k] * gg_38[k];

        t_133[k] = f_2 * gf_29[k]
                   + pb_y[k] * hf_61[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_y, pa_z, pb_y, pb_z, fg_5, \
                         gf_23, gg_24, gg_29, hd_26, hf_62, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_y[k] * gg_29[k];

        t_135[k] = f_3 * fg_5[k]
                   + pa_z[k] * gg_24[k];

        t_136[k] = pb_y[k] * hf_62[k];

        t_137[k] = f_1 * gf_23[k]
                   + pb_z[k] * hf_62[k];

        t_138[k] = f_2 * hd_26[k]
                   + pb_y[k] * hf_63[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, gf_47, gf_48, gf_49, \
                         hd_29, hf_64, hf_65, hf_66, hf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * hf_64[k];

        t_140[k] = f_3 * gf_47[k]
                   + f_2 * hd_29[k]
                   + pb_x[k] * hf_65[k];

        t_141[k] = f_3 * gf_48[k]
                   + pb_x[k] * hf_66[k];

        t_142[k] = f_3 * gf_49[k]
                   + pb_x[k] * hf_67[k];

        t_143[k] = pb_y[k] * hf_65[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pb_x, pb_y, gf_50, hd_27, hd_28, \
                         hd_29, hf_66, hf_67, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_3 * gf_50[k]
                   + pb_x[k] * hf_69[k];

        t_145[k] = f_1 * hd_27[k]
                   + pb_y[k] * hf_66[k];

        t_146[k] = f_3 * hd_28[k]
                   + pb_y[k] * hf_67[k];

        t_147[k] = f_2 * hd_29[k]
                   + pb_y[k] * hf_68[k];

        t_148[k] = pb_y[k] * hf_69[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pa_x, pb_y, pb_z, fg_26, gf_30, \
                         gf_51, gf_53, gg_43, gg_44, gg_46, hf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_2 * fg_26[k]
                   + pa_x[k] * gg_43[k];

        t_150[k] = f_4 * gf_51[k]
                   + pa_x[k] * gg_44[k];

        t_151[k] = f_4 * gf_30[k]
                   + pb_y[k] * hf_70[k];

        t_152[k] = pb_z[k] * hf_70[k];

        t_153[k] = f_3 * gf_53[k]
                   + pa_x[k] * gg_46[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, pb_x, pb_z, gf_55, gf_57, hd_30, \
                         hf_71, hf_72, hf_73, hf_74, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * hf_71[k];

        t_155[k] = f_2 * hd_30[k]
                   + pb_z[k] * hf_72[k];

        t_156[k] = f_2 * gf_55[k]
                   + pb_x[k] * hf_74[k];

        t_157[k] = pb_z[k] * hf_73[k];

        t_158[k] = f_2 * gf_57[k]
                   + pb_x[k] * hf_75[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, t_164, pa_x, pb_x, pb_z, gf_58, \
                         gg_49, gg_50, gg_51, gg_52, hf_74, hf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_2 * gf_58[k]
                   + pb_x[k] * hf_76[k];

        t_160[k] = pa_x[k] * gg_49[k];

        t_161[k] = pb_z[k] * hf_74[k];

        t_162[k] = pa_x[k] * gg_50[k];

        t_163[k] = pa_x[k] * gg_51[k];

        t_164[k] = pa_x[k] * gg_52[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pa_z, pb_y, pb_z, gf_30, gf_36, \
                         gg_30, gg_31, gg_32, hf_77, hf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_z[k] * gg_30[k];

        t_166[k] = pa_z[k] * gg_31[k];

        t_167[k] = f_2 * gf_30[k]
                   + pb_z[k] * hf_77[k];

        t_168[k] = pa_z[k] * gg_32[k];

        t_169[k] = f_1 * gf_36[k]
                   + pb_y[k] * hf_78[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, pa_z, pb_x, gf_59, gf_61, gf_62, \
                         gg_33, gg_53, hf_79, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_3 * gf_59[k]
                   + pa_x[k] * gg_53[k];

        t_171[k] = pa_z[k] * gg_33[k];

        t_172[k] = f_2 * gf_61[k]
                   + pb_x[k] * hf_79[k];

        t_173[k] = f_2 * gf_62[k]
                   + pb_x[k] * hf_80[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, pa_x, pb_x, gf_63, gg_54, \
                         gg_55, gg_56, gg_57, gg_58, hf_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_2 * gf_63[k]
                   + pb_x[k] * hf_81[k];

        t_175[k] = pa_x[k] * gg_54[k];

        t_176[k] = pa_x[k] * gg_55[k];

        t_177[k] = pa_x[k] * gg_56[k];

        t_178[k] = pa_x[k] * gg_57[k];

        t_179[k] = pa_x[k] * gg_58[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pb_y, pb_z, gf_35, gf_40, gf_64, \
                         gf_65, gg_59, gg_60, hf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * gf_64[k]
                   + pa_x[k] * gg_59[k];

        t_181[k] = f_3 * gf_40[k]
                   + pb_y[k] * hf_82[k];

        t_182[k] = f_3 * gf_35[k]
                   + pb_z[k] * hf_82[k];

        t_183[k] = f_3 * gf_65[k]
                   + pa_x[k] * gg_60[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pb_x, pb_y, gf_41, gf_66, gf_67, \
                         gf_68, gg_61, hf_83, hf_84, hf_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * gf_41[k]
                   + pb_y[k] * hf_83[k];

        t_185[k] = f_3 * gf_66[k]
                   + pa_x[k] * gg_61[k];

        t_186[k] = f_2 * gf_67[k]
                   + pb_x[k] * hf_84[k];

        t_187[k] = f_2 * gf_68[k]
                   + pb_x[k] * hf_85[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, t_193, pa_x, pb_x, gf_69, gf_70, \
                         gg_62, gg_63, gg_64, gg_65, hf_86, hf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_2 * gf_69[k]
                   + pb_x[k] * hf_86[k];

        t_189[k] = f_2 * gf_70[k]
                   + pb_x[k] * hf_87[k];

        t_190[k] = pa_x[k] * gg_62[k];

        t_191[k] = pa_x[k] * gg_63[k];

        t_192[k] = pa_x[k] * gg_64[k];

        t_193[k] = pa_x[k] * gg_65[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, pa_x, pa_y, pb_y, gf_45, gf_71, \
                         gg_39, gg_40, gg_66, gg_67, hf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * gg_66[k];

        t_195[k] = pa_y[k] * gg_39[k];

        t_196[k] = f_2 * gf_45[k]
                   + pb_y[k] * hf_88[k];

        t_197[k] = pa_y[k] * gg_40[k];

        t_198[k] = f_3 * gf_71[k]
                   + pa_x[k] * gg_67[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pb_x, pb_y, gf_46, gf_72, gf_73, \
                         gg_41, hf_89, hf_90, hf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_2 * gf_46[k]
                   + pb_y[k] * hf_89[k];

        t_200[k] = pa_y[k] * gg_41[k];

        t_201[k] = f_2 * gf_72[k]
                   + pb_x[k] * hf_90[k];

        t_202[k] = f_2 * gf_73[k]
                   + pb_x[k] * hf_91[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, pa_x, pa_y, pb_x, gf_74, \
                         gg_42, gg_68, gg_69, gg_70, gg_71, hf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_2 * gf_74[k]
                   + pb_x[k] * hf_92[k];

        t_204[k] = pa_y[k] * gg_42[k];

        t_205[k] = pa_x[k] * gg_68[k];

        t_206[k] = pa_x[k] * gg_69[k];

        t_207[k] = pa_x[k] * gg_70[k];

        t_208[k] = pa_x[k] * gg_71[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pa_x, pb_y, pb_z, gf_45, gf_76, \
                         gg_72, gg_73, hd_34, hf_93, hf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_x[k] * gg_72[k];

        t_210[k] = f_4 * gf_76[k]
                   + pa_x[k] * gg_73[k];

        t_211[k] = pb_y[k] * hf_93[k];

        t_212[k] = f_4 * gf_45[k]
                   + pb_z[k] * hf_93[k];

        t_213[k] = f_2 * hd_34[k]
                   + pb_y[k] * hf_94[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pa_x, pb_x, pb_y, gf_80, gf_81, \
                         gf_82, gg_78, hf_95, hf_96, hf_97, hf_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * hf_95[k];

        t_215[k] = f_3 * gf_80[k]
                   + pa_x[k] * gg_78[k];

        t_216[k] = f_2 * gf_81[k]
                   + pb_x[k] * hf_97[k];

        t_217[k] = f_2 * gf_82[k]
                   + pb_x[k] * hf_98[k];

        t_218[k] = pb_y[k] * hf_96[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, pa_x, pb_x, pb_y, gf_84, \
                         gg_79, gg_80, gg_81, gg_82, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_2 * gf_84[k]
                   + pb_x[k] * hf_99[k];

        t_220[k] = pa_x[k] * gg_79[k];

        t_221[k] = pa_x[k] * gg_80[k];

        t_222[k] = pa_x[k] * gg_81[k];

        t_223[k] = pb_y[k] * hf_99[k];

        t_224[k] = pa_x[k] * gg_82[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, t_230, pb_x, pb_z, hd_36, hd_37, \
                         hd_38, hd_39, hf_100, hf_101, hf_102, hf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_1 * hd_36[k]
                   + pb_x[k] * hf_100[k];

        t_226[k] = f_3 * hd_37[k]
                   + pb_x[k] * hf_101[k];

        t_227[k] = pb_z[k] * hf_100[k];

        t_228[k] = f_2 * hd_38[k]
                   + pb_x[k] * hf_102[k];

        t_229[k] = pb_z[k] * hf_101[k];

        t_230[k] = f_2 * hd_39[k]
                   + pb_x[k] * hf_103[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, pb_x, pb_y, pb_z, \
                         gf_55, hd_38, hf_104, hf_105, hf_106, hf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pb_x[k] * hf_104[k];

        t_232[k] = pb_x[k] * hf_105[k];

        t_233[k] = pb_x[k] * hf_106[k];

        t_234[k] = pb_x[k] * hf_107[k];

        t_235[k] = f_0 * gf_55[k]
                   + f_1 * hd_38[k]
                   + pb_y[k] * hf_104[k];

        t_236[k] = pb_z[k] * hf_104[k];

        t_237[k] = f_2 * hd_38[k]
                   + pb_z[k] * hf_105[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pa_z, pb_x, pb_y, pb_z, gf_58, \
                         gg_44, gg_45, hd_39, hd_40, hf_107, hf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * gf_58[k]
                   + pb_y[k] * hf_107[k];

        t_239[k] = f_1 * hd_39[k]
                   + pb_z[k] * hf_107[k];

        t_240[k] = pa_z[k] * gg_44[k];

        t_241[k] = pa_z[k] * gg_45[k];

        t_242[k] = f_3 * hd_40[k]
                   + pb_x[k] * hf_108[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, t_248, pa_z, pb_x, gg_46, hd_42, \
                         hd_43, hf_109, hf_110, hf_111, hf_112, \
                         hf_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * gg_46[k];

        t_244[k] = f_2 * hd_42[k]
                   + pb_x[k] * hf_109[k];

        t_245[k] = f_2 * hd_43[k]
                   + pb_x[k] * hf_110[k];

        t_246[k] = pb_x[k] * hf_111[k];

        t_247[k] = pb_x[k] * hf_112[k];

        t_248[k] = pb_x[k] * hf_113[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, pa_z, pb_x, pb_y, pb_z, gf_55, \
                         gf_56, gf_63, gg_49, gg_50, hf_111, hf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pb_x[k] * hf_114[k];

        t_250[k] = pa_z[k] * gg_49[k];

        t_251[k] = f_2 * gf_55[k]
                   + pb_z[k] * hf_111[k];

        t_252[k] = f_3 * gf_56[k]
                   + pa_z[k] * gg_50[k];

        t_253[k] = f_4 * gf_63[k]
                   + pb_y[k] * hf_114[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pa_y, pb_x, fg_17, gg_58, hd_44, hd_45, \
                         hd_46, hf_115, hf_116, hf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_1 * fg_17[k]
                   + pa_y[k] * gg_58[k];

        t_255[k] = f_1 * hd_44[k]
                   + pb_x[k] * hf_115[k];

        t_256[k] = f_3 * hd_45[k]
                   + pb_x[k] * hf_116[k];

        t_257[k] = f_3 * hd_46[k]
                   + pb_x[k] * hf_117[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, t_263, pb_x, hd_47, hd_48, hd_49, \
                         hf_118, hf_119, hf_120, hf_121, hf_122, \
                         hf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_2 * hd_47[k]
                   + pb_x[k] * hf_118[k];

        t_259[k] = f_2 * hd_48[k]
                   + pb_x[k] * hf_119[k];

        t_260[k] = f_2 * hd_49[k]
                   + pb_x[k] * hf_120[k];

        t_261[k] = pb_x[k] * hf_121[k];

        t_262[k] = pb_x[k] * hf_122[k];

        t_263[k] = pb_x[k] * hf_123[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_z, pb_x, pb_y, pb_z, fg_13, gf_60, \
                         gf_69, gg_54, hd_49, hf_121, hf_123, hf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = pb_x[k] * hf_124[k];

        t_265[k] = f_2 * fg_13[k]
                   + pa_z[k] * gg_54[k];

        t_266[k] = f_3 * gf_60[k]
                   + pb_z[k] * hf_121[k];

        t_267[k] = f_1 * gf_69[k]
                   + f_2 * hd_49[k]
                   + pb_y[k] * hf_123[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pb_x, pb_y, fg_20, gf_70, gg_66, \
                         hd_50, hd_51, hf_124, hf_125, hf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_1 * gf_70[k]
                   + pb_y[k] * hf_124[k];

        t_269[k] = f_3 * fg_20[k]
                   + pa_y[k] * gg_66[k];

        t_270[k] = f_1 * hd_50[k]
                   + pb_x[k] * hf_125[k];

        t_271[k] = f_3 * hd_51[k]
                   + pb_x[k] * hf_126[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, pb_x, hd_52, hd_53, hd_54, hd_55, \
                         hf_127, hf_128, hf_129, hf_130, hf_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_3 * hd_52[k]
                   + pb_x[k] * hf_127[k];

        t_273[k] = f_2 * hd_53[k]
                   + pb_x[k] * hf_128[k];

        t_274[k] = f_2 * hd_54[k]
                   + pb_x[k] * hf_129[k];

        t_275[k] = f_2 * hd_55[k]
                   + pb_x[k] * hf_130[k];

        t_276[k] = pb_x[k] * hf_131[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, pa_z, pb_x, pb_z, fg_15, gf_67, \
                         gg_62, hf_131, hf_132, hf_133, hf_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = pb_x[k] * hf_132[k];

        t_278[k] = pb_x[k] * hf_133[k];

        t_279[k] = pb_x[k] * hf_134[k];

        t_280[k] = f_3 * fg_15[k]
                   + pa_z[k] * gg_62[k];

        t_281[k] = f_1 * gf_67[k]
                   + pb_z[k] * hf_131[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_y, pb_y, fg_26, gf_74, gf_75, gg_72, \
                         gg_73, hd_55, hf_133, hf_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_3 * gf_74[k]
                   + f_2 * hd_55[k]
                   + pb_y[k] * hf_133[k];

        t_283[k] = f_3 * gf_75[k]
                   + pb_y[k] * hf_134[k];

        t_284[k] = f_2 * fg_26[k]
                   + pa_y[k] * gg_72[k];

        t_285[k] = pa_y[k] * gg_73[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, pa_y, pb_x, gg_75, gg_78, hd_56, \
                         hd_57, hd_58, hf_135, hf_136, hf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_3 * hd_56[k]
                   + pb_x[k] * hf_135[k];

        t_287[k] = pa_y[k] * gg_75[k];

        t_288[k] = f_2 * hd_57[k]
                   + pb_x[k] * hf_136[k];

        t_289[k] = f_2 * hd_58[k]
                   + pb_x[k] * hf_137[k];

        t_290[k] = pa_y[k] * gg_78[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, t_295, t_296, pa_y, pb_x, pb_z, gf_72, \
                         gf_81, gg_79, hf_138, hf_139, hf_140, hf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = pb_x[k] * hf_138[k];

        t_292[k] = pb_x[k] * hf_139[k];

        t_293[k] = pb_x[k] * hf_140[k];

        t_294[k] = pb_x[k] * hf_141[k];

        t_295[k] = f_4 * gf_81[k]
                   + pa_y[k] * gg_79[k];

        t_296[k] = f_4 * gf_72[k]
                   + pb_z[k] * hf_138[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, pa_y, pb_x, pb_y, gf_83, gf_84, \
                         gg_81, gg_82, hd_60, hf_141, hf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_3 * gf_83[k]
                   + pa_y[k] * gg_81[k];

        t_298[k] = f_2 * gf_84[k]
                   + pb_y[k] * hf_141[k];

        t_299[k] = pa_y[k] * gg_82[k];

        t_300[k] = f_1 * hd_60[k]
                   + pb_x[k] * hf_142[k];

        t_301[k] = pb_y[k] * hf_142[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, pb_x, pb_y, hd_61, hd_62, \
                         hd_64, hf_143, hf_144, hf_145, hf_146, \
                         hf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_3 * hd_61[k]
                   + pb_x[k] * hf_143[k];

        t_303[k] = f_2 * hd_62[k]
                   + pb_x[k] * hf_144[k];

        t_304[k] = pb_y[k] * hf_143[k];

        t_305[k] = f_2 * hd_64[k]
                   + pb_x[k] * hf_145[k];

        t_306[k] = pb_x[k] * hf_146[k];

        t_307[k] = pb_x[k] * hf_147[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, pb_x, pb_y, hd_62, hd_63, \
                         hd_64, hf_146, hf_147, hf_148, hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_x[k] * hf_148[k];

        t_309[k] = pb_x[k] * hf_149[k];

        t_310[k] = f_1 * hd_62[k]
                   + pb_y[k] * hf_146[k];

        t_311[k] = f_3 * hd_63[k]
                   + pb_y[k] * hf_147[k];

        t_312[k] = f_2 * hd_64[k]
                   + pb_y[k] * hf_148[k];

        t_313[k] = pb_y[k] * hf_149[k];
    }

#pragma omp simd aligned(t_314, pb_z, gf_84, hd_64, hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_0 * gf_84[k]
                   + f_1 * hd_64[k]
                   + pb_z[k] * hf_149[k];
    }
}

auto
compute_prim_hg_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fg, const size_t gf, const size_t gg,
                          const size_t hd, const size_t hf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_55 = buffer.data(fg + 55);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_41 = buffer.data(gf + 41);
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
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_129 = buffer.data(gg + 129);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gf_0, gf_3, hd_0, \
                         hf_0, hf_1, hf_2, hf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_2 * hd_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_2 * hd_0[k]
                 + pb_z[k] * hf_2[k];

        t_5[k] = f_0 * gf_3[k]
                 + pb_x[k] * hf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, gf_5, hd_1, hd_2, hf_3, \
                         hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * gf_5[k]
                 + pb_x[k] * hf_5[k];

        t_7[k] = f_1 * hd_1[k]
                 + pb_y[k] * hf_3[k];

        t_8[k] = f_2 * hd_2[k]
                 + pb_y[k] * hf_4[k];

        t_9[k] = pb_y[k] * hf_5[k];

        t_10[k] = f_1 * hd_2[k]
                  + pb_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, gf_0, gf_1, gf_7, \
                         gg_0, gg_3, gg_4, hf_6, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * gg_0[k];

        t_12[k] = f_2 * gf_0[k]
                  + pb_y[k] * hf_6[k];

        t_13[k] = f_3 * gf_1[k]
                  + pa_y[k] * gg_3[k];

        t_14[k] = pa_y[k] * gg_4[k];

        t_15[k] = f_4 * gf_7[k]
                  + pb_x[k] * hf_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, fg_9, gf_5, gg_11, hd_3, \
                         hf_7, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * fg_9[k]
                  + pa_x[k] * gg_11[k];

        t_17[k] = pb_z[k] * hf_7[k];

        t_18[k] = f_2 * hd_3[k]
                  + pb_z[k] * hf_8[k];

        t_19[k] = f_2 * gf_5[k]
                  + pb_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_y, pb_z, gf_0, gf_2, \
                         gg_0, gg_4, gg_7, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * gg_7[k];

        t_21[k] = pa_z[k] * gg_0[k];

        t_22[k] = f_2 * gf_0[k]
                  + pb_z[k] * hf_10[k];

        t_23[k] = pb_y[k] * hf_11[k];

        t_24[k] = f_3 * gf_2[k]
                  + pa_z[k] * gg_4[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, pb_y, fg_13, gf_13, gg_20, \
                         hd_5, hd_6, hf_12, hf_13, hf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * gf_13[k]
                  + pb_x[k] * hf_14[k];

        t_26[k] = f_3 * hd_5[k]
                  + pb_y[k] * hf_12[k];

        t_27[k] = f_2 * hd_6[k]
                  + pb_y[k] * hf_13[k];

        t_28[k] = pb_y[k] * hf_14[k];

        t_29[k] = f_1 * fg_13[k]
                  + pa_x[k] * gg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, fg_0, gf_6, gf_16, \
                         gg_8, hd_8, hf_15, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * fg_0[k]
                  + pa_y[k] * gg_8[k];

        t_31[k] = f_3 * gf_6[k]
                  + pb_y[k] * hf_15[k];

        t_32[k] = pb_z[k] * hf_15[k];

        t_33[k] = f_1 * gf_16[k]
                  + f_2 * hd_8[k]
                  + pb_x[k] * hf_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, pb_z, fg_17, gf_17, gg_25, \
                         hd_7, hd_8, hf_16, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * hd_7[k]
                  + pb_z[k] * hf_16[k];

        t_35[k] = f_1 * gf_17[k]
                  + pb_x[k] * hf_18[k];

        t_36[k] = f_3 * fg_17[k]
                  + pa_x[k] * gg_25[k];

        t_37[k] = pb_z[k] * hf_18[k];

        t_38[k] = f_2 * hd_8[k]
                  + pb_z[k] * hf_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pa_z, pb_y, pb_z, gf_9, gg_9, \
                         gg_16, gg_17, hd_9, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * gf_9[k]
                  + pb_y[k] * hf_20[k];

        t_40[k] = f_1 * hd_9[k]
                  + pb_z[k] * hf_20[k];

        t_41[k] = pa_y[k] * gg_16[k];

        t_42[k] = pa_z[k] * gg_9[k];

        t_43[k] = pa_y[k] * gg_17[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_z, pb_y, pb_z, fg_18, gf_7, gf_13, \
                         gg_11, gg_34, hf_21, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_z[k] * gg_11[k];

        t_45[k] = f_2 * gf_7[k]
                  + pb_z[k] * hf_21[k];

        t_46[k] = f_3 * fg_18[k]
                  + pa_x[k] * gg_34[k];

        t_47[k] = f_2 * gf_13[k]
                  + pb_y[k] * hf_22[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, pb_z, fg_0, gf_10, \
                         gg_15, gg_20, hd_10, hf_23, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * gg_20[k];

        t_49[k] = f_2 * fg_0[k]
                  + pa_z[k] * gg_15[k];

        t_50[k] = pb_y[k] * hf_23[k];

        t_51[k] = f_3 * gf_10[k]
                  + pb_z[k] * hf_23[k];

        t_52[k] = f_2 * hd_10[k]
                  + pb_y[k] * hf_24[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, gf_24, gf_28, hd_11, hd_13, \
                         hf_25, hf_26, hf_27, hf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * hf_25[k];

        t_54[k] = f_1 * gf_24[k]
                  + f_2 * hd_13[k]
                  + pb_x[k] * hf_26[k];

        t_55[k] = f_1 * gf_28[k]
                  + pb_x[k] * hf_30[k];

        t_56[k] = f_1 * hd_11[k]
                  + pb_y[k] * hf_27[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_y, fg_23, gg_46, hd_12, hd_13, \
                         hf_28, hf_29, hf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * hd_12[k]
                  + pb_y[k] * hf_28[k];

        t_58[k] = f_2 * hd_13[k]
                  + pb_y[k] * hf_29[k];

        t_59[k] = pb_y[k] * hf_30[k];

        t_60[k] = f_3 * fg_23[k]
                  + pa_x[k] * gg_46[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pb_x, pb_y, pb_z, fg_7, gf_14, gf_31, \
                         gg_21, hd_15, hf_31, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * fg_7[k]
                  + pa_y[k] * gg_21[k];

        t_62[k] = f_1 * gf_14[k]
                  + pb_y[k] * hf_31[k];

        t_63[k] = pb_z[k] * hf_31[k];

        t_64[k] = f_3 * gf_31[k]
                  + f_2 * hd_15[k]
                  + pb_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pb_x, pb_z, fg_28, gf_32, gg_51, \
                         hd_14, hd_15, hf_32, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * hd_14[k]
                  + pb_z[k] * hf_32[k];

        t_66[k] = f_3 * gf_32[k]
                  + pb_x[k] * hf_34[k];

        t_67[k] = f_2 * fg_28[k]
                  + pa_x[k] * gg_51[k];

        t_68[k] = pb_z[k] * hf_34[k];

        t_69[k] = f_2 * hd_15[k]
                  + pb_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_y, pb_z, gf_14, gf_19, gg_21, \
                         gg_22, hd_16, hf_36, hf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * gf_19[k]
                  + pb_y[k] * hf_36[k];

        t_71[k] = f_1 * hd_16[k]
                  + pb_z[k] * hf_36[k];

        t_72[k] = pa_z[k] * gg_21[k];

        t_73[k] = f_2 * gf_14[k]
                  + pb_z[k] * hf_37[k];

        t_74[k] = pa_z[k] * gg_22[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_x, pa_y, pa_z, pb_z, fg_12, fg_36, gf_17, \
                         gg_25, gg_31, gg_60, hf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * fg_12[k]
                  + pa_y[k] * gg_31[k];

        t_76[k] = pa_z[k] * gg_25[k];

        t_77[k] = f_2 * gf_17[k]
                  + pb_z[k] * hf_38[k];

        t_78[k] = f_2 * fg_36[k]
                  + pa_x[k] * gg_60[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_x, pa_y, pb_y, fg_38, gf_21, gf_23, \
                         gg_37, gg_39, gg_40, gg_62, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * gf_21[k]
                  + pb_y[k] * hf_39[k];

        t_80[k] = f_2 * fg_38[k]
                  + pa_x[k] * gg_62[k];

        t_81[k] = pa_y[k] * gg_37[k];

        t_82[k] = pa_y[k] * gg_39[k];

        t_83[k] = f_3 * gf_23[k]
                  + pa_y[k] * gg_40[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pa_y, pb_z, fg_40, fg_42, gf_20, gg_41, \
                         gg_67, gg_69, hf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_y[k] * gg_41[k];

        t_85[k] = f_2 * fg_40[k]
                  + pa_x[k] * gg_67[k];

        t_86[k] = f_3 * gf_20[k]
                  + pb_z[k] * hf_40[k];

        t_87[k] = f_2 * fg_42[k]
                  + pa_x[k] * gg_69[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, fg_10, gf_22, \
                         gf_28, gg_37, gg_46, hf_41, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * gf_28[k]
                  + pb_y[k] * hf_41[k];

        t_89[k] = pa_y[k] * gg_46[k];

        t_90[k] = f_3 * fg_10[k]
                  + pa_z[k] * gg_37[k];

        t_91[k] = pb_y[k] * hf_42[k];

        t_92[k] = f_1 * gf_22[k]
                  + pb_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_x, pb_y, gf_37, gf_38, hd_17, hd_20, \
                         hf_43, hf_44, hf_45, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * hd_17[k]
                  + pb_y[k] * hf_43[k];

        t_94[k] = pb_y[k] * hf_44[k];

        t_95[k] = f_3 * gf_37[k]
                  + f_2 * hd_20[k]
                  + pb_x[k] * hf_45[k];

        t_96[k] = f_3 * gf_38[k]
                  + pb_x[k] * hf_49[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pa_x, pb_y, fg_55, gg_78, hd_18, \
                         hd_19, hd_20, hf_46, hf_47, hf_48, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * hd_18[k]
                  + pb_y[k] * hf_46[k];

        t_98[k] = f_3 * hd_19[k]
                  + pb_y[k] * hf_47[k];

        t_99[k] = f_2 * hd_20[k]
                  + pb_y[k] * hf_48[k];

        t_100[k] = pb_y[k] * hf_49[k];

        t_101[k] = f_2 * fg_55[k]
                   + pa_x[k] * gg_78[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_x, pb_y, pb_z, gf_29, gf_39, \
                         gf_41, gg_79, gg_81, hd_21, hf_50, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * gf_39[k]
                   + pa_x[k] * gg_79[k];

        t_103[k] = f_4 * gf_29[k]
                   + pb_y[k] * hf_50[k];

        t_104[k] = pb_z[k] * hf_50[k];

        t_105[k] = f_3 * gf_41[k]
                   + pa_x[k] * gg_81[k];

        t_106[k] = f_2 * hd_21[k]
                   + pb_z[k] * hf_51[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, pa_x, pa_z, pb_x, gf_43, \
                         gg_47, gg_87, gg_89, gg_90, gg_91, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_2 * gf_43[k]
                   + pb_x[k] * hf_52[k];

        t_108[k] = pa_x[k] * gg_87[k];

        t_109[k] = pa_x[k] * gg_89[k];

        t_110[k] = pa_x[k] * gg_90[k];

        t_111[k] = pa_x[k] * gg_91[k];

        t_112[k] = pa_z[k] * gg_47[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_x, pa_z, pb_z, gf_29, gf_47, \
                         gg_48, gg_92, gg_95, gg_96, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_2 * gf_29[k]
                   + pb_z[k] * hf_53[k];

        t_114[k] = pa_z[k] * gg_48[k];

        t_115[k] = f_3 * gf_47[k]
                   + pa_x[k] * gg_92[k];

        t_116[k] = pa_x[k] * gg_95[k];

        t_117[k] = pa_x[k] * gg_96[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_x, pb_z, gf_33, gf_50, gf_51, \
                         gg_97, gg_98, gg_99, gg_100, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_x[k] * gg_97[k];

        t_119[k] = pa_x[k] * gg_98[k];

        t_120[k] = f_4 * gf_50[k]
                   + pa_x[k] * gg_99[k];

        t_121[k] = f_3 * gf_33[k]
                   + pb_z[k] * hf_54[k];

        t_122[k] = f_3 * gf_51[k]
                   + pa_x[k] * gg_100[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, pa_x, gf_52, gg_101, \
                         gg_104, gg_105, gg_106, gg_107, gg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_3 * gf_52[k]
                   + pa_x[k] * gg_101[k];

        t_124[k] = pa_x[k] * gg_104[k];

        t_125[k] = pa_x[k] * gg_105[k];

        t_126[k] = pa_x[k] * gg_106[k];

        t_127[k] = pa_x[k] * gg_107[k];

        t_128[k] = pa_x[k] * gg_108[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, t_134, pa_x, pa_y, gf_56, gg_71, \
                         gg_72, gg_73, gg_109, gg_111, gg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pa_y[k] * gg_71[k];

        t_130[k] = pa_y[k] * gg_72[k];

        t_131[k] = f_3 * gf_56[k]
                   + pa_x[k] * gg_109[k];

        t_132[k] = pa_y[k] * gg_73[k];

        t_133[k] = pa_x[k] * gg_111[k];

        t_134[k] = pa_x[k] * gg_112[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pa_x, pb_y, pb_z, gf_36, gf_60, \
                         gg_113, gg_114, gg_116, hf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pa_x[k] * gg_113[k];

        t_136[k] = pa_x[k] * gg_114[k];

        t_137[k] = f_4 * gf_60[k]
                   + pa_x[k] * gg_116[k];

        t_138[k] = pb_y[k] * hf_55[k];

        t_139[k] = f_4 * gf_36[k]
                   + pb_z[k] * hf_55[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pa_x, pb_x, pb_y, gf_63, gf_67, \
                         gg_121, gg_125, hd_22, hf_56, hf_57, hf_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_2 * hd_22[k]
                   + pb_y[k] * hf_56[k];

        t_141[k] = pb_y[k] * hf_57[k];

        t_142[k] = f_3 * gf_63[k]
                   + pa_x[k] * gg_121[k];

        t_143[k] = f_2 * gf_67[k]
                   + pb_x[k] * hf_58[k];

        t_144[k] = pa_x[k] * gg_125[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_x, pb_x, gg_126, gg_127, \
                         gg_129, hd_23, hd_24, hf_59, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_x[k] * gg_126[k];

        t_146[k] = pa_x[k] * gg_127[k];

        t_147[k] = pa_x[k] * gg_129[k];

        t_148[k] = f_1 * hd_23[k]
                   + pb_x[k] * hf_59[k];

        t_149[k] = f_3 * hd_24[k]
                   + pb_x[k] * hf_60[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pb_x, pb_y, gf_43, hd_25, \
                         hd_26, hf_61, hf_62, hf_63, hf_65, hf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_2 * hd_25[k]
                   + pb_x[k] * hf_61[k];

        t_151[k] = f_2 * hd_26[k]
                   + pb_x[k] * hf_62[k];

        t_152[k] = pb_x[k] * hf_63[k];

        t_153[k] = pb_x[k] * hf_65[k];

        t_154[k] = pb_x[k] * hf_66[k];

        t_155[k] = f_0 * gf_43[k]
                   + f_1 * hd_25[k]
                   + pb_y[k] * hf_63[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pb_x, pb_y, pb_z, gf_46, hd_25, \
                         hd_26, hd_27, hf_63, hf_64, hf_66, hf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_z[k] * hf_63[k];

        t_157[k] = f_2 * hd_25[k]
                   + pb_z[k] * hf_64[k];

        t_158[k] = f_0 * gf_46[k]
                   + pb_y[k] * hf_66[k];

        t_159[k] = f_1 * hd_26[k]
                   + pb_z[k] * hf_66[k];

        t_160[k] = f_3 * hd_27[k]
                   + pb_x[k] * hf_67[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, pa_z, pb_x, gg_87, hd_29, \
                         hd_30, hf_68, hf_69, hf_71, hf_72, hf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_2 * hd_29[k]
                   + pb_x[k] * hf_68[k];

        t_162[k] = f_2 * hd_30[k]
                   + pb_x[k] * hf_69[k];

        t_163[k] = pb_x[k] * hf_71[k];

        t_164[k] = pb_x[k] * hf_72[k];

        t_165[k] = pb_x[k] * hf_73[k];

        t_166[k] = pa_z[k] * gg_87[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pa_z, pb_y, pb_z, fg_38, gf_43, \
                         gf_44, gf_49, gg_89, gg_98, hf_70, hf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_2 * gf_43[k]
                   + pb_z[k] * hf_70[k];

        t_168[k] = f_3 * gf_44[k]
                   + pa_z[k] * gg_89[k];

        t_169[k] = f_4 * gf_49[k]
                   + pb_y[k] * hf_73[k];

        t_170[k] = f_1 * fg_38[k]
                   + pa_y[k] * gg_98[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, t_175, pb_x, hd_31, hd_32, hd_33, hd_34, \
                         hd_35, hf_74, hf_75, hf_76, hf_77, hf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_1 * hd_31[k]
                   + pb_x[k] * hf_74[k];

        t_172[k] = f_3 * hd_32[k]
                   + pb_x[k] * hf_75[k];

        t_173[k] = f_3 * hd_33[k]
                   + pb_x[k] * hf_76[k];

        t_174[k] = f_2 * hd_34[k]
                   + pb_x[k] * hf_77[k];

        t_175[k] = f_2 * hd_35[k]
                   + pb_x[k] * hf_78[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, t_181, pa_z, pb_x, fg_28, gg_94, \
                         hd_36, hf_79, hf_80, hf_81, hf_82, hf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_2 * hd_36[k]
                   + pb_x[k] * hf_79[k];

        t_177[k] = pb_x[k] * hf_80[k];

        t_178[k] = pb_x[k] * hf_81[k];

        t_179[k] = pb_x[k] * hf_82[k];

        t_180[k] = pb_x[k] * hf_83[k];

        t_181[k] = f_2 * fg_28[k]
                   + pa_z[k] * gg_94[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pa_y, pb_y, pb_z, fg_44, gf_48, gf_54, \
                         gf_55, gg_108, hd_36, hf_80, hf_82, hf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_3 * gf_48[k]
                   + pb_z[k] * hf_80[k];

        t_183[k] = f_1 * gf_54[k]
                   + f_2 * hd_36[k]
                   + pb_y[k] * hf_82[k];

        t_184[k] = f_1 * gf_55[k]
                   + pb_y[k] * hf_83[k];

        t_185[k] = f_3 * fg_44[k]
                   + pa_y[k] * gg_108[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pb_x, hd_37, hd_38, hd_39, hd_40, \
                         hd_41, hf_84, hf_85, hf_86, hf_87, hf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_1 * hd_37[k]
                   + pb_x[k] * hf_84[k];

        t_187[k] = f_3 * hd_38[k]
                   + pb_x[k] * hf_85[k];

        t_188[k] = f_3 * hd_39[k]
                   + pb_x[k] * hf_86[k];

        t_189[k] = f_2 * hd_40[k]
                   + pb_x[k] * hf_87[k];

        t_190[k] = f_2 * hd_41[k]
                   + pb_x[k] * hf_88[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, pa_z, pb_x, fg_34, gg_104, \
                         hd_42, hf_89, hf_90, hf_91, hf_92, hf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_2 * hd_42[k]
                   + pb_x[k] * hf_89[k];

        t_192[k] = pb_x[k] * hf_90[k];

        t_193[k] = pb_x[k] * hf_91[k];

        t_194[k] = pb_x[k] * hf_92[k];

        t_195[k] = pb_x[k] * hf_93[k];

        t_196[k] = f_3 * fg_34[k]
                   + pa_z[k] * gg_104[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_y, pb_y, pb_z, fg_55, gf_53, gf_58, \
                         gf_59, gg_115, hd_42, hf_90, hf_92, hf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_1 * gf_53[k]
                   + pb_z[k] * hf_90[k];

        t_198[k] = f_3 * gf_58[k]
                   + f_2 * hd_42[k]
                   + pb_y[k] * hf_92[k];

        t_199[k] = f_3 * gf_59[k]
                   + pb_y[k] * hf_93[k];

        t_200[k] = f_2 * fg_55[k]
                   + pa_y[k] * gg_115[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pb_x, hd_43, hd_44, hd_45, \
                         hf_94, hf_95, hf_96, hf_97, hf_98, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_3 * hd_43[k]
                   + pb_x[k] * hf_94[k];

        t_202[k] = f_2 * hd_44[k]
                   + pb_x[k] * hf_95[k];

        t_203[k] = f_2 * hd_45[k]
                   + pb_x[k] * hf_96[k];

        t_204[k] = pb_x[k] * hf_97[k];

        t_205[k] = pb_x[k] * hf_98[k];

        t_206[k] = pb_x[k] * hf_99[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pa_y, pb_y, pb_z, gf_57, gf_64, gf_66, \
                         gf_67, gg_125, gg_127, hf_97, hf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_4 * gf_64[k]
                   + pa_y[k] * gg_125[k];

        t_208[k] = f_4 * gf_57[k]
                   + pb_z[k] * hf_97[k];

        t_209[k] = f_3 * gf_66[k]
                   + pa_y[k] * gg_127[k];

        t_210[k] = f_2 * gf_67[k]
                   + pb_y[k] * hf_100[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pa_y, pb_x, gg_129, hd_47, hd_48, \
                         hd_49, hd_51, hf_101, hf_102, hf_103, hf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = pa_y[k] * gg_129[k];

        t_212[k] = f_1 * hd_47[k]
                   + pb_x[k] * hf_101[k];

        t_213[k] = f_3 * hd_48[k]
                   + pb_x[k] * hf_102[k];

        t_214[k] = f_2 * hd_49[k]
                   + pb_x[k] * hf_103[k];

        t_215[k] = f_2 * hd_51[k]
                   + pb_x[k] * hf_104[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, t_222, pb_x, pb_y, hd_49, \
                         hd_50, hd_51, hf_105, hf_106, hf_107, hf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * hf_105[k];

        t_217[k] = pb_x[k] * hf_106[k];

        t_218[k] = pb_x[k] * hf_108[k];

        t_219[k] = f_1 * hd_49[k]
                   + pb_y[k] * hf_105[k];

        t_220[k] = f_3 * hd_50[k]
                   + pb_y[k] * hf_106[k];

        t_221[k] = f_2 * hd_51[k]
                   + pb_y[k] * hf_107[k];

        t_222[k] = pb_y[k] * hf_108[k];
    }

#pragma omp simd aligned(t_223, pb_z, gf_67, hd_51, hf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_0 * gf_67[k]
                   + f_1 * hd_51[k]
                   + pb_z[k] * hf_108[k];
    }
}

auto
compute_prim_hg_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fg, const size_t gf, const size_t gg,
                          const size_t hd, const size_t hf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_49 = buffer.data(fg + 49);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_97 = buffer.data(gg + 97);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gf_0, hd_0, hd_1, \
                         hf_0, hf_1, hf_2, hf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_2 * hd_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_2 * hd_0[k]
                 + pb_z[k] * hf_2[k];

        t_5[k] = f_1 * hd_1[k]
                 + pb_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_y, pb_z, gf_1, gg_0, gg_3, hd_2, \
                         hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hd_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd_2[k]
                 + pb_z[k] * hf_5[k];

        t_9[k] = pa_y[k] * gg_0[k];

        t_10[k] = f_3 * gf_1[k]
                  + pa_y[k] * gg_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, pb_z, fg_8, gg_0, \
                         gg_8, gg_10, hd_3, hf_6, hf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * fg_8[k]
                  + pa_x[k] * gg_10[k];

        t_12[k] = pb_z[k] * hf_6[k];

        t_13[k] = f_2 * hd_3[k]
                  + pb_z[k] * hf_7[k];

        t_14[k] = pa_y[k] * gg_8[k];

        t_15[k] = pa_z[k] * gg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_y, pb_z, gf_0, gf_2, gg_4, hd_5, \
                         hf_8, hf_9, hf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * gf_0[k]
                  + pb_z[k] * hf_8[k];

        t_17[k] = pb_y[k] * hf_9[k];

        t_18[k] = f_3 * gf_2[k]
                  + pa_z[k] * gg_4[k];

        t_19[k] = f_3 * hd_5[k]
                  + pb_y[k] * hf_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, fg_0, fg_11, gg_9, gg_19, \
                         hd_6, hf_11, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * hd_6[k]
                  + pb_y[k] * hf_11[k];

        t_21[k] = pb_y[k] * hf_12[k];

        t_22[k] = f_1 * fg_11[k]
                  + pa_x[k] * gg_19[k];

        t_23[k] = f_2 * fg_0[k]
                  + pa_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_z, gf_15, gf_16, hd_7, hd_8, hf_13, \
                         hf_14, hf_15, hf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_z[k] * hf_13[k];

        t_25[k] = f_1 * gf_15[k]
                  + f_2 * hd_8[k]
                  + pb_x[k] * hf_15[k];

        t_26[k] = f_2 * hd_7[k]
                  + pb_z[k] * hf_14[k];

        t_27[k] = f_1 * gf_16[k]
                  + pb_x[k] * hf_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, fg_15, gg_10, gg_25, \
                         hd_8, hd_9, hf_16, hf_17, hf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * fg_15[k]
                  + pa_x[k] * gg_25[k];

        t_29[k] = pb_z[k] * hf_16[k];

        t_30[k] = f_2 * hd_8[k]
                  + pb_z[k] * hf_17[k];

        t_31[k] = f_1 * hd_9[k]
                  + pb_z[k] * hf_18[k];

        t_32[k] = pa_z[k] * gg_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, pb_y, pb_z, fg_0, gf_9, \
                         gg_14, gg_19, hd_10, hf_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * gg_19[k];

        t_34[k] = f_2 * fg_0[k]
                  + pa_z[k] * gg_14[k];

        t_35[k] = pb_y[k] * hf_19[k];

        t_36[k] = f_3 * gf_9[k]
                  + pb_z[k] * hf_19[k];

        t_37[k] = f_2 * hd_10[k]
                  + pb_y[k] * hf_20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, gf_21, gf_25, hd_11, hd_13, \
                         hf_21, hf_22, hf_23, hf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * hf_21[k];

        t_39[k] = f_1 * gf_21[k]
                  + f_2 * hd_13[k]
                  + pb_x[k] * hf_22[k];

        t_40[k] = f_1 * gf_25[k]
                  + pb_x[k] * hf_26[k];

        t_41[k] = f_1 * hd_11[k]
                  + pb_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_y, fg_19, gg_40, hd_12, hd_13, \
                         hf_24, hf_25, hf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hd_12[k]
                  + pb_y[k] * hf_24[k];

        t_43[k] = f_2 * hd_13[k]
                  + pb_y[k] * hf_25[k];

        t_44[k] = pb_y[k] * hf_26[k];

        t_45[k] = f_3 * fg_19[k]
                  + pa_x[k] * gg_40[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pb_x, pb_z, fg_7, gf_28, gg_20, hd_14, \
                         hd_15, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * fg_7[k]
                  + pa_y[k] * gg_20[k];

        t_47[k] = pb_z[k] * hf_27[k];

        t_48[k] = f_3 * gf_28[k]
                  + f_2 * hd_15[k]
                  + pb_x[k] * hf_29[k];

        t_49[k] = f_2 * hd_14[k]
                  + pb_z[k] * hf_28[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_x, pb_x, pb_z, fg_25, gf_29, gg_46, \
                         hd_15, hd_16, hf_30, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * gf_29[k]
                  + pb_x[k] * hf_30[k];

        t_51[k] = f_2 * fg_25[k]
                  + pa_x[k] * gg_46[k];

        t_52[k] = pb_z[k] * hf_30[k];

        t_53[k] = f_2 * hd_15[k]
                  + pb_z[k] * hf_31[k];

        t_54[k] = f_1 * hd_16[k]
                  + pb_z[k] * hf_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pa_y, pa_z, fg_32, fg_35, gg_20, \
                         gg_25, gg_40, gg_48, gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * gg_20[k];

        t_56[k] = pa_z[k] * gg_25[k];

        t_57[k] = f_2 * fg_32[k]
                  + pa_x[k] * gg_48[k];

        t_58[k] = f_2 * fg_35[k]
                  + pa_x[k] * gg_49[k];

        t_59[k] = pa_y[k] * gg_40[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_z, pb_y, pb_z, fg_10, gf_19, gg_31, \
                         hd_17, hf_33, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * fg_10[k]
                  + pa_z[k] * gg_31[k];

        t_61[k] = pb_y[k] * hf_33[k];

        t_62[k] = f_1 * gf_19[k]
                  + pb_z[k] * hf_33[k];

        t_63[k] = f_2 * hd_17[k]
                  + pb_y[k] * hf_34[k];

        t_64[k] = pb_y[k] * hf_35[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_x, pb_y, gf_31, gf_32, hd_18, hd_19, \
                         hd_20, hf_36, hf_37, hf_38, hf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * gf_31[k]
                  + f_2 * hd_20[k]
                  + pb_x[k] * hf_36[k];

        t_66[k] = f_3 * gf_32[k]
                  + pb_x[k] * hf_40[k];

        t_67[k] = f_1 * hd_18[k]
                  + pb_y[k] * hf_37[k];

        t_68[k] = f_3 * hd_19[k]
                  + pb_y[k] * hf_38[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_x, pb_y, pb_z, fg_49, gf_33, gg_53, \
                         gg_54, hd_20, hf_39, hf_40, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * hd_20[k]
                  + pb_y[k] * hf_39[k];

        t_70[k] = pb_y[k] * hf_40[k];

        t_71[k] = f_2 * fg_49[k]
                  + pa_x[k] * gg_53[k];

        t_72[k] = f_4 * gf_33[k]
                  + pa_x[k] * gg_54[k];

        t_73[k] = pb_z[k] * hf_41[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pa_x, pa_z, pb_z, gf_35, gf_44, gg_41, \
                         gg_56, gg_61, gg_71, hd_21, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * gf_35[k]
                  + pa_x[k] * gg_56[k];

        t_75[k] = f_2 * hd_21[k]
                  + pb_z[k] * hf_42[k];

        t_76[k] = pa_x[k] * gg_61[k];

        t_77[k] = pa_z[k] * gg_41[k];

        t_78[k] = f_4 * gf_44[k]
                  + pa_x[k] * gg_71[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_x, pb_y, pb_z, gf_30, gf_54, gg_86, \
                         hd_22, hf_43, hf_44, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_4 * gf_54[k]
                  + pa_x[k] * gg_86[k];

        t_80[k] = pb_y[k] * hf_43[k];

        t_81[k] = f_4 * gf_30[k]
                  + pb_z[k] * hf_43[k];

        t_82[k] = f_2 * hd_22[k]
                  + pb_y[k] * hf_44[k];

        t_83[k] = pb_y[k] * hf_45[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_x, pb_x, gf_57, gg_89, gg_97, hd_23, \
                         hd_24, hd_25, hf_46, hf_47, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_3 * gf_57[k]
                  + pa_x[k] * gg_89[k];

        t_85[k] = pa_x[k] * gg_97[k];

        t_86[k] = f_1 * hd_23[k]
                  + pb_x[k] * hf_46[k];

        t_87[k] = f_3 * hd_24[k]
                  + pb_x[k] * hf_47[k];

        t_88[k] = f_2 * hd_25[k]
                  + pb_x[k] * hf_48[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, t_94, pb_x, pb_y, pb_z, gf_37, hd_25, \
                         hd_26, hf_49, hf_50, hf_52, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_2 * hd_26[k]
                  + pb_x[k] * hf_49[k];

        t_90[k] = pb_x[k] * hf_50[k];

        t_91[k] = pb_x[k] * hf_52[k];

        t_92[k] = pb_x[k] * hf_53[k];

        t_93[k] = f_0 * gf_37[k]
                  + f_1 * hd_25[k]
                  + pb_y[k] * hf_50[k];

        t_94[k] = pb_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_x, pb_y, pb_z, gf_40, hd_25, hd_26, hd_27, \
                         hf_51, hf_53, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_2 * hd_25[k]
                  + pb_z[k] * hf_51[k];

        t_96[k] = f_0 * gf_40[k]
                  + pb_y[k] * hf_53[k];

        t_97[k] = f_1 * hd_26[k]
                  + pb_z[k] * hf_53[k];

        t_98[k] = f_3 * hd_27[k]
                  + pb_x[k] * hf_54[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, pa_z, pb_x, gg_61, hd_29, \
                         hd_30, hf_55, hf_56, hf_58, hf_59, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_2 * hd_29[k]
                  + pb_x[k] * hf_55[k];

        t_100[k] = f_2 * hd_30[k]
                   + pb_x[k] * hf_56[k];

        t_101[k] = pb_x[k] * hf_58[k];

        t_102[k] = pb_x[k] * hf_59[k];

        t_103[k] = pb_x[k] * hf_60[k];

        t_104[k] = pa_z[k] * gg_61[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_y, pa_z, pb_y, pb_z, fg_32, gf_37, \
                         gf_38, gf_43, gg_63, gg_70, hf_57, hf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_2 * gf_37[k]
                   + pb_z[k] * hf_57[k];

        t_106[k] = f_3 * gf_38[k]
                   + pa_z[k] * gg_63[k];

        t_107[k] = f_4 * gf_43[k]
                   + pb_y[k] * hf_60[k];

        t_108[k] = f_1 * fg_32[k]
                   + pa_y[k] * gg_70[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, hd_31, hd_32, hd_33, hd_34, \
                         hd_35, hf_61, hf_62, hf_63, hf_64, hf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_1 * hd_31[k]
                   + pb_x[k] * hf_61[k];

        t_110[k] = f_3 * hd_32[k]
                   + pb_x[k] * hf_62[k];

        t_111[k] = f_3 * hd_33[k]
                   + pb_x[k] * hf_63[k];

        t_112[k] = f_2 * hd_34[k]
                   + pb_x[k] * hf_64[k];

        t_113[k] = f_2 * hd_35[k]
                   + pb_x[k] * hf_65[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pa_z, pb_x, fg_25, gg_68, \
                         hd_36, hf_66, hf_67, hf_68, hf_69, hf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * hd_36[k]
                   + pb_x[k] * hf_66[k];

        t_115[k] = pb_x[k] * hf_67[k];

        t_116[k] = pb_x[k] * hf_68[k];

        t_117[k] = pb_x[k] * hf_69[k];

        t_118[k] = pb_x[k] * hf_70[k];

        t_119[k] = f_2 * fg_25[k]
                   + pa_z[k] * gg_68[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pb_y, pb_z, fg_38, gf_42, gf_48, \
                         gf_49, gg_79, hd_36, hf_67, hf_69, hf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_3 * gf_42[k]
                   + pb_z[k] * hf_67[k];

        t_121[k] = f_1 * gf_48[k]
                   + f_2 * hd_36[k]
                   + pb_y[k] * hf_69[k];

        t_122[k] = f_1 * gf_49[k]
                   + pb_y[k] * hf_70[k];

        t_123[k] = f_3 * fg_38[k]
                   + pa_y[k] * gg_79[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pb_x, hd_37, hd_38, hd_39, hd_40, \
                         hd_41, hf_71, hf_72, hf_73, hf_74, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * hd_37[k]
                   + pb_x[k] * hf_71[k];

        t_125[k] = f_3 * hd_38[k]
                   + pb_x[k] * hf_72[k];

        t_126[k] = f_3 * hd_39[k]
                   + pb_x[k] * hf_73[k];

        t_127[k] = f_2 * hd_40[k]
                   + pb_x[k] * hf_74[k];

        t_128[k] = f_2 * hd_41[k]
                   + pb_x[k] * hf_75[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, t_134, pa_z, pb_x, fg_31, gg_76, \
                         hd_42, hf_76, hf_77, hf_78, hf_79, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_2 * hd_42[k]
                   + pb_x[k] * hf_76[k];

        t_130[k] = pb_x[k] * hf_77[k];

        t_131[k] = pb_x[k] * hf_78[k];

        t_132[k] = pb_x[k] * hf_79[k];

        t_133[k] = pb_x[k] * hf_80[k];

        t_134[k] = f_3 * fg_31[k]
                   + pa_z[k] * gg_76[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_y, pb_y, pb_z, fg_49, gf_47, gf_52, \
                         gf_53, gg_85, hd_42, hf_77, hf_79, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * gf_47[k]
                   + pb_z[k] * hf_77[k];

        t_136[k] = f_3 * gf_52[k]
                   + f_2 * hd_42[k]
                   + pb_y[k] * hf_79[k];

        t_137[k] = f_3 * gf_53[k]
                   + pb_y[k] * hf_80[k];

        t_138[k] = f_2 * fg_49[k]
                   + pa_y[k] * gg_85[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, pb_x, hd_43, hd_44, hd_45, \
                         hf_81, hf_82, hf_83, hf_84, hf_85, hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * hd_43[k]
                   + pb_x[k] * hf_81[k];

        t_140[k] = f_2 * hd_44[k]
                   + pb_x[k] * hf_82[k];

        t_141[k] = f_2 * hd_45[k]
                   + pb_x[k] * hf_83[k];

        t_142[k] = pb_x[k] * hf_84[k];

        t_143[k] = pb_x[k] * hf_85[k];

        t_144[k] = pb_x[k] * hf_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_y, pb_y, pb_z, gf_51, gf_58, gf_60, \
                         gf_61, gg_93, gg_95, hf_84, hf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_4 * gf_58[k]
                   + pa_y[k] * gg_93[k];

        t_146[k] = f_4 * gf_51[k]
                   + pb_z[k] * hf_84[k];

        t_147[k] = f_3 * gf_60[k]
                   + pa_y[k] * gg_95[k];

        t_148[k] = f_2 * gf_61[k]
                   + pb_y[k] * hf_87[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pa_y, pb_x, gg_97, hd_47, hd_48, \
                         hd_49, hd_51, hf_88, hf_89, hf_90, hf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_y[k] * gg_97[k];

        t_150[k] = f_1 * hd_47[k]
                   + pb_x[k] * hf_88[k];

        t_151[k] = f_3 * hd_48[k]
                   + pb_x[k] * hf_89[k];

        t_152[k] = f_2 * hd_49[k]
                   + pb_x[k] * hf_90[k];

        t_153[k] = f_2 * hd_51[k]
                   + pb_x[k] * hf_91[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, t_159, t_160, pb_x, pb_y, hd_49, \
                         hd_50, hd_51, hf_92, hf_93, hf_94, hf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_x[k] * hf_92[k];

        t_155[k] = pb_x[k] * hf_93[k];

        t_156[k] = pb_x[k] * hf_95[k];

        t_157[k] = f_1 * hd_49[k]
                   + pb_y[k] * hf_92[k];

        t_158[k] = f_3 * hd_50[k]
                   + pb_y[k] * hf_93[k];

        t_159[k] = f_2 * hd_51[k]
                   + pb_y[k] * hf_94[k];

        t_160[k] = pb_y[k] * hf_95[k];
    }

#pragma omp simd aligned(t_161, pb_z, gf_61, hd_51, hf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * gf_61[k]
                   + f_1 * hd_51[k]
                   + pb_z[k] * hf_95[k];
    }
}

auto
compute_prim_hg_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fg, const size_t gf, const size_t gg,
                          const size_t hd, const size_t hf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_20 = buffer.data(fg + 20);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);

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
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gf_0, gf_3, hd_0, hf_0, hf_1, \
                         hf_2, hf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = f_2 * hd_0[k]
                 + pb_y[k] * hf_1[k];

        t_2[k] = f_2 * hd_0[k]
                 + pb_z[k] * hf_2[k];

        t_3[k] = f_0 * gf_3[k]
                 + pb_x[k] * hf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, gf_4, gg_0, hd_1, hd_2, \
                         hf_3, hf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * gf_4[k]
                 + pb_x[k] * hf_4[k];

        t_5[k] = f_1 * hd_1[k]
                 + pb_y[k] * hf_3[k];

        t_6[k] = f_1 * hd_2[k]
                 + pb_z[k] * hf_4[k];

        t_7[k] = pa_y[k] * gg_0[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, pb_x, pb_y, fg_3, gf_0, gf_1, gf_6, \
                         gg_1, gg_4, hf_5, hf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * gf_0[k]
                 + pb_y[k] * hf_5[k];

        t_9[k] = f_3 * gf_1[k]
                 + pa_y[k] * gg_1[k];

        t_10[k] = f_4 * gf_6[k]
                  + pb_x[k] * hf_6[k];

        t_11[k] = f_1 * fg_3[k]
                  + pa_x[k] * gg_4[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_z, gf_0, gf_2, gf_9, gg_0, \
                         gg_2, hf_9, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * gg_0[k];

        t_13[k] = f_2 * gf_0[k]
                  + pb_z[k] * hf_9[k];

        t_14[k] = f_3 * gf_2[k]
                  + pa_z[k] * gg_2[k];

        t_15[k] = f_4 * gf_9[k]
                  + pb_x[k] * hf_12[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_y, fg_0, fg_6, gf_5, gg_3, gg_7, \
                         hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * fg_6[k]
                  + pa_x[k] * gg_7[k];

        t_17[k] = f_2 * fg_0[k]
                  + pa_y[k] * gg_3[k];

        t_18[k] = f_3 * gf_5[k]
                  + pb_y[k] * hf_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_y, pb_x, fg_7, gf_11, gf_12, gg_6, \
                         gg_10, hd_10, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * gf_11[k]
                  + f_2 * hd_10[k]
                  + pb_x[k] * hf_14[k];

        t_20[k] = f_1 * gf_12[k]
                  + pb_x[k] * hf_15[k];

        t_21[k] = f_3 * fg_7[k]
                  + pa_x[k] * gg_10[k];

        t_22[k] = pa_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_z, pb_y, pb_z, fg_0, fg_8, gf_7, \
                         gg_5, gg_12, hd_14, hf_23, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * fg_8[k]
                  + pa_x[k] * gg_12[k];

        t_24[k] = f_2 * fg_0[k]
                  + pa_z[k] * gg_5[k];

        t_25[k] = f_3 * gf_7[k]
                  + pb_z[k] * hf_23[k];

        t_26[k] = f_2 * hd_14[k]
                  + pb_y[k] * hf_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pb_x, fg_2, fg_9, gf_15, gf_16, \
                         gg_8, gg_15, hd_17, hf_26, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * gf_15[k]
                  + f_2 * hd_17[k]
                  + pb_x[k] * hf_26[k];

        t_28[k] = f_1 * gf_16[k]
                  + pb_x[k] * hf_29[k];

        t_29[k] = f_3 * fg_9[k]
                  + pa_x[k] * gg_15[k];

        t_30[k] = f_3 * fg_2[k]
                  + pa_y[k] * gg_8[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_x, pb_y, fg_10, gf_10, gf_18, gf_19, \
                         gg_16, hd_19, hf_30, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * gf_10[k]
                  + pb_y[k] * hf_30[k];

        t_32[k] = f_3 * gf_18[k]
                  + f_2 * hd_19[k]
                  + pb_x[k] * hf_31[k];

        t_33[k] = f_3 * gf_19[k]
                  + pb_x[k] * hf_32[k];

        t_34[k] = f_2 * fg_10[k]
                  + pa_x[k] * gg_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pa_y, pa_z, fg_5, fg_13, fg_14, \
                         gg_9, gg_11, gg_13, gg_17, gg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * gg_9[k];

        t_36[k] = f_2 * fg_5[k]
                  + pa_y[k] * gg_11[k];

        t_37[k] = f_2 * fg_13[k]
                  + pa_x[k] * gg_17[k];

        t_38[k] = f_2 * fg_14[k]
                  + pa_x[k] * gg_18[k];

        t_39[k] = pa_y[k] * gg_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pa_z, fg_4, fg_15, fg_16, gg_13, \
                         gg_14, gg_19, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * gg_14[k];

        t_41[k] = f_2 * fg_15[k]
                  + pa_x[k] * gg_19[k];

        t_42[k] = f_2 * fg_16[k]
                  + pa_x[k] * gg_20[k];

        t_43[k] = f_3 * fg_4[k]
                  + pa_z[k] * gg_13[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, gf_13, gf_24, gf_25, hd_27, \
                         hd_30, hf_47, hf_48, hf_50, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * gf_13[k]
                  + pb_z[k] * hf_47[k];

        t_45[k] = f_2 * hd_27[k]
                  + pb_y[k] * hf_48[k];

        t_46[k] = f_3 * gf_24[k]
                  + f_2 * hd_30[k]
                  + pb_x[k] * hf_50[k];

        t_47[k] = f_3 * gf_25[k]
                  + pb_x[k] * hf_53[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_y, fg_20, gf_17, gf_26, gf_27, \
                         gg_21, gg_22, gg_23, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * fg_20[k]
                  + pa_x[k] * gg_21[k];

        t_49[k] = f_4 * gf_26[k]
                  + pa_x[k] * gg_22[k];

        t_50[k] = f_4 * gf_17[k]
                  + pb_y[k] * hf_54[k];

        t_51[k] = f_3 * gf_27[k]
                  + pa_x[k] * gg_23[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, pa_x, pb_x, gf_28, gg_24, gg_27, \
                         gg_28, gg_29, gg_30, hf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * gf_28[k]
                  + pb_x[k] * hf_56[k];

        t_53[k] = pa_x[k] * gg_24[k];

        t_54[k] = pa_x[k] * gg_27[k];

        t_55[k] = pa_x[k] * gg_28[k];

        t_56[k] = pa_x[k] * gg_29[k];

        t_57[k] = pa_x[k] * gg_30[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pb_z, gf_22, gf_45, gg_31, gg_32, \
                         gg_33, gg_35, hf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_x[k] * gg_31[k];

        t_59[k] = pa_x[k] * gg_32[k];

        t_60[k] = pa_x[k] * gg_33[k];

        t_61[k] = f_4 * gf_45[k]
                  + pa_x[k] * gg_35[k];

        t_62[k] = f_4 * gf_22[k]
                  + pb_z[k] * hf_75[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_x, pb_x, gf_48, gf_52, gg_37, gg_40, \
                         hd_41, hd_42, hf_80, hf_81, hf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * gf_48[k]
                  + pa_x[k] * gg_37[k];

        t_64[k] = f_2 * gf_52[k]
                  + pb_x[k] * hf_80[k];

        t_65[k] = pa_x[k] * gg_40[k];

        t_66[k] = f_1 * hd_41[k]
                  + pb_x[k] * hf_81[k];

        t_67[k] = f_3 * hd_42[k]
                  + pb_x[k] * hf_82[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, gf_28, hd_43, hd_44, \
                         hf_82, hf_83, hf_84, hf_85, hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * hd_43[k]
                  + pb_x[k] * hf_83[k];

        t_69[k] = pb_z[k] * hf_82[k];

        t_70[k] = f_2 * hd_44[k]
                  + pb_x[k] * hf_84[k];

        t_71[k] = f_0 * gf_28[k]
                  + f_1 * hd_43[k]
                  + pb_y[k] * hf_85[k];

        t_72[k] = f_2 * hd_43[k]
                  + pb_z[k] * hf_86[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pb_x, pb_y, pb_z, gf_31, gg_24, hd_44, \
                         hd_46, hf_88, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * gf_31[k]
                  + pb_y[k] * hf_88[k];

        t_74[k] = f_1 * hd_44[k]
                  + pb_z[k] * hf_88[k];

        t_75[k] = f_2 * hd_46[k]
                  + pb_x[k] * hf_89[k];

        t_76[k] = pa_z[k] * gg_24[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pa_z, pb_y, pb_z, fg_14, gf_28, gf_29, \
                         gf_35, gg_25, gg_28, hf_90, hf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_2 * gf_28[k]
                  + pb_z[k] * hf_90[k];

        t_78[k] = f_3 * gf_29[k]
                  + pa_z[k] * gg_25[k];

        t_79[k] = f_4 * gf_35[k]
                  + pb_y[k] * hf_93[k];

        t_80[k] = f_1 * fg_14[k]
                  + pa_y[k] * gg_28[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_z, pb_x, fg_10, gg_26, hd_47, hd_48, \
                         hd_49, hf_94, hf_95, hf_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_1 * hd_47[k]
                  + pb_x[k] * hf_94[k];

        t_82[k] = f_2 * hd_48[k]
                  + pb_x[k] * hf_95[k];

        t_83[k] = f_2 * hd_49[k]
                  + pb_x[k] * hf_96[k];

        t_84[k] = f_2 * fg_10[k]
                  + pa_z[k] * gg_26[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pb_y, pb_z, fg_17, gf_32, gf_39, gf_40, \
                         gg_31, hd_49, hf_97, hf_99, hf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_3 * gf_32[k]
                  + pb_z[k] * hf_97[k];

        t_86[k] = f_1 * gf_39[k]
                  + f_2 * hd_49[k]
                  + pb_y[k] * hf_99[k];

        t_87[k] = f_1 * gf_40[k]
                  + pb_y[k] * hf_100[k];

        t_88[k] = f_3 * fg_17[k]
                  + pa_y[k] * gg_31[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_z, pb_x, fg_12, gg_29, hd_50, hd_51, \
                         hd_52, hf_101, hf_102, hf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_1 * hd_50[k]
                  + pb_x[k] * hf_101[k];

        t_90[k] = f_2 * hd_51[k]
                  + pb_x[k] * hf_102[k];

        t_91[k] = f_2 * hd_52[k]
                  + pb_x[k] * hf_103[k];

        t_92[k] = f_3 * fg_12[k]
                  + pa_z[k] * gg_29[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_y, pb_z, fg_20, gf_37, gf_43, gf_44, \
                         gg_34, hd_52, hf_104, hf_106, hf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_1 * gf_37[k]
                  + pb_z[k] * hf_104[k];

        t_94[k] = f_3 * gf_43[k]
                  + f_2 * hd_52[k]
                  + pb_y[k] * hf_106[k];

        t_95[k] = f_3 * gf_44[k]
                  + pb_y[k] * hf_107[k];

        t_96[k] = f_2 * fg_20[k]
                  + pa_y[k] * gg_34[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pb_x, pb_z, gf_41, gf_49, gf_51, \
                         gg_38, gg_39, hd_53, hf_108, hf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * hd_53[k]
                  + pb_x[k] * hf_108[k];

        t_98[k] = f_4 * gf_49[k]
                  + pa_y[k] * gg_38[k];

        t_99[k] = f_4 * gf_41[k]
                  + pb_z[k] * hf_109[k];

        t_100[k] = f_3 * gf_51[k]
                   + pa_y[k] * gg_39[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pa_y, pb_x, pb_y, gf_52, gg_40, \
                         hd_55, hd_56, hf_112, hf_113, hf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * gf_52[k]
                   + pb_y[k] * hf_112[k];

        t_102[k] = pa_y[k] * gg_40[k];

        t_103[k] = f_1 * hd_55[k]
                   + pb_x[k] * hf_113[k];

        t_104[k] = pb_y[k] * hf_113[k];

        t_105[k] = f_3 * hd_56[k]
                   + pb_x[k] * hf_115[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pb_x, pb_y, hd_57, hd_58, hd_59, \
                         hf_115, hf_116, hf_117, hf_118, hf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_2 * hd_57[k]
                   + pb_x[k] * hf_116[k];

        t_107[k] = pb_y[k] * hf_115[k];

        t_108[k] = f_2 * hd_59[k]
                   + pb_x[k] * hf_117[k];

        t_109[k] = f_1 * hd_57[k]
                   + pb_y[k] * hf_118[k];

        t_110[k] = f_3 * hd_58[k]
                   + pb_y[k] * hf_119[k];
    }

#pragma omp simd aligned(t_111, t_112, pb_y, pb_z, gf_52, hd_59, hf_120, \
                         hf_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_2 * hd_59[k]
                   + pb_y[k] * hf_120[k];

        t_112[k] = f_0 * gf_52[k]
                   + f_1 * hd_59[k]
                   + pb_z[k] * hf_121[k];
    }
}

auto
compute_prim_hg_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fg, const size_t gf, const size_t gg,
                          const size_t hd, const size_t hf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_34 = buffer.data(fg + 34);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_4 = buffer.data(gf + 4);
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
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
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
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
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
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
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
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_43 = buffer.data(gg + 43);
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
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_78 = buffer.data(gg + 78);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);

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
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_58 = buffer.data(hf + 58);
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
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gf_0, hd_0, hd_1, \
                         hf_0, hf_1, hf_2, hf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_2 * hd_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_2 * hd_0[k]
                 + pb_z[k] * hf_2[k];

        t_5[k] = f_1 * hd_1[k]
                 + pb_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_y, pb_z, gf_1, gg_0, gg_3, gg_4, \
                         hd_2, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hd_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = f_1 * hd_2[k]
                 + pb_z[k] * hf_5[k];

        t_8[k] = pa_y[k] * gg_0[k];

        t_9[k] = f_3 * gf_1[k]
                 + pa_y[k] * gg_3[k];

        t_10[k] = pa_y[k] * gg_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pa_y, pb_y, pb_z, fg_5, gf_4, gg_6, \
                         gg_9, hd_4, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * fg_5[k]
                  + pa_x[k] * gg_9[k];

        t_12[k] = f_2 * hd_4[k]
                  + pb_z[k] * hf_8[k];

        t_13[k] = f_2 * gf_4[k]
                  + pb_y[k] * hf_9[k];

        t_14[k] = pa_y[k] * gg_6[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_z, pb_y, pb_z, gf_0, gf_2, gg_0, gg_4, \
                         hd_6, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * gg_0[k];

        t_16[k] = f_2 * gf_0[k]
                  + pb_z[k] * hf_10[k];

        t_17[k] = f_3 * gf_2[k]
                  + pa_z[k] * gg_4[k];

        t_18[k] = f_3 * hd_6[k]
                  + pb_y[k] * hf_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pb_y, fg_0, fg_8, gg_7, gg_13, hd_7, \
                         hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * hd_7[k]
                  + pb_y[k] * hf_12[k];

        t_20[k] = f_1 * fg_8[k]
                  + pa_x[k] * gg_13[k];

        t_21[k] = f_2 * fg_0[k]
                  + pa_y[k] * gg_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, fg_11, gf_11, gf_12, gg_17, \
                         hd_8, hd_9, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * gf_11[k]
                  + f_2 * hd_9[k]
                  + pb_x[k] * hf_16[k];

        t_23[k] = f_2 * hd_8[k]
                  + pb_z[k] * hf_15[k];

        t_24[k] = f_1 * gf_12[k]
                  + pb_x[k] * hf_17[k];

        t_25[k] = f_3 * fg_11[k]
                  + pa_x[k] * gg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, gf_7, gg_8, \
                         gg_11, hd_9, hd_10, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * hd_9[k]
                  + pb_z[k] * hf_18[k];

        t_27[k] = f_3 * gf_7[k]
                  + pb_y[k] * hf_19[k];

        t_28[k] = f_1 * hd_10[k]
                  + pb_z[k] * hf_19[k];

        t_29[k] = pa_y[k] * gg_11[k];

        t_30[k] = pa_z[k] * gg_8[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_y, pa_z, pb_z, fg_12, gf_6, gg_9, \
                         gg_12, gg_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_y[k] * gg_12[k];

        t_32[k] = pa_z[k] * gg_9[k];

        t_33[k] = f_2 * gf_6[k]
                  + pb_z[k] * hf_20[k];

        t_34[k] = f_3 * fg_12[k]
                  + pa_x[k] * gg_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, pb_y, pb_z, fg_0, gf_8, \
                         gf_9, gg_10, gg_13, hf_21, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_2 * gf_9[k]
                  + pb_y[k] * hf_21[k];

        t_36[k] = pa_y[k] * gg_13[k];

        t_37[k] = f_2 * fg_0[k]
                  + pa_z[k] * gg_10[k];

        t_38[k] = pb_y[k] * hf_22[k];

        t_39[k] = f_3 * gf_8[k]
                  + pb_z[k] * hf_22[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pb_y, gf_18, gf_19, hd_11, hd_12, \
                         hd_14, hf_23, hf_24, hf_25, hf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * hd_11[k]
                  + pb_y[k] * hf_23[k];

        t_41[k] = f_1 * gf_18[k]
                  + f_2 * hd_14[k]
                  + pb_x[k] * hf_24[k];

        t_42[k] = f_1 * gf_19[k]
                  + pb_x[k] * hf_28[k];

        t_43[k] = f_1 * hd_12[k]
                  + pb_y[k] * hf_25[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, pb_y, fg_4, fg_15, gg_14, gg_25, \
                         hd_13, hd_14, hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * hd_13[k]
                  + pb_y[k] * hf_26[k];

        t_45[k] = f_2 * hd_14[k]
                  + pb_y[k] * hf_27[k];

        t_46[k] = f_3 * fg_15[k]
                  + pa_x[k] * gg_25[k];

        t_47[k] = f_3 * fg_4[k]
                  + pa_y[k] * gg_14[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, fg_18, gf_21, gf_22, gg_29, \
                         hd_15, hd_16, hf_30, hf_31, hf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * gf_21[k]
                  + f_2 * hd_16[k]
                  + pb_x[k] * hf_31[k];

        t_49[k] = f_2 * hd_15[k]
                  + pb_z[k] * hf_30[k];

        t_50[k] = f_3 * gf_22[k]
                  + pb_x[k] * hf_32[k];

        t_51[k] = f_2 * fg_18[k]
                  + pa_x[k] * gg_29[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_z, pb_y, pb_z, gf_10, gf_13, gg_14, \
                         hd_16, hd_17, hf_33, hf_34, hf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * hd_16[k]
                  + pb_z[k] * hf_33[k];

        t_53[k] = f_1 * gf_13[k]
                  + pb_y[k] * hf_34[k];

        t_54[k] = f_1 * hd_17[k]
                  + pb_z[k] * hf_34[k];

        t_55[k] = pa_z[k] * gg_14[k];

        t_56[k] = f_2 * gf_10[k]
                  + pb_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_y, pa_z, pb_z, fg_7, gf_12, gg_15, gg_17, \
                         gg_18, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * gg_15[k];

        t_58[k] = f_2 * fg_7[k]
                  + pa_y[k] * gg_18[k];

        t_59[k] = pa_z[k] * gg_17[k];

        t_60[k] = f_2 * gf_12[k]
                  + pb_z[k] * hf_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, fg_22, fg_23, gf_15, \
                         gg_20, gg_21, gg_30, gg_31, hf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_2 * fg_22[k]
                  + pa_x[k] * gg_30[k];

        t_62[k] = f_3 * gf_15[k]
                  + pb_y[k] * hf_37[k];

        t_63[k] = f_2 * fg_23[k]
                  + pa_x[k] * gg_31[k];

        t_64[k] = pa_y[k] * gg_20[k];

        t_65[k] = pa_y[k] * gg_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pa_y, pb_z, fg_24, gf_14, gf_17, gg_22, \
                         gg_23, gg_32, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * gf_17[k]
                  + pa_y[k] * gg_22[k];

        t_67[k] = pa_y[k] * gg_23[k];

        t_68[k] = f_2 * fg_24[k]
                  + pa_x[k] * gg_32[k];

        t_69[k] = f_3 * gf_14[k]
                  + pb_z[k] * hf_39[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pa_y, pa_z, pb_y, fg_6, fg_25, gf_19, \
                         gg_20, gg_25, gg_33, hf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * fg_25[k]
                  + pa_x[k] * gg_33[k];

        t_71[k] = f_2 * gf_19[k]
                  + pb_y[k] * hf_40[k];

        t_72[k] = pa_y[k] * gg_25[k];

        t_73[k] = f_3 * fg_6[k]
                  + pa_z[k] * gg_20[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_x, pb_y, pb_z, gf_16, gf_27, hd_18, hd_21, \
                         hf_41, hf_42, hf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_y[k] * hf_41[k];

        t_75[k] = f_1 * gf_16[k]
                  + pb_z[k] * hf_41[k];

        t_76[k] = f_2 * hd_18[k]
                  + pb_y[k] * hf_42[k];

        t_77[k] = f_3 * gf_27[k]
                  + f_2 * hd_21[k]
                  + pb_x[k] * hf_43[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pb_y, gf_28, hd_19, hd_20, hd_21, \
                         hf_44, hf_45, hf_46, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_3 * gf_28[k]
                  + pb_x[k] * hf_47[k];

        t_79[k] = f_1 * hd_19[k]
                  + pb_y[k] * hf_44[k];

        t_80[k] = f_3 * hd_20[k]
                  + pb_y[k] * hf_45[k];

        t_81[k] = f_2 * hd_21[k]
                  + pb_y[k] * hf_46[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pb_z, fg_34, gf_29, gf_30, gg_38, \
                         gg_39, gg_40, hd_22, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_2 * fg_34[k]
                  + pa_x[k] * gg_38[k];

        t_83[k] = f_4 * gf_29[k]
                  + pa_x[k] * gg_39[k];

        t_84[k] = f_3 * gf_30[k]
                  + pa_x[k] * gg_40[k];

        t_85[k] = f_2 * hd_22[k]
                  + pb_z[k] * hf_49[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pa_x, pa_z, pb_x, gf_32, gg_26, \
                         gg_43, gg_45, gg_46, gg_47, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_2 * gf_32[k]
                  + pb_x[k] * hf_51[k];

        t_87[k] = pa_x[k] * gg_43[k];

        t_88[k] = pa_x[k] * gg_45[k];

        t_89[k] = pa_x[k] * gg_46[k];

        t_90[k] = pa_x[k] * gg_47[k];

        t_91[k] = pa_z[k] * gg_26[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pa_x, pa_z, pb_z, gf_20, gf_35, gg_27, \
                         gg_48, gg_50, gg_51, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * gf_20[k]
                  + pb_z[k] * hf_52[k];

        t_93[k] = pa_z[k] * gg_27[k];

        t_94[k] = f_3 * gf_35[k]
                  + pa_x[k] * gg_48[k];

        t_95[k] = pa_x[k] * gg_50[k];

        t_96[k] = pa_x[k] * gg_51[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pa_x, pb_z, gf_23, gf_38, gf_39, \
                         gg_52, gg_53, gg_54, gg_55, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pa_x[k] * gg_52[k];

        t_98[k] = pa_x[k] * gg_53[k];

        t_99[k] = f_4 * gf_38[k]
                  + pa_x[k] * gg_54[k];

        t_100[k] = f_3 * gf_23[k]
                   + pb_z[k] * hf_54[k];

        t_101[k] = f_3 * gf_39[k]
                   + pa_x[k] * gg_55[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pa_x, gf_40, gg_56, gg_57, \
                         gg_58, gg_59, gg_60, gg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_3 * gf_40[k]
                   + pa_x[k] * gg_56[k];

        t_103[k] = pa_x[k] * gg_57[k];

        t_104[k] = pa_x[k] * gg_58[k];

        t_105[k] = pa_x[k] * gg_59[k];

        t_106[k] = pa_x[k] * gg_60[k];

        t_107[k] = pa_x[k] * gg_61[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pa_x, pa_y, gf_44, gg_34, \
                         gg_35, gg_36, gg_62, gg_63, gg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_y[k] * gg_34[k];

        t_109[k] = pa_y[k] * gg_35[k];

        t_110[k] = f_3 * gf_44[k]
                   + pa_x[k] * gg_62[k];

        t_111[k] = pa_y[k] * gg_36[k];

        t_112[k] = pa_x[k] * gg_63[k];

        t_113[k] = pa_x[k] * gg_64[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pa_x, pb_z, gf_26, gf_48, gf_51, \
                         gg_65, gg_66, gg_68, gg_71, hf_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_x[k] * gg_65[k];

        t_115[k] = pa_x[k] * gg_66[k];

        t_116[k] = f_4 * gf_48[k]
                   + pa_x[k] * gg_68[k];

        t_117[k] = f_4 * gf_26[k]
                   + pb_z[k] * hf_58[k];

        t_118[k] = f_3 * gf_51[k]
                   + pa_x[k] * gg_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, pa_x, pb_x, gf_55, gg_74, \
                         gg_75, gg_76, gg_78, hd_26, hf_60, hf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * gf_55[k]
                   + pb_x[k] * hf_60[k];

        t_120[k] = pa_x[k] * gg_74[k];

        t_121[k] = pa_x[k] * gg_75[k];

        t_122[k] = pa_x[k] * gg_76[k];

        t_123[k] = pa_x[k] * gg_78[k];

        t_124[k] = f_1 * hd_26[k]
                   + pb_x[k] * hf_61[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, pb_x, pb_z, hd_27, hd_28, \
                         hd_29, hf_62, hf_63, hf_64, hf_65, hf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * hd_27[k]
                   + pb_x[k] * hf_62[k];

        t_126[k] = f_2 * hd_28[k]
                   + pb_x[k] * hf_63[k];

        t_127[k] = pb_z[k] * hf_62[k];

        t_128[k] = f_2 * hd_29[k]
                   + pb_x[k] * hf_64[k];

        t_129[k] = pb_x[k] * hf_65[k];

        t_130[k] = pb_x[k] * hf_67[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, pb_x, pb_y, pb_z, gf_32, \
                         gf_34, hd_28, hd_29, hf_65, hf_66, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_x[k] * hf_68[k];

        t_132[k] = f_0 * gf_32[k]
                   + f_1 * hd_28[k]
                   + pb_y[k] * hf_65[k];

        t_133[k] = pb_z[k] * hf_65[k];

        t_134[k] = f_2 * hd_28[k]
                   + pb_z[k] * hf_66[k];

        t_135[k] = f_0 * gf_34[k]
                   + pb_y[k] * hf_68[k];

        t_136[k] = f_1 * hd_29[k]
                   + pb_z[k] * hf_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pa_z, pb_x, pb_z, gf_32, gf_33, \
                         gg_43, gg_45, hd_31, hf_69, hf_70, hf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * hd_31[k]
                   + pb_x[k] * hf_69[k];

        t_138[k] = pb_x[k] * hf_71[k];

        t_139[k] = pa_z[k] * gg_43[k];

        t_140[k] = f_2 * gf_32[k]
                   + pb_z[k] * hf_70[k];

        t_141[k] = f_3 * gf_33[k]
                   + pa_z[k] * gg_45[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_y, pb_x, pb_y, fg_23, gf_37, gg_53, \
                         hd_32, hd_33, hf_71, hf_72, hf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_4 * gf_37[k]
                   + pb_y[k] * hf_71[k];

        t_143[k] = f_1 * fg_23[k]
                   + pa_y[k] * gg_53[k];

        t_144[k] = f_1 * hd_32[k]
                   + pb_x[k] * hf_72[k];

        t_145[k] = f_2 * hd_33[k]
                   + pb_x[k] * hf_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, pa_z, pb_x, pb_z, fg_18, gf_36, \
                         gg_49, hd_34, hf_74, hf_75, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_2 * hd_34[k]
                   + pb_x[k] * hf_74[k];

        t_147[k] = pb_x[k] * hf_75[k];

        t_148[k] = pb_x[k] * hf_77[k];

        t_149[k] = f_2 * fg_18[k]
                   + pa_z[k] * gg_49[k];

        t_150[k] = f_3 * gf_36[k]
                   + pb_z[k] * hf_75[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pb_x, pb_y, fg_27, gf_42, gf_43, \
                         gg_61, hd_34, hd_35, hf_76, hf_77, hf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_1 * gf_42[k]
                   + f_2 * hd_34[k]
                   + pb_y[k] * hf_76[k];

        t_152[k] = f_1 * gf_43[k]
                   + pb_y[k] * hf_77[k];

        t_153[k] = f_3 * fg_27[k]
                   + pa_y[k] * gg_61[k];

        t_154[k] = f_1 * hd_35[k]
                   + pb_x[k] * hf_78[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_z, pb_x, fg_21, gg_57, hd_36, \
                         hd_37, hf_79, hf_80, hf_81, hf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_2 * hd_36[k]
                   + pb_x[k] * hf_79[k];

        t_156[k] = f_2 * hd_37[k]
                   + pb_x[k] * hf_80[k];

        t_157[k] = pb_x[k] * hf_81[k];

        t_158[k] = pb_x[k] * hf_83[k];

        t_159[k] = f_3 * fg_21[k]
                   + pa_z[k] * gg_57[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_y, pb_y, pb_z, fg_34, gf_41, gf_46, \
                         gf_47, gg_67, hd_37, hf_81, hf_82, hf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_1 * gf_41[k]
                   + pb_z[k] * hf_81[k];

        t_161[k] = f_3 * gf_46[k]
                   + f_2 * hd_37[k]
                   + pb_y[k] * hf_82[k];

        t_162[k] = f_3 * gf_47[k]
                   + pb_y[k] * hf_83[k];

        t_163[k] = f_2 * fg_34[k]
                   + pa_y[k] * gg_67[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_y, pb_x, pb_z, gf_45, gf_52, \
                         gf_54, gg_74, gg_76, hd_38, hf_84, hf_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_2 * hd_38[k]
                   + pb_x[k] * hf_84[k];

        t_165[k] = pb_x[k] * hf_85[k];

        t_166[k] = f_4 * gf_52[k]
                   + pa_y[k] * gg_74[k];

        t_167[k] = f_4 * gf_45[k]
                   + pb_z[k] * hf_85[k];

        t_168[k] = f_3 * gf_54[k]
                   + pa_y[k] * gg_76[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_y, pb_x, pb_y, gf_55, gg_78, \
                         hd_40, hd_41, hf_87, hf_88, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_2 * gf_55[k]
                   + pb_y[k] * hf_87[k];

        t_170[k] = pa_y[k] * gg_78[k];

        t_171[k] = f_1 * hd_40[k]
                   + pb_x[k] * hf_88[k];

        t_172[k] = pb_y[k] * hf_88[k];

        t_173[k] = f_3 * hd_41[k]
                   + pb_x[k] * hf_89[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, pb_x, pb_y, hd_42, hd_44, \
                         hf_89, hf_90, hf_91, hf_92, hf_93, hf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_2 * hd_42[k]
                   + pb_x[k] * hf_90[k];

        t_175[k] = pb_y[k] * hf_89[k];

        t_176[k] = f_2 * hd_44[k]
                   + pb_x[k] * hf_91[k];

        t_177[k] = pb_x[k] * hf_92[k];

        t_178[k] = pb_x[k] * hf_93[k];

        t_179[k] = pb_x[k] * hf_95[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pb_y, pb_z, gf_55, hd_42, hd_43, \
                         hd_44, hf_92, hf_93, hf_94, hf_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * hd_42[k]
                   + pb_y[k] * hf_92[k];

        t_181[k] = f_3 * hd_43[k]
                   + pb_y[k] * hf_93[k];

        t_182[k] = f_2 * hd_44[k]
                   + pb_y[k] * hf_94[k];

        t_183[k] = pb_y[k] * hf_95[k];

        t_184[k] = f_0 * gf_55[k]
                   + f_1 * hd_44[k]
                   + pb_z[k] * hf_95[k];
    }
}

auto
compute_prim_hg_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fg, const size_t gf, const size_t gg,
                          const size_t hd, const size_t hf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_40 = buffer.data(fg + 40);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_74 = buffer.data(gg + 74);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
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
    const auto *hf_42 = buffer.data(hf + 42);
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
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gf_0, hd_0, hd_1, \
                         hf_0, hf_1, hf_2, hf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 + f_1 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = pb_y[k] * hf_0[k];

        t_2[k] = pb_z[k] * hf_0[k];

        t_3[k] = f_2 * hd_0[k]
                 + pb_y[k] * hf_1[k];

        t_4[k] = f_2 * hd_0[k]
                 + pb_z[k] * hf_2[k];

        t_5[k] = f_1 * hd_1[k]
                 + pb_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_y, pb_z, fg_7, gg_0, gg_8, \
                         hd_2, hf_4, hf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hd_2[k]
                 + pb_y[k] * hf_4[k];

        t_7[k] = pb_y[k] * hf_5[k];

        t_8[k] = f_1 * hd_2[k]
                 + pb_z[k] * hf_5[k];

        t_9[k] = pa_y[k] * gg_0[k];

        t_10[k] = f_1 * fg_7[k]
                  + pa_x[k] * gg_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, pb_z, gf_2, gg_0, gg_4, \
                         gg_6, hd_4, hf_7, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * hf_7[k];

        t_12[k] = f_2 * hd_4[k]
                  + pb_z[k] * hf_8[k];

        t_13[k] = pa_y[k] * gg_6[k];

        t_14[k] = pa_z[k] * gg_0[k];

        t_15[k] = f_3 * gf_2[k]
                  + pa_z[k] * gg_4[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, fg_9, gg_11, hd_6, hd_7, hf_10, \
                         hf_11, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * hd_6[k]
                  + pb_y[k] * hf_10[k];

        t_17[k] = f_2 * hd_7[k]
                  + pb_y[k] * hf_11[k];

        t_18[k] = pb_y[k] * hf_12[k];

        t_19[k] = f_1 * fg_9[k]
                  + pa_x[k] * gg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pb_x, pb_z, fg_0, gf_10, gg_7, hd_8, \
                         hd_9, hf_13, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * fg_0[k]
                  + pa_y[k] * gg_7[k];

        t_21[k] = pb_z[k] * hf_13[k];

        t_22[k] = f_1 * gf_10[k]
                  + f_2 * hd_9[k]
                  + pb_x[k] * hf_15[k];

        t_23[k] = f_2 * hd_8[k]
                  + pb_z[k] * hf_14[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, pb_z, fg_13, gf_11, gg_15, \
                         hd_9, hd_10, hf_16, hf_17, hf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * gf_11[k]
                  + pb_x[k] * hf_16[k];

        t_25[k] = f_3 * fg_13[k]
                  + pa_x[k] * gg_15[k];

        t_26[k] = pb_z[k] * hf_16[k];

        t_27[k] = f_2 * hd_9[k]
                  + pb_z[k] * hf_17[k];

        t_28[k] = f_1 * hd_10[k]
                  + pb_z[k] * hf_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_y, pa_z, pb_y, fg_0, gg_8, gg_10, \
                         gg_11, hd_11, hf_19, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * gg_8[k];

        t_30[k] = pa_y[k] * gg_11[k];

        t_31[k] = f_2 * fg_0[k]
                  + pa_z[k] * gg_10[k];

        t_32[k] = pb_y[k] * hf_19[k];

        t_33[k] = f_2 * hd_11[k]
                  + pb_y[k] * hf_20[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, pb_y, gf_15, gf_16, hd_12, hd_13, \
                         hd_14, hf_21, hf_22, hf_23, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * gf_15[k]
                  + f_2 * hd_14[k]
                  + pb_x[k] * hf_21[k];

        t_35[k] = f_1 * gf_16[k]
                  + pb_x[k] * hf_25[k];

        t_36[k] = f_1 * hd_12[k]
                  + pb_y[k] * hf_22[k];

        t_37[k] = f_3 * hd_13[k]
                  + pb_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pa_y, pb_y, fg_6, fg_17, gg_12, gg_24, \
                         hd_14, hf_24, hf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * hd_14[k]
                  + pb_y[k] * hf_24[k];

        t_39[k] = pb_y[k] * hf_25[k];

        t_40[k] = f_3 * fg_17[k]
                  + pa_x[k] * gg_24[k];

        t_41[k] = f_3 * fg_6[k]
                  + pa_y[k] * gg_12[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_x, pb_z, gf_17, gf_18, hd_15, hd_16, \
                         hf_26, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * hf_26[k];

        t_43[k] = f_3 * gf_17[k]
                  + f_2 * hd_16[k]
                  + pb_x[k] * hf_28[k];

        t_44[k] = f_2 * hd_15[k]
                  + pb_z[k] * hf_27[k];

        t_45[k] = f_3 * gf_18[k]
                  + pb_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pa_z, pb_z, fg_21, gg_12, gg_28, \
                         hd_16, hd_17, hf_29, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * fg_21[k]
                  + pa_x[k] * gg_28[k];

        t_47[k] = pb_z[k] * hf_29[k];

        t_48[k] = f_2 * hd_16[k]
                  + pb_z[k] * hf_30[k];

        t_49[k] = f_1 * hd_17[k]
                  + pb_z[k] * hf_31[k];

        t_50[k] = pa_z[k] * gg_12[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_x, pa_y, pa_z, fg_26, fg_27, gg_15, \
                         gg_19, gg_24, gg_30, gg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * gg_15[k];

        t_52[k] = f_2 * fg_26[k]
                  + pa_x[k] * gg_30[k];

        t_53[k] = pa_y[k] * gg_19[k];

        t_54[k] = f_2 * fg_27[k]
                  + pa_x[k] * gg_31[k];

        t_55[k] = pa_y[k] * gg_24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, pb_y, fg_8, gf_19, gg_19, hd_18, \
                         hd_21, hf_32, hf_33, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * fg_8[k]
                  + pa_z[k] * gg_19[k];

        t_57[k] = pb_y[k] * hf_32[k];

        t_58[k] = f_2 * hd_18[k]
                  + pb_y[k] * hf_33[k];

        t_59[k] = f_3 * gf_19[k]
                  + f_2 * hd_21[k]
                  + pb_x[k] * hf_34[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pb_x, pb_y, gf_20, hd_19, hd_20, hd_21, \
                         hf_35, hf_36, hf_37, hf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * gf_20[k]
                  + pb_x[k] * hf_38[k];

        t_61[k] = f_1 * hd_19[k]
                  + pb_y[k] * hf_35[k];

        t_62[k] = f_3 * hd_20[k]
                  + pb_y[k] * hf_36[k];

        t_63[k] = f_2 * hd_21[k]
                  + pb_y[k] * hf_37[k];

        t_64[k] = pb_y[k] * hf_38[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pb_z, fg_40, gf_21, gf_22, gg_35, \
                         gg_36, gg_37, hd_22, hf_39, hf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * fg_40[k]
                  + pa_x[k] * gg_35[k];

        t_66[k] = f_4 * gf_21[k]
                  + pa_x[k] * gg_36[k];

        t_67[k] = pb_z[k] * hf_39[k];

        t_68[k] = f_3 * gf_22[k]
                  + pa_x[k] * gg_37[k];

        t_69[k] = f_2 * hd_22[k]
                  + pb_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_x, pa_z, pb_x, gf_24, gf_30, gg_25, \
                         gg_41, gg_48, gg_49, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * gf_24[k]
                  + pb_x[k] * hf_42[k];

        t_71[k] = pa_x[k] * gg_41[k];

        t_72[k] = pa_z[k] * gg_25[k];

        t_73[k] = pa_x[k] * gg_48[k];

        t_74[k] = f_4 * gf_30[k]
                  + pa_x[k] * gg_49[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, pa_x, gf_40, gf_43, gg_54, gg_55, \
                         gg_57, gg_60, gg_64, gg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pa_x[k] * gg_54[k];

        t_76[k] = pa_x[k] * gg_55[k];

        t_77[k] = pa_x[k] * gg_57[k];

        t_78[k] = pa_x[k] * gg_60[k];

        t_79[k] = f_4 * gf_40[k]
                  + pa_x[k] * gg_64[k];

        t_80[k] = f_3 * gf_43[k]
                  + pa_x[k] * gg_67[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_x, pb_x, gf_47, gg_74, hd_26, hd_27, \
                         hd_28, hf_45, hf_46, hf_47, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_2 * gf_47[k]
                  + pb_x[k] * hf_45[k];

        t_82[k] = pa_x[k] * gg_74[k];

        t_83[k] = f_1 * hd_26[k]
                  + pb_x[k] * hf_46[k];

        t_84[k] = f_3 * hd_27[k]
                  + pb_x[k] * hf_47[k];

        t_85[k] = f_2 * hd_28[k]
                  + pb_x[k] * hf_48[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pb_x, pb_y, pb_z, gf_24, hd_28, \
                         hd_29, hf_49, hf_50, hf_52, hf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_2 * hd_29[k]
                  + pb_x[k] * hf_49[k];

        t_87[k] = pb_x[k] * hf_50[k];

        t_88[k] = pb_x[k] * hf_52[k];

        t_89[k] = pb_x[k] * hf_53[k];

        t_90[k] = f_0 * gf_24[k]
                  + f_1 * hd_28[k]
                  + pb_y[k] * hf_50[k];

        t_91[k] = pb_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pb_x, pb_y, pb_z, gf_26, hd_28, hd_29, \
                         hd_31, hf_51, hf_53, hf_54, hf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * hd_28[k]
                  + pb_z[k] * hf_51[k];

        t_93[k] = f_0 * gf_26[k]
                  + pb_y[k] * hf_53[k];

        t_94[k] = f_1 * hd_29[k]
                  + pb_z[k] * hf_53[k];

        t_95[k] = f_2 * hd_31[k]
                  + pb_x[k] * hf_54[k];

        t_96[k] = pb_x[k] * hf_56[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pa_z, pb_x, fg_26, gf_25, gg_41, \
                         gg_43, gg_48, hd_32, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pa_z[k] * gg_41[k];

        t_98[k] = f_3 * gf_25[k]
                  + pa_z[k] * gg_43[k];

        t_99[k] = f_1 * fg_26[k]
                  + pa_y[k] * gg_48[k];

        t_100[k] = f_1 * hd_32[k]
                   + pb_x[k] * hf_57[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pa_z, pb_x, fg_21, gg_47, hd_33, \
                         hd_34, hf_58, hf_59, hf_60, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * hd_33[k]
                   + pb_x[k] * hf_58[k];

        t_102[k] = f_2 * hd_34[k]
                   + pb_x[k] * hf_59[k];

        t_103[k] = pb_x[k] * hf_60[k];

        t_104[k] = pb_x[k] * hf_62[k];

        t_105[k] = f_2 * fg_21[k]
                   + pa_z[k] * gg_47[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_y, pb_x, pb_y, fg_30, gf_34, gf_35, \
                         gg_57, hd_34, hd_35, hf_61, hf_62, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_1 * gf_34[k]
                   + f_2 * hd_34[k]
                   + pb_y[k] * hf_61[k];

        t_107[k] = f_1 * gf_35[k]
                   + pb_y[k] * hf_62[k];

        t_108[k] = f_3 * fg_30[k]
                   + pa_y[k] * gg_57[k];

        t_109[k] = f_1 * hd_35[k]
                   + pb_x[k] * hf_63[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_z, pb_x, fg_25, gg_54, hd_36, \
                         hd_37, hf_64, hf_65, hf_66, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_2 * hd_36[k]
                   + pb_x[k] * hf_64[k];

        t_111[k] = f_2 * hd_37[k]
                   + pb_x[k] * hf_65[k];

        t_112[k] = pb_x[k] * hf_66[k];

        t_113[k] = pb_x[k] * hf_68[k];

        t_114[k] = f_3 * fg_25[k]
                   + pa_z[k] * gg_54[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_y, pb_x, pb_y, fg_40, gf_38, gf_39, \
                         gg_63, hd_37, hd_38, hf_67, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * gf_38[k]
                   + f_2 * hd_37[k]
                   + pb_y[k] * hf_67[k];

        t_116[k] = f_3 * gf_39[k]
                   + pb_y[k] * hf_68[k];

        t_117[k] = f_2 * fg_40[k]
                   + pa_y[k] * gg_63[k];

        t_118[k] = f_2 * hd_38[k]
                   + pb_x[k] * hf_69[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pa_y, pb_x, pb_y, gf_44, gf_46, \
                         gf_47, gg_70, gg_72, gg_74, hf_70, hf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_x[k] * hf_70[k];

        t_120[k] = f_4 * gf_44[k]
                   + pa_y[k] * gg_70[k];

        t_121[k] = f_3 * gf_46[k]
                   + pa_y[k] * gg_72[k];

        t_122[k] = f_2 * gf_47[k]
                   + pb_y[k] * hf_72[k];

        t_123[k] = pa_y[k] * gg_74[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pb_x, hd_40, hd_41, hd_42, hd_44, \
                         hf_73, hf_74, hf_75, hf_76, hf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * hd_40[k]
                   + pb_x[k] * hf_73[k];

        t_125[k] = f_3 * hd_41[k]
                   + pb_x[k] * hf_74[k];

        t_126[k] = f_2 * hd_42[k]
                   + pb_x[k] * hf_75[k];

        t_127[k] = f_2 * hd_44[k]
                   + pb_x[k] * hf_76[k];

        t_128[k] = pb_x[k] * hf_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, t_134, pb_x, pb_y, hd_42, hd_43, \
                         hd_44, hf_77, hf_78, hf_79, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pb_x[k] * hf_78[k];

        t_130[k] = pb_x[k] * hf_80[k];

        t_131[k] = f_1 * hd_42[k]
                   + pb_y[k] * hf_77[k];

        t_132[k] = f_3 * hd_43[k]
                   + pb_y[k] * hf_78[k];

        t_133[k] = f_2 * hd_44[k]
                   + pb_y[k] * hf_79[k];

        t_134[k] = pb_y[k] * hf_80[k];
    }

#pragma omp simd aligned(t_135, pb_z, gf_47, hd_44, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * gf_47[k]
                   + f_1 * hd_44[k]
                   + pb_z[k] * hf_80[k];
    }
}

}  // namespace simdovl
