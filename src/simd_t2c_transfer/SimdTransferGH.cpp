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


#include "SimdTransferGH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gh(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fh, const size_t fi, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

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
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_172 = buffer.data(fi + 172);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_176 = buffer.data(fi + 176);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_181 = buffer.data(fi + 181);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_184 = buffer.data(fi + 184);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_187 = buffer.data(fi + 187);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_197 = buffer.data(fi + 197);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_200 = buffer.data(fi + 200);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_203 = buffer.data(fi + 203);
    const auto *fi_204 = buffer.data(fi + 204);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_207 = buffer.data(fi + 207);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_209 = buffer.data(fi + 209);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_211 = buffer.data(fi + 211);
    const auto *fi_212 = buffer.data(fi + 212);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_215 = buffer.data(fi + 215);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_225 = buffer.data(fi + 225);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_228 = buffer.data(fi + 228);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_231 = buffer.data(fi + 231);
    const auto *fi_232 = buffer.data(fi + 232);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_235 = buffer.data(fi + 235);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_237 = buffer.data(fi + 237);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_240 = buffer.data(fi + 240);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_243 = buffer.data(fi + 243);
    const auto *fi_244 = buffer.data(fi + 244);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_256 = buffer.data(fi + 256);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_259 = buffer.data(fi + 259);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_263 = buffer.data(fi + 263);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_268 = buffer.data(fi + 268);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_271 = buffer.data(fi + 271);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fh_0, fh_1, fh_2, fh_3, fh_4, fi_0, \
                         fi_1, fi_2, fi_3, fi_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * fh_0[k]
                 + fi_0[k];

        t_1[k] = -ab_x[k] * fh_1[k]
                 + fi_1[k];

        t_2[k] = -ab_x[k] * fh_2[k]
                 + fi_2[k];

        t_3[k] = -ab_x[k] * fh_3[k]
                 + fi_3[k];

        t_4[k] = -ab_x[k] * fh_4[k]
                 + fi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fh_5, fh_6, fh_7, fh_8, fh_9, fi_5, \
                         fi_6, fi_7, fi_8, fi_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * fh_5[k]
                 + fi_5[k];

        t_6[k] = -ab_x[k] * fh_6[k]
                 + fi_6[k];

        t_7[k] = -ab_x[k] * fh_7[k]
                 + fi_7[k];

        t_8[k] = -ab_x[k] * fh_8[k]
                 + fi_8[k];

        t_9[k] = -ab_x[k] * fh_9[k]
                 + fi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fh_10, fh_11, fh_12, fh_13, \
                         fh_14, fi_10, fi_11, fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * fh_10[k]
                  + fi_10[k];

        t_11[k] = -ab_x[k] * fh_11[k]
                  + fi_11[k];

        t_12[k] = -ab_x[k] * fh_12[k]
                  + fi_12[k];

        t_13[k] = -ab_x[k] * fh_13[k]
                  + fi_13[k];

        t_14[k] = -ab_x[k] * fh_14[k]
                  + fi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fh_15, fh_16, fh_17, fh_18, \
                         fh_19, fi_15, fi_16, fi_17, fi_18, fi_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * fh_15[k]
                  + fi_15[k];

        t_16[k] = -ab_x[k] * fh_16[k]
                  + fi_16[k];

        t_17[k] = -ab_x[k] * fh_17[k]
                  + fi_17[k];

        t_18[k] = -ab_x[k] * fh_18[k]
                  + fi_18[k];

        t_19[k] = -ab_x[k] * fh_19[k]
                  + fi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fh_20, fh_21, fh_22, fh_23, \
                         fh_24, fi_20, fi_28, fi_29, fi_30, fi_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * fh_20[k]
                  + fi_20[k];

        t_21[k] = -ab_x[k] * fh_21[k]
                  + fi_28[k];

        t_22[k] = -ab_x[k] * fh_22[k]
                  + fi_29[k];

        t_23[k] = -ab_x[k] * fh_23[k]
                  + fi_30[k];

        t_24[k] = -ab_x[k] * fh_24[k]
                  + fi_31[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fh_25, fh_26, fh_27, fh_28, \
                         fh_29, fi_32, fi_33, fi_34, fi_35, fi_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * fh_25[k]
                  + fi_32[k];

        t_26[k] = -ab_x[k] * fh_26[k]
                  + fi_33[k];

        t_27[k] = -ab_x[k] * fh_27[k]
                  + fi_34[k];

        t_28[k] = -ab_x[k] * fh_28[k]
                  + fi_35[k];

        t_29[k] = -ab_x[k] * fh_29[k]
                  + fi_36[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fh_30, fh_31, fh_32, fh_33, \
                         fh_34, fi_37, fi_38, fi_39, fi_40, fi_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * fh_30[k]
                  + fi_37[k];

        t_31[k] = -ab_x[k] * fh_31[k]
                  + fi_38[k];

        t_32[k] = -ab_x[k] * fh_32[k]
                  + fi_39[k];

        t_33[k] = -ab_x[k] * fh_33[k]
                  + fi_40[k];

        t_34[k] = -ab_x[k] * fh_34[k]
                  + fi_41[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, fh_35, fh_36, fh_37, fh_38, \
                         fh_39, fi_42, fi_43, fi_44, fi_45, fi_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * fh_35[k]
                  + fi_42[k];

        t_36[k] = -ab_x[k] * fh_36[k]
                  + fi_43[k];

        t_37[k] = -ab_x[k] * fh_37[k]
                  + fi_44[k];

        t_38[k] = -ab_x[k] * fh_38[k]
                  + fi_45[k];

        t_39[k] = -ab_x[k] * fh_39[k]
                  + fi_46[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fh_40, fh_41, fh_42, fh_43, \
                         fh_44, fi_47, fi_48, fi_56, fi_57, fi_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * fh_40[k]
                  + fi_47[k];

        t_41[k] = -ab_x[k] * fh_41[k]
                  + fi_48[k];

        t_42[k] = -ab_x[k] * fh_42[k]
                  + fi_56[k];

        t_43[k] = -ab_x[k] * fh_43[k]
                  + fi_57[k];

        t_44[k] = -ab_x[k] * fh_44[k]
                  + fi_58[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, fh_45, fh_46, fh_47, fh_48, \
                         fh_49, fi_59, fi_60, fi_61, fi_62, fi_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * fh_45[k]
                  + fi_59[k];

        t_46[k] = -ab_x[k] * fh_46[k]
                  + fi_60[k];

        t_47[k] = -ab_x[k] * fh_47[k]
                  + fi_61[k];

        t_48[k] = -ab_x[k] * fh_48[k]
                  + fi_62[k];

        t_49[k] = -ab_x[k] * fh_49[k]
                  + fi_63[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fh_50, fh_51, fh_52, fh_53, \
                         fh_54, fi_64, fi_65, fi_66, fi_67, fi_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * fh_50[k]
                  + fi_64[k];

        t_51[k] = -ab_x[k] * fh_51[k]
                  + fi_65[k];

        t_52[k] = -ab_x[k] * fh_52[k]
                  + fi_66[k];

        t_53[k] = -ab_x[k] * fh_53[k]
                  + fi_67[k];

        t_54[k] = -ab_x[k] * fh_54[k]
                  + fi_68[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, fh_55, fh_56, fh_57, fh_58, \
                         fh_59, fi_69, fi_70, fi_71, fi_72, fi_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * fh_55[k]
                  + fi_69[k];

        t_56[k] = -ab_x[k] * fh_56[k]
                  + fi_70[k];

        t_57[k] = -ab_x[k] * fh_57[k]
                  + fi_71[k];

        t_58[k] = -ab_x[k] * fh_58[k]
                  + fi_72[k];

        t_59[k] = -ab_x[k] * fh_59[k]
                  + fi_73[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fh_60, fh_61, fh_62, fh_63, \
                         fh_64, fi_74, fi_75, fi_76, fi_84, fi_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * fh_60[k]
                  + fi_74[k];

        t_61[k] = -ab_x[k] * fh_61[k]
                  + fi_75[k];

        t_62[k] = -ab_x[k] * fh_62[k]
                  + fi_76[k];

        t_63[k] = -ab_x[k] * fh_63[k]
                  + fi_84[k];

        t_64[k] = -ab_x[k] * fh_64[k]
                  + fi_85[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, fh_65, fh_66, fh_67, fh_68, \
                         fh_69, fi_86, fi_87, fi_88, fi_89, fi_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * fh_65[k]
                  + fi_86[k];

        t_66[k] = -ab_x[k] * fh_66[k]
                  + fi_87[k];

        t_67[k] = -ab_x[k] * fh_67[k]
                  + fi_88[k];

        t_68[k] = -ab_x[k] * fh_68[k]
                  + fi_89[k];

        t_69[k] = -ab_x[k] * fh_69[k]
                  + fi_90[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fh_70, fh_71, fh_72, fh_73, \
                         fh_74, fi_91, fi_92, fi_93, fi_94, fi_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * fh_70[k]
                  + fi_91[k];

        t_71[k] = -ab_x[k] * fh_71[k]
                  + fi_92[k];

        t_72[k] = -ab_x[k] * fh_72[k]
                  + fi_93[k];

        t_73[k] = -ab_x[k] * fh_73[k]
                  + fi_94[k];

        t_74[k] = -ab_x[k] * fh_74[k]
                  + fi_95[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, fh_75, fh_76, fh_77, fh_78, \
                         fh_79, fi_96, fi_97, fi_98, fi_99, fi_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * fh_75[k]
                  + fi_96[k];

        t_76[k] = -ab_x[k] * fh_76[k]
                  + fi_97[k];

        t_77[k] = -ab_x[k] * fh_77[k]
                  + fi_98[k];

        t_78[k] = -ab_x[k] * fh_78[k]
                  + fi_99[k];

        t_79[k] = -ab_x[k] * fh_79[k]
                  + fi_100[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fh_80, fh_81, fh_82, fh_83, \
                         fh_84, fi_101, fi_102, fi_103, fi_104, \
                         fi_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * fh_80[k]
                  + fi_101[k];

        t_81[k] = -ab_x[k] * fh_81[k]
                  + fi_102[k];

        t_82[k] = -ab_x[k] * fh_82[k]
                  + fi_103[k];

        t_83[k] = -ab_x[k] * fh_83[k]
                  + fi_104[k];

        t_84[k] = -ab_x[k] * fh_84[k]
                  + fi_112[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, fh_85, fh_86, fh_87, fh_88, \
                         fh_89, fi_113, fi_114, fi_115, fi_116, \
                         fi_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * fh_85[k]
                  + fi_113[k];

        t_86[k] = -ab_x[k] * fh_86[k]
                  + fi_114[k];

        t_87[k] = -ab_x[k] * fh_87[k]
                  + fi_115[k];

        t_88[k] = -ab_x[k] * fh_88[k]
                  + fi_116[k];

        t_89[k] = -ab_x[k] * fh_89[k]
                  + fi_117[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fh_90, fh_91, fh_92, fh_93, \
                         fh_94, fi_118, fi_119, fi_120, fi_121, \
                         fi_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * fh_90[k]
                  + fi_118[k];

        t_91[k] = -ab_x[k] * fh_91[k]
                  + fi_119[k];

        t_92[k] = -ab_x[k] * fh_92[k]
                  + fi_120[k];

        t_93[k] = -ab_x[k] * fh_93[k]
                  + fi_121[k];

        t_94[k] = -ab_x[k] * fh_94[k]
                  + fi_122[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, fh_95, fh_96, fh_97, fh_98, \
                         fh_99, fi_123, fi_124, fi_125, fi_126, \
                         fi_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * fh_95[k]
                  + fi_123[k];

        t_96[k] = -ab_x[k] * fh_96[k]
                  + fi_124[k];

        t_97[k] = -ab_x[k] * fh_97[k]
                  + fi_125[k];

        t_98[k] = -ab_x[k] * fh_98[k]
                  + fi_126[k];

        t_99[k] = -ab_x[k] * fh_99[k]
                  + fi_127[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, fh_100, fh_101, fh_102, \
                         fh_103, fh_104, fi_128, fi_129, fi_130, fi_131, \
                         fi_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * fh_100[k]
                   + fi_128[k];

        t_101[k] = -ab_x[k] * fh_101[k]
                   + fi_129[k];

        t_102[k] = -ab_x[k] * fh_102[k]
                   + fi_130[k];

        t_103[k] = -ab_x[k] * fh_103[k]
                   + fi_131[k];

        t_104[k] = -ab_x[k] * fh_104[k]
                   + fi_132[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, fh_105, fh_106, fh_107, \
                         fh_108, fh_109, fi_140, fi_141, fi_142, fi_143, \
                         fi_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * fh_105[k]
                   + fi_140[k];

        t_106[k] = -ab_x[k] * fh_106[k]
                   + fi_141[k];

        t_107[k] = -ab_x[k] * fh_107[k]
                   + fi_142[k];

        t_108[k] = -ab_x[k] * fh_108[k]
                   + fi_143[k];

        t_109[k] = -ab_x[k] * fh_109[k]
                   + fi_144[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, fh_110, fh_111, fh_112, \
                         fh_113, fh_114, fi_145, fi_146, fi_147, fi_148, \
                         fi_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * fh_110[k]
                   + fi_145[k];

        t_111[k] = -ab_x[k] * fh_111[k]
                   + fi_146[k];

        t_112[k] = -ab_x[k] * fh_112[k]
                   + fi_147[k];

        t_113[k] = -ab_x[k] * fh_113[k]
                   + fi_148[k];

        t_114[k] = -ab_x[k] * fh_114[k]
                   + fi_149[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, fh_115, fh_116, fh_117, \
                         fh_118, fh_119, fi_150, fi_151, fi_152, fi_153, \
                         fi_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * fh_115[k]
                   + fi_150[k];

        t_116[k] = -ab_x[k] * fh_116[k]
                   + fi_151[k];

        t_117[k] = -ab_x[k] * fh_117[k]
                   + fi_152[k];

        t_118[k] = -ab_x[k] * fh_118[k]
                   + fi_153[k];

        t_119[k] = -ab_x[k] * fh_119[k]
                   + fi_154[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, fh_120, fh_121, fh_122, \
                         fh_123, fh_124, fi_155, fi_156, fi_157, fi_158, \
                         fi_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * fh_120[k]
                   + fi_155[k];

        t_121[k] = -ab_x[k] * fh_121[k]
                   + fi_156[k];

        t_122[k] = -ab_x[k] * fh_122[k]
                   + fi_157[k];

        t_123[k] = -ab_x[k] * fh_123[k]
                   + fi_158[k];

        t_124[k] = -ab_x[k] * fh_124[k]
                   + fi_159[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, fh_125, fh_126, fh_127, \
                         fh_128, fh_129, fi_160, fi_168, fi_169, fi_170, \
                         fi_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * fh_125[k]
                   + fi_160[k];

        t_126[k] = -ab_x[k] * fh_126[k]
                   + fi_168[k];

        t_127[k] = -ab_x[k] * fh_127[k]
                   + fi_169[k];

        t_128[k] = -ab_x[k] * fh_128[k]
                   + fi_170[k];

        t_129[k] = -ab_x[k] * fh_129[k]
                   + fi_171[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, fh_130, fh_131, fh_132, \
                         fh_133, fh_134, fi_172, fi_173, fi_174, fi_175, \
                         fi_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_x[k] * fh_130[k]
                   + fi_172[k];

        t_131[k] = -ab_x[k] * fh_131[k]
                   + fi_173[k];

        t_132[k] = -ab_x[k] * fh_132[k]
                   + fi_174[k];

        t_133[k] = -ab_x[k] * fh_133[k]
                   + fi_175[k];

        t_134[k] = -ab_x[k] * fh_134[k]
                   + fi_176[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, fh_135, fh_136, fh_137, \
                         fh_138, fh_139, fi_177, fi_178, fi_179, fi_180, \
                         fi_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_x[k] * fh_135[k]
                   + fi_177[k];

        t_136[k] = -ab_x[k] * fh_136[k]
                   + fi_178[k];

        t_137[k] = -ab_x[k] * fh_137[k]
                   + fi_179[k];

        t_138[k] = -ab_x[k] * fh_138[k]
                   + fi_180[k];

        t_139[k] = -ab_x[k] * fh_139[k]
                   + fi_181[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, fh_140, fh_141, fh_142, \
                         fh_143, fh_144, fi_182, fi_183, fi_184, fi_185, \
                         fi_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_x[k] * fh_140[k]
                   + fi_182[k];

        t_141[k] = -ab_x[k] * fh_141[k]
                   + fi_183[k];

        t_142[k] = -ab_x[k] * fh_142[k]
                   + fi_184[k];

        t_143[k] = -ab_x[k] * fh_143[k]
                   + fi_185[k];

        t_144[k] = -ab_x[k] * fh_144[k]
                   + fi_186[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, fh_145, fh_146, fh_147, \
                         fh_148, fh_149, fi_187, fi_188, fi_196, fi_197, \
                         fi_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_x[k] * fh_145[k]
                   + fi_187[k];

        t_146[k] = -ab_x[k] * fh_146[k]
                   + fi_188[k];

        t_147[k] = -ab_x[k] * fh_147[k]
                   + fi_196[k];

        t_148[k] = -ab_x[k] * fh_148[k]
                   + fi_197[k];

        t_149[k] = -ab_x[k] * fh_149[k]
                   + fi_198[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, fh_150, fh_151, fh_152, \
                         fh_153, fh_154, fi_199, fi_200, fi_201, fi_202, \
                         fi_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_x[k] * fh_150[k]
                   + fi_199[k];

        t_151[k] = -ab_x[k] * fh_151[k]
                   + fi_200[k];

        t_152[k] = -ab_x[k] * fh_152[k]
                   + fi_201[k];

        t_153[k] = -ab_x[k] * fh_153[k]
                   + fi_202[k];

        t_154[k] = -ab_x[k] * fh_154[k]
                   + fi_203[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, fh_155, fh_156, fh_157, \
                         fh_158, fh_159, fi_204, fi_205, fi_206, fi_207, \
                         fi_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_x[k] * fh_155[k]
                   + fi_204[k];

        t_156[k] = -ab_x[k] * fh_156[k]
                   + fi_205[k];

        t_157[k] = -ab_x[k] * fh_157[k]
                   + fi_206[k];

        t_158[k] = -ab_x[k] * fh_158[k]
                   + fi_207[k];

        t_159[k] = -ab_x[k] * fh_159[k]
                   + fi_208[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, fh_160, fh_161, fh_162, \
                         fh_163, fh_164, fi_209, fi_210, fi_211, fi_212, \
                         fi_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_x[k] * fh_160[k]
                   + fi_209[k];

        t_161[k] = -ab_x[k] * fh_161[k]
                   + fi_210[k];

        t_162[k] = -ab_x[k] * fh_162[k]
                   + fi_211[k];

        t_163[k] = -ab_x[k] * fh_163[k]
                   + fi_212[k];

        t_164[k] = -ab_x[k] * fh_164[k]
                   + fi_213[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, fh_165, fh_166, fh_167, \
                         fh_168, fh_169, fi_214, fi_215, fi_216, fi_224, \
                         fi_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = -ab_x[k] * fh_165[k]
                   + fi_214[k];

        t_166[k] = -ab_x[k] * fh_166[k]
                   + fi_215[k];

        t_167[k] = -ab_x[k] * fh_167[k]
                   + fi_216[k];

        t_168[k] = -ab_x[k] * fh_168[k]
                   + fi_224[k];

        t_169[k] = -ab_x[k] * fh_169[k]
                   + fi_225[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, fh_170, fh_171, fh_172, \
                         fh_173, fh_174, fi_226, fi_227, fi_228, fi_229, \
                         fi_230 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = -ab_x[k] * fh_170[k]
                   + fi_226[k];

        t_171[k] = -ab_x[k] * fh_171[k]
                   + fi_227[k];

        t_172[k] = -ab_x[k] * fh_172[k]
                   + fi_228[k];

        t_173[k] = -ab_x[k] * fh_173[k]
                   + fi_229[k];

        t_174[k] = -ab_x[k] * fh_174[k]
                   + fi_230[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, fh_175, fh_176, fh_177, \
                         fh_178, fh_179, fi_231, fi_232, fi_233, fi_234, \
                         fi_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = -ab_x[k] * fh_175[k]
                   + fi_231[k];

        t_176[k] = -ab_x[k] * fh_176[k]
                   + fi_232[k];

        t_177[k] = -ab_x[k] * fh_177[k]
                   + fi_233[k];

        t_178[k] = -ab_x[k] * fh_178[k]
                   + fi_234[k];

        t_179[k] = -ab_x[k] * fh_179[k]
                   + fi_235[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, fh_180, fh_181, fh_182, \
                         fh_183, fh_184, fi_236, fi_237, fi_238, fi_239, \
                         fi_240 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = -ab_x[k] * fh_180[k]
                   + fi_236[k];

        t_181[k] = -ab_x[k] * fh_181[k]
                   + fi_237[k];

        t_182[k] = -ab_x[k] * fh_182[k]
                   + fi_238[k];

        t_183[k] = -ab_x[k] * fh_183[k]
                   + fi_239[k];

        t_184[k] = -ab_x[k] * fh_184[k]
                   + fi_240[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, fh_185, fh_186, fh_187, \
                         fh_188, fh_189, fi_241, fi_242, fi_243, fi_244, \
                         fi_252 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = -ab_x[k] * fh_185[k]
                   + fi_241[k];

        t_186[k] = -ab_x[k] * fh_186[k]
                   + fi_242[k];

        t_187[k] = -ab_x[k] * fh_187[k]
                   + fi_243[k];

        t_188[k] = -ab_x[k] * fh_188[k]
                   + fi_244[k];

        t_189[k] = -ab_x[k] * fh_189[k]
                   + fi_252[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, fh_190, fh_191, fh_192, \
                         fh_193, fh_194, fi_253, fi_254, fi_255, fi_256, \
                         fi_257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = -ab_x[k] * fh_190[k]
                   + fi_253[k];

        t_191[k] = -ab_x[k] * fh_191[k]
                   + fi_254[k];

        t_192[k] = -ab_x[k] * fh_192[k]
                   + fi_255[k];

        t_193[k] = -ab_x[k] * fh_193[k]
                   + fi_256[k];

        t_194[k] = -ab_x[k] * fh_194[k]
                   + fi_257[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, fh_195, fh_196, fh_197, \
                         fh_198, fh_199, fi_258, fi_259, fi_260, fi_261, \
                         fi_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = -ab_x[k] * fh_195[k]
                   + fi_258[k];

        t_196[k] = -ab_x[k] * fh_196[k]
                   + fi_259[k];

        t_197[k] = -ab_x[k] * fh_197[k]
                   + fi_260[k];

        t_198[k] = -ab_x[k] * fh_198[k]
                   + fi_261[k];

        t_199[k] = -ab_x[k] * fh_199[k]
                   + fi_262[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, fh_200, fh_201, fh_202, \
                         fh_203, fh_204, fi_263, fi_264, fi_265, fi_266, \
                         fi_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = -ab_x[k] * fh_200[k]
                   + fi_263[k];

        t_201[k] = -ab_x[k] * fh_201[k]
                   + fi_264[k];

        t_202[k] = -ab_x[k] * fh_202[k]
                   + fi_265[k];

        t_203[k] = -ab_x[k] * fh_203[k]
                   + fi_266[k];

        t_204[k] = -ab_x[k] * fh_204[k]
                   + fi_267[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, fh_205, fh_206, fh_207, \
                         fh_208, fh_209, fi_268, fi_269, fi_270, fi_271, \
                         fi_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = -ab_x[k] * fh_205[k]
                   + fi_268[k];

        t_206[k] = -ab_x[k] * fh_206[k]
                   + fi_269[k];

        t_207[k] = -ab_x[k] * fh_207[k]
                   + fi_270[k];

        t_208[k] = -ab_x[k] * fh_208[k]
                   + fi_271[k];

        t_209[k] = -ab_x[k] * fh_209[k]
                   + fi_272[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_y, fh_126, fh_127, fh_128, \
                         fh_129, fh_130, fi_169, fi_171, fi_172, fi_174, \
                         fi_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = -ab_y[k] * fh_126[k]
                   + fi_169[k];

        t_211[k] = -ab_y[k] * fh_127[k]
                   + fi_171[k];

        t_212[k] = -ab_y[k] * fh_128[k]
                   + fi_172[k];

        t_213[k] = -ab_y[k] * fh_129[k]
                   + fi_174[k];

        t_214[k] = -ab_y[k] * fh_130[k]
                   + fi_175[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_y, fh_131, fh_132, fh_133, \
                         fh_134, fh_135, fi_176, fi_178, fi_179, fi_180, \
                         fi_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = -ab_y[k] * fh_131[k]
                   + fi_176[k];

        t_216[k] = -ab_y[k] * fh_132[k]
                   + fi_178[k];

        t_217[k] = -ab_y[k] * fh_133[k]
                   + fi_179[k];

        t_218[k] = -ab_y[k] * fh_134[k]
                   + fi_180[k];

        t_219[k] = -ab_y[k] * fh_135[k]
                   + fi_181[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, fh_136, fh_137, fh_138, \
                         fh_139, fh_140, fi_183, fi_184, fi_185, fi_186, \
                         fi_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = -ab_y[k] * fh_136[k]
                   + fi_183[k];

        t_221[k] = -ab_y[k] * fh_137[k]
                   + fi_184[k];

        t_222[k] = -ab_y[k] * fh_138[k]
                   + fi_185[k];

        t_223[k] = -ab_y[k] * fh_139[k]
                   + fi_186[k];

        t_224[k] = -ab_y[k] * fh_140[k]
                   + fi_187[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_y, fh_141, fh_142, fh_143, \
                         fh_144, fh_145, fi_189, fi_190, fi_191, fi_192, \
                         fi_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = -ab_y[k] * fh_141[k]
                   + fi_189[k];

        t_226[k] = -ab_y[k] * fh_142[k]
                   + fi_190[k];

        t_227[k] = -ab_y[k] * fh_143[k]
                   + fi_191[k];

        t_228[k] = -ab_y[k] * fh_144[k]
                   + fi_192[k];

        t_229[k] = -ab_y[k] * fh_145[k]
                   + fi_193[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_y, fh_146, fh_147, fh_148, \
                         fh_149, fh_150, fi_194, fi_197, fi_199, fi_200, \
                         fi_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = -ab_y[k] * fh_146[k]
                   + fi_194[k];

        t_231[k] = -ab_y[k] * fh_147[k]
                   + fi_197[k];

        t_232[k] = -ab_y[k] * fh_148[k]
                   + fi_199[k];

        t_233[k] = -ab_y[k] * fh_149[k]
                   + fi_200[k];

        t_234[k] = -ab_y[k] * fh_150[k]
                   + fi_202[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, fh_151, fh_152, fh_153, \
                         fh_154, fh_155, fi_203, fi_204, fi_206, fi_207, \
                         fi_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = -ab_y[k] * fh_151[k]
                   + fi_203[k];

        t_236[k] = -ab_y[k] * fh_152[k]
                   + fi_204[k];

        t_237[k] = -ab_y[k] * fh_153[k]
                   + fi_206[k];

        t_238[k] = -ab_y[k] * fh_154[k]
                   + fi_207[k];

        t_239[k] = -ab_y[k] * fh_155[k]
                   + fi_208[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_y, fh_156, fh_157, fh_158, \
                         fh_159, fh_160, fi_209, fi_211, fi_212, fi_213, \
                         fi_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = -ab_y[k] * fh_156[k]
                   + fi_209[k];

        t_241[k] = -ab_y[k] * fh_157[k]
                   + fi_211[k];

        t_242[k] = -ab_y[k] * fh_158[k]
                   + fi_212[k];

        t_243[k] = -ab_y[k] * fh_159[k]
                   + fi_213[k];

        t_244[k] = -ab_y[k] * fh_160[k]
                   + fi_214[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_y, fh_161, fh_162, fh_163, \
                         fh_164, fh_165, fi_215, fi_217, fi_218, fi_219, \
                         fi_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = -ab_y[k] * fh_161[k]
                   + fi_215[k];

        t_246[k] = -ab_y[k] * fh_162[k]
                   + fi_217[k];

        t_247[k] = -ab_y[k] * fh_163[k]
                   + fi_218[k];

        t_248[k] = -ab_y[k] * fh_164[k]
                   + fi_219[k];

        t_249[k] = -ab_y[k] * fh_165[k]
                   + fi_220[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, fh_166, fh_167, fh_168, \
                         fh_169, fh_170, fi_221, fi_222, fi_225, fi_227, \
                         fi_228 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = -ab_y[k] * fh_166[k]
                   + fi_221[k];

        t_251[k] = -ab_y[k] * fh_167[k]
                   + fi_222[k];

        t_252[k] = -ab_y[k] * fh_168[k]
                   + fi_225[k];

        t_253[k] = -ab_y[k] * fh_169[k]
                   + fi_227[k];

        t_254[k] = -ab_y[k] * fh_170[k]
                   + fi_228[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_y, fh_171, fh_172, fh_173, \
                         fh_174, fh_175, fi_230, fi_231, fi_232, fi_234, \
                         fi_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = -ab_y[k] * fh_171[k]
                   + fi_230[k];

        t_256[k] = -ab_y[k] * fh_172[k]
                   + fi_231[k];

        t_257[k] = -ab_y[k] * fh_173[k]
                   + fi_232[k];

        t_258[k] = -ab_y[k] * fh_174[k]
                   + fi_234[k];

        t_259[k] = -ab_y[k] * fh_175[k]
                   + fi_235[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_y, fh_176, fh_177, fh_178, \
                         fh_179, fh_180, fi_236, fi_237, fi_239, fi_240, \
                         fi_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = -ab_y[k] * fh_176[k]
                   + fi_236[k];

        t_261[k] = -ab_y[k] * fh_177[k]
                   + fi_237[k];

        t_262[k] = -ab_y[k] * fh_178[k]
                   + fi_239[k];

        t_263[k] = -ab_y[k] * fh_179[k]
                   + fi_240[k];

        t_264[k] = -ab_y[k] * fh_180[k]
                   + fi_241[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, fh_181, fh_182, fh_183, \
                         fh_184, fh_185, fi_242, fi_243, fi_245, fi_246, \
                         fi_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = -ab_y[k] * fh_181[k]
                   + fi_242[k];

        t_266[k] = -ab_y[k] * fh_182[k]
                   + fi_243[k];

        t_267[k] = -ab_y[k] * fh_183[k]
                   + fi_245[k];

        t_268[k] = -ab_y[k] * fh_184[k]
                   + fi_246[k];

        t_269[k] = -ab_y[k] * fh_185[k]
                   + fi_247[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_y, fh_186, fh_187, fh_188, \
                         fh_189, fh_190, fi_248, fi_249, fi_250, fi_253, \
                         fi_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = -ab_y[k] * fh_186[k]
                   + fi_248[k];

        t_271[k] = -ab_y[k] * fh_187[k]
                   + fi_249[k];

        t_272[k] = -ab_y[k] * fh_188[k]
                   + fi_250[k];

        t_273[k] = -ab_y[k] * fh_189[k]
                   + fi_253[k];

        t_274[k] = -ab_y[k] * fh_190[k]
                   + fi_255[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_y, fh_191, fh_192, fh_193, \
                         fh_194, fh_195, fi_256, fi_258, fi_259, fi_260, \
                         fi_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = -ab_y[k] * fh_191[k]
                   + fi_256[k];

        t_276[k] = -ab_y[k] * fh_192[k]
                   + fi_258[k];

        t_277[k] = -ab_y[k] * fh_193[k]
                   + fi_259[k];

        t_278[k] = -ab_y[k] * fh_194[k]
                   + fi_260[k];

        t_279[k] = -ab_y[k] * fh_195[k]
                   + fi_262[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, fh_196, fh_197, fh_198, \
                         fh_199, fh_200, fi_263, fi_264, fi_265, fi_267, \
                         fi_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_280[k] = -ab_y[k] * fh_196[k]
                   + fi_263[k];

        t_281[k] = -ab_y[k] * fh_197[k]
                   + fi_264[k];

        t_282[k] = -ab_y[k] * fh_198[k]
                   + fi_265[k];

        t_283[k] = -ab_y[k] * fh_199[k]
                   + fi_267[k];

        t_284[k] = -ab_y[k] * fh_200[k]
                   + fi_268[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_y, fh_201, fh_202, fh_203, \
                         fh_204, fh_205, fi_269, fi_270, fi_271, fi_273, \
                         fi_274 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_285[k] = -ab_y[k] * fh_201[k]
                   + fi_269[k];

        t_286[k] = -ab_y[k] * fh_202[k]
                   + fi_270[k];

        t_287[k] = -ab_y[k] * fh_203[k]
                   + fi_271[k];

        t_288[k] = -ab_y[k] * fh_204[k]
                   + fi_273[k];

        t_289[k] = -ab_y[k] * fh_205[k]
                   + fi_274[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, ab_y, fh_206, fh_207, fh_208, fh_209, \
                         fi_275, fi_276, fi_277, fi_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_290[k] = -ab_y[k] * fh_206[k]
                   + fi_275[k];

        t_291[k] = -ab_y[k] * fh_207[k]
                   + fi_276[k];

        t_292[k] = -ab_y[k] * fh_208[k]
                   + fi_277[k];

        t_293[k] = -ab_y[k] * fh_209[k]
                   + fi_278[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_z, fh_189, fh_190, fh_191, \
                         fh_192, fh_193, fi_254, fi_256, fi_257, fi_259, \
                         fi_260 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_294[k] = -ab_z[k] * fh_189[k]
                   + fi_254[k];

        t_295[k] = -ab_z[k] * fh_190[k]
                   + fi_256[k];

        t_296[k] = -ab_z[k] * fh_191[k]
                   + fi_257[k];

        t_297[k] = -ab_z[k] * fh_192[k]
                   + fi_259[k];

        t_298[k] = -ab_z[k] * fh_193[k]
                   + fi_260[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_z, fh_194, fh_195, fh_196, \
                         fh_197, fh_198, fi_261, fi_263, fi_264, fi_265, \
                         fi_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_299[k] = -ab_z[k] * fh_194[k]
                   + fi_261[k];

        t_300[k] = -ab_z[k] * fh_195[k]
                   + fi_263[k];

        t_301[k] = -ab_z[k] * fh_196[k]
                   + fi_264[k];

        t_302[k] = -ab_z[k] * fh_197[k]
                   + fi_265[k];

        t_303[k] = -ab_z[k] * fh_198[k]
                   + fi_266[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_z, fh_199, fh_200, fh_201, \
                         fh_202, fh_203, fi_268, fi_269, fi_270, fi_271, \
                         fi_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_304[k] = -ab_z[k] * fh_199[k]
                   + fi_268[k];

        t_305[k] = -ab_z[k] * fh_200[k]
                   + fi_269[k];

        t_306[k] = -ab_z[k] * fh_201[k]
                   + fi_270[k];

        t_307[k] = -ab_z[k] * fh_202[k]
                   + fi_271[k];

        t_308[k] = -ab_z[k] * fh_203[k]
                   + fi_272[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_z, fh_204, fh_205, fh_206, \
                         fh_207, fh_208, fi_274, fi_275, fi_276, fi_277, \
                         fi_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_309[k] = -ab_z[k] * fh_204[k]
                   + fi_274[k];

        t_310[k] = -ab_z[k] * fh_205[k]
                   + fi_275[k];

        t_311[k] = -ab_z[k] * fh_206[k]
                   + fi_276[k];

        t_312[k] = -ab_z[k] * fh_207[k]
                   + fi_277[k];

        t_313[k] = -ab_z[k] * fh_208[k]
                   + fi_278[k];
    }

#pragma omp simd aligned(t_314, ab_z, fh_209, fi_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_314[k] = -ab_z[k] * fh_209[k]
                   + fi_279[k];
    }
}

}  // namespace simdtrf
