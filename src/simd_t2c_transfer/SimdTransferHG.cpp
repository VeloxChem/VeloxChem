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


#include "SimdTransferHG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t hf, const size_t if_, const size_t nmax) -> void
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
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_152 = buffer.data(hf + 152);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_154 = buffer.data(hf + 154);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_161 = buffer.data(hf + 161);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_163 = buffer.data(hf + 163);
    const auto *hf_164 = buffer.data(hf + 164);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_191 = buffer.data(hf + 191);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_194 = buffer.data(hf + 194);
    const auto *hf_195 = buffer.data(hf + 195);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_204 = buffer.data(hf + 204);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__4 = buffer.data(if_ + 4);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__143 = buffer.data(if_ + 143);
    const auto *if__144 = buffer.data(if_ + 144);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__152 = buffer.data(if_ + 152);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__154 = buffer.data(if_ + 154);
    const auto *if__155 = buffer.data(if_ + 155);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__157 = buffer.data(if_ + 157);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__201 = buffer.data(if_ + 201);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__203 = buffer.data(if_ + 203);
    const auto *if__204 = buffer.data(if_ + 204);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__208 = buffer.data(if_ + 208);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, hf_0, hf_1, hf_2, hf_3, hf_4, if__0, \
                         if__1, if__2, if__3, if__4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * hf_0[k]
                 + if__0[k];

        t_1[k] = ab_x[k] * hf_1[k]
                 + if__1[k];

        t_2[k] = ab_x[k] * hf_2[k]
                 + if__2[k];

        t_3[k] = ab_x[k] * hf_3[k]
                 + if__3[k];

        t_4[k] = ab_x[k] * hf_4[k]
                 + if__4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, hf_5, hf_6, hf_7, hf_8, hf_9, if__5, \
                         if__6, if__7, if__8, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_x[k] * hf_5[k]
                 + if__5[k];

        t_6[k] = ab_x[k] * hf_6[k]
                 + if__6[k];

        t_7[k] = ab_x[k] * hf_7[k]
                 + if__7[k];

        t_8[k] = ab_x[k] * hf_8[k]
                 + if__8[k];

        t_9[k] = ab_x[k] * hf_9[k]
                 + if__9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, hf_6, hf_7, hf_8, hf_9, \
                         if__16, if__17, if__18, if__19, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = ab_y[k] * hf_6[k]
                  + if__16[k];

        t_11[k] = ab_y[k] * hf_7[k]
                  + if__17[k];

        t_12[k] = ab_y[k] * hf_8[k]
                  + if__18[k];

        t_13[k] = ab_y[k] * hf_9[k]
                  + if__19[k];

        t_14[k] = ab_z[k] * hf_9[k]
                  + if__29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, hf_10, hf_11, hf_12, hf_13, \
                         hf_14, if__10, if__11, if__12, if__13, \
                         if__14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = ab_x[k] * hf_10[k]
                  + if__10[k];

        t_16[k] = ab_x[k] * hf_11[k]
                  + if__11[k];

        t_17[k] = ab_x[k] * hf_12[k]
                  + if__12[k];

        t_18[k] = ab_x[k] * hf_13[k]
                  + if__13[k];

        t_19[k] = ab_x[k] * hf_14[k]
                  + if__14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, hf_15, hf_16, hf_17, hf_18, \
                         hf_19, if__15, if__16, if__17, if__18, \
                         if__19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = ab_x[k] * hf_15[k]
                  + if__15[k];

        t_21[k] = ab_x[k] * hf_16[k]
                  + if__16[k];

        t_22[k] = ab_x[k] * hf_17[k]
                  + if__17[k];

        t_23[k] = ab_x[k] * hf_18[k]
                  + if__18[k];

        t_24[k] = ab_x[k] * hf_19[k]
                  + if__19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, hf_16, hf_17, hf_18, hf_19, \
                         if__36, if__37, if__38, if__39, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = ab_y[k] * hf_16[k]
                  + if__36[k];

        t_26[k] = ab_y[k] * hf_17[k]
                  + if__37[k];

        t_27[k] = ab_y[k] * hf_18[k]
                  + if__38[k];

        t_28[k] = ab_y[k] * hf_19[k]
                  + if__39[k];

        t_29[k] = ab_z[k] * hf_19[k]
                  + if__49[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, hf_20, hf_21, hf_22, hf_23, \
                         hf_24, if__20, if__21, if__22, if__23, \
                         if__24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = ab_x[k] * hf_20[k]
                  + if__20[k];

        t_31[k] = ab_x[k] * hf_21[k]
                  + if__21[k];

        t_32[k] = ab_x[k] * hf_22[k]
                  + if__22[k];

        t_33[k] = ab_x[k] * hf_23[k]
                  + if__23[k];

        t_34[k] = ab_x[k] * hf_24[k]
                  + if__24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, hf_25, hf_26, hf_27, hf_28, \
                         hf_29, if__25, if__26, if__27, if__28, \
                         if__29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = ab_x[k] * hf_25[k]
                  + if__25[k];

        t_36[k] = ab_x[k] * hf_26[k]
                  + if__26[k];

        t_37[k] = ab_x[k] * hf_27[k]
                  + if__27[k];

        t_38[k] = ab_x[k] * hf_28[k]
                  + if__28[k];

        t_39[k] = ab_x[k] * hf_29[k]
                  + if__29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, hf_26, hf_27, hf_28, hf_29, \
                         if__46, if__47, if__48, if__49, if__59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = ab_y[k] * hf_26[k]
                  + if__46[k];

        t_41[k] = ab_y[k] * hf_27[k]
                  + if__47[k];

        t_42[k] = ab_y[k] * hf_28[k]
                  + if__48[k];

        t_43[k] = ab_y[k] * hf_29[k]
                  + if__49[k];

        t_44[k] = ab_z[k] * hf_29[k]
                  + if__59[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, hf_30, hf_31, hf_32, hf_33, \
                         hf_34, if__30, if__31, if__32, if__33, \
                         if__34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_x[k] * hf_30[k]
                  + if__30[k];

        t_46[k] = ab_x[k] * hf_31[k]
                  + if__31[k];

        t_47[k] = ab_x[k] * hf_32[k]
                  + if__32[k];

        t_48[k] = ab_x[k] * hf_33[k]
                  + if__33[k];

        t_49[k] = ab_x[k] * hf_34[k]
                  + if__34[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, hf_35, hf_36, hf_37, hf_38, \
                         hf_39, if__35, if__36, if__37, if__38, \
                         if__39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = ab_x[k] * hf_35[k]
                  + if__35[k];

        t_51[k] = ab_x[k] * hf_36[k]
                  + if__36[k];

        t_52[k] = ab_x[k] * hf_37[k]
                  + if__37[k];

        t_53[k] = ab_x[k] * hf_38[k]
                  + if__38[k];

        t_54[k] = ab_x[k] * hf_39[k]
                  + if__39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, hf_36, hf_37, hf_38, hf_39, \
                         if__66, if__67, if__68, if__69, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = ab_y[k] * hf_36[k]
                  + if__66[k];

        t_56[k] = ab_y[k] * hf_37[k]
                  + if__67[k];

        t_57[k] = ab_y[k] * hf_38[k]
                  + if__68[k];

        t_58[k] = ab_y[k] * hf_39[k]
                  + if__69[k];

        t_59[k] = ab_z[k] * hf_39[k]
                  + if__79[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, hf_40, hf_41, hf_42, hf_43, \
                         hf_44, if__40, if__41, if__42, if__43, \
                         if__44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * hf_40[k]
                  + if__40[k];

        t_61[k] = ab_x[k] * hf_41[k]
                  + if__41[k];

        t_62[k] = ab_x[k] * hf_42[k]
                  + if__42[k];

        t_63[k] = ab_x[k] * hf_43[k]
                  + if__43[k];

        t_64[k] = ab_x[k] * hf_44[k]
                  + if__44[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, hf_45, hf_46, hf_47, hf_48, \
                         hf_49, if__45, if__46, if__47, if__48, \
                         if__49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_x[k] * hf_45[k]
                  + if__45[k];

        t_66[k] = ab_x[k] * hf_46[k]
                  + if__46[k];

        t_67[k] = ab_x[k] * hf_47[k]
                  + if__47[k];

        t_68[k] = ab_x[k] * hf_48[k]
                  + if__48[k];

        t_69[k] = ab_x[k] * hf_49[k]
                  + if__49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, hf_46, hf_47, hf_48, hf_49, \
                         if__76, if__77, if__78, if__79, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = ab_y[k] * hf_46[k]
                  + if__76[k];

        t_71[k] = ab_y[k] * hf_47[k]
                  + if__77[k];

        t_72[k] = ab_y[k] * hf_48[k]
                  + if__78[k];

        t_73[k] = ab_y[k] * hf_49[k]
                  + if__79[k];

        t_74[k] = ab_z[k] * hf_49[k]
                  + if__89[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, hf_50, hf_51, hf_52, hf_53, \
                         hf_54, if__50, if__51, if__52, if__53, \
                         if__54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = ab_x[k] * hf_50[k]
                  + if__50[k];

        t_76[k] = ab_x[k] * hf_51[k]
                  + if__51[k];

        t_77[k] = ab_x[k] * hf_52[k]
                  + if__52[k];

        t_78[k] = ab_x[k] * hf_53[k]
                  + if__53[k];

        t_79[k] = ab_x[k] * hf_54[k]
                  + if__54[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, hf_55, hf_56, hf_57, hf_58, \
                         hf_59, if__55, if__56, if__57, if__58, \
                         if__59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = ab_x[k] * hf_55[k]
                  + if__55[k];

        t_81[k] = ab_x[k] * hf_56[k]
                  + if__56[k];

        t_82[k] = ab_x[k] * hf_57[k]
                  + if__57[k];

        t_83[k] = ab_x[k] * hf_58[k]
                  + if__58[k];

        t_84[k] = ab_x[k] * hf_59[k]
                  + if__59[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, hf_56, hf_57, hf_58, hf_59, \
                         if__86, if__87, if__88, if__89, if__99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_y[k] * hf_56[k]
                  + if__86[k];

        t_86[k] = ab_y[k] * hf_57[k]
                  + if__87[k];

        t_87[k] = ab_y[k] * hf_58[k]
                  + if__88[k];

        t_88[k] = ab_y[k] * hf_59[k]
                  + if__89[k];

        t_89[k] = ab_z[k] * hf_59[k]
                  + if__99[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, hf_60, hf_61, hf_62, hf_63, \
                         hf_64, if__60, if__61, if__62, if__63, \
                         if__64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * hf_60[k]
                  + if__60[k];

        t_91[k] = ab_x[k] * hf_61[k]
                  + if__61[k];

        t_92[k] = ab_x[k] * hf_62[k]
                  + if__62[k];

        t_93[k] = ab_x[k] * hf_63[k]
                  + if__63[k];

        t_94[k] = ab_x[k] * hf_64[k]
                  + if__64[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, hf_65, hf_66, hf_67, hf_68, \
                         hf_69, if__65, if__66, if__67, if__68, \
                         if__69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_x[k] * hf_65[k]
                  + if__65[k];

        t_96[k] = ab_x[k] * hf_66[k]
                  + if__66[k];

        t_97[k] = ab_x[k] * hf_67[k]
                  + if__67[k];

        t_98[k] = ab_x[k] * hf_68[k]
                  + if__68[k];

        t_99[k] = ab_x[k] * hf_69[k]
                  + if__69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, hf_66, hf_67, hf_68, \
                         hf_69, if__106, if__107, if__108, if__109, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = ab_y[k] * hf_66[k]
                   + if__106[k];

        t_101[k] = ab_y[k] * hf_67[k]
                   + if__107[k];

        t_102[k] = ab_y[k] * hf_68[k]
                   + if__108[k];

        t_103[k] = ab_y[k] * hf_69[k]
                   + if__109[k];

        t_104[k] = ab_z[k] * hf_69[k]
                   + if__119[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, hf_70, hf_71, hf_72, hf_73, \
                         hf_74, if__70, if__71, if__72, if__73, \
                         if__74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = ab_x[k] * hf_70[k]
                   + if__70[k];

        t_106[k] = ab_x[k] * hf_71[k]
                   + if__71[k];

        t_107[k] = ab_x[k] * hf_72[k]
                   + if__72[k];

        t_108[k] = ab_x[k] * hf_73[k]
                   + if__73[k];

        t_109[k] = ab_x[k] * hf_74[k]
                   + if__74[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, hf_75, hf_76, hf_77, hf_78, \
                         hf_79, if__75, if__76, if__77, if__78, \
                         if__79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = ab_x[k] * hf_75[k]
                   + if__75[k];

        t_111[k] = ab_x[k] * hf_76[k]
                   + if__76[k];

        t_112[k] = ab_x[k] * hf_77[k]
                   + if__77[k];

        t_113[k] = ab_x[k] * hf_78[k]
                   + if__78[k];

        t_114[k] = ab_x[k] * hf_79[k]
                   + if__79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, hf_76, hf_77, hf_78, \
                         hf_79, if__116, if__117, if__118, if__119, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = ab_y[k] * hf_76[k]
                   + if__116[k];

        t_116[k] = ab_y[k] * hf_77[k]
                   + if__117[k];

        t_117[k] = ab_y[k] * hf_78[k]
                   + if__118[k];

        t_118[k] = ab_y[k] * hf_79[k]
                   + if__119[k];

        t_119[k] = ab_z[k] * hf_79[k]
                   + if__129[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, hf_80, hf_81, hf_82, hf_83, \
                         hf_84, if__80, if__81, if__82, if__83, \
                         if__84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = ab_x[k] * hf_80[k]
                   + if__80[k];

        t_121[k] = ab_x[k] * hf_81[k]
                   + if__81[k];

        t_122[k] = ab_x[k] * hf_82[k]
                   + if__82[k];

        t_123[k] = ab_x[k] * hf_83[k]
                   + if__83[k];

        t_124[k] = ab_x[k] * hf_84[k]
                   + if__84[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, hf_85, hf_86, hf_87, hf_88, \
                         hf_89, if__85, if__86, if__87, if__88, \
                         if__89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = ab_x[k] * hf_85[k]
                   + if__85[k];

        t_126[k] = ab_x[k] * hf_86[k]
                   + if__86[k];

        t_127[k] = ab_x[k] * hf_87[k]
                   + if__87[k];

        t_128[k] = ab_x[k] * hf_88[k]
                   + if__88[k];

        t_129[k] = ab_x[k] * hf_89[k]
                   + if__89[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, hf_86, hf_87, hf_88, \
                         hf_89, if__126, if__127, if__128, if__129, \
                         if__139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = ab_y[k] * hf_86[k]
                   + if__126[k];

        t_131[k] = ab_y[k] * hf_87[k]
                   + if__127[k];

        t_132[k] = ab_y[k] * hf_88[k]
                   + if__128[k];

        t_133[k] = ab_y[k] * hf_89[k]
                   + if__129[k];

        t_134[k] = ab_z[k] * hf_89[k]
                   + if__139[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, hf_90, hf_91, hf_92, hf_93, \
                         hf_94, if__90, if__91, if__92, if__93, \
                         if__94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_x[k] * hf_90[k]
                   + if__90[k];

        t_136[k] = ab_x[k] * hf_91[k]
                   + if__91[k];

        t_137[k] = ab_x[k] * hf_92[k]
                   + if__92[k];

        t_138[k] = ab_x[k] * hf_93[k]
                   + if__93[k];

        t_139[k] = ab_x[k] * hf_94[k]
                   + if__94[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, hf_95, hf_96, hf_97, hf_98, \
                         hf_99, if__95, if__96, if__97, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = ab_x[k] * hf_95[k]
                   + if__95[k];

        t_141[k] = ab_x[k] * hf_96[k]
                   + if__96[k];

        t_142[k] = ab_x[k] * hf_97[k]
                   + if__97[k];

        t_143[k] = ab_x[k] * hf_98[k]
                   + if__98[k];

        t_144[k] = ab_x[k] * hf_99[k]
                   + if__99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, hf_96, hf_97, hf_98, \
                         hf_99, if__136, if__137, if__138, if__139, \
                         if__149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = ab_y[k] * hf_96[k]
                   + if__136[k];

        t_146[k] = ab_y[k] * hf_97[k]
                   + if__137[k];

        t_147[k] = ab_y[k] * hf_98[k]
                   + if__138[k];

        t_148[k] = ab_y[k] * hf_99[k]
                   + if__139[k];

        t_149[k] = ab_z[k] * hf_99[k]
                   + if__149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, hf_100, hf_101, hf_102, \
                         hf_103, hf_104, if__100, if__101, if__102, if__103, \
                         if__104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = ab_x[k] * hf_100[k]
                   + if__100[k];

        t_151[k] = ab_x[k] * hf_101[k]
                   + if__101[k];

        t_152[k] = ab_x[k] * hf_102[k]
                   + if__102[k];

        t_153[k] = ab_x[k] * hf_103[k]
                   + if__103[k];

        t_154[k] = ab_x[k] * hf_104[k]
                   + if__104[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, hf_105, hf_106, hf_107, \
                         hf_108, hf_109, if__105, if__106, if__107, if__108, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = ab_x[k] * hf_105[k]
                   + if__105[k];

        t_156[k] = ab_x[k] * hf_106[k]
                   + if__106[k];

        t_157[k] = ab_x[k] * hf_107[k]
                   + if__107[k];

        t_158[k] = ab_x[k] * hf_108[k]
                   + if__108[k];

        t_159[k] = ab_x[k] * hf_109[k]
                   + if__109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, hf_106, hf_107, \
                         hf_108, hf_109, if__156, if__157, if__158, if__159, \
                         if__169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = ab_y[k] * hf_106[k]
                   + if__156[k];

        t_161[k] = ab_y[k] * hf_107[k]
                   + if__157[k];

        t_162[k] = ab_y[k] * hf_108[k]
                   + if__158[k];

        t_163[k] = ab_y[k] * hf_109[k]
                   + if__159[k];

        t_164[k] = ab_z[k] * hf_109[k]
                   + if__169[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, hf_110, hf_111, hf_112, \
                         hf_113, hf_114, if__110, if__111, if__112, if__113, \
                         if__114 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = ab_x[k] * hf_110[k]
                   + if__110[k];

        t_166[k] = ab_x[k] * hf_111[k]
                   + if__111[k];

        t_167[k] = ab_x[k] * hf_112[k]
                   + if__112[k];

        t_168[k] = ab_x[k] * hf_113[k]
                   + if__113[k];

        t_169[k] = ab_x[k] * hf_114[k]
                   + if__114[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, hf_115, hf_116, hf_117, \
                         hf_118, hf_119, if__115, if__116, if__117, if__118, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = ab_x[k] * hf_115[k]
                   + if__115[k];

        t_171[k] = ab_x[k] * hf_116[k]
                   + if__116[k];

        t_172[k] = ab_x[k] * hf_117[k]
                   + if__117[k];

        t_173[k] = ab_x[k] * hf_118[k]
                   + if__118[k];

        t_174[k] = ab_x[k] * hf_119[k]
                   + if__119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, hf_116, hf_117, \
                         hf_118, hf_119, if__166, if__167, if__168, if__169, \
                         if__179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = ab_y[k] * hf_116[k]
                   + if__166[k];

        t_176[k] = ab_y[k] * hf_117[k]
                   + if__167[k];

        t_177[k] = ab_y[k] * hf_118[k]
                   + if__168[k];

        t_178[k] = ab_y[k] * hf_119[k]
                   + if__169[k];

        t_179[k] = ab_z[k] * hf_119[k]
                   + if__179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, hf_120, hf_121, hf_122, \
                         hf_123, hf_124, if__120, if__121, if__122, if__123, \
                         if__124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = ab_x[k] * hf_120[k]
                   + if__120[k];

        t_181[k] = ab_x[k] * hf_121[k]
                   + if__121[k];

        t_182[k] = ab_x[k] * hf_122[k]
                   + if__122[k];

        t_183[k] = ab_x[k] * hf_123[k]
                   + if__123[k];

        t_184[k] = ab_x[k] * hf_124[k]
                   + if__124[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, hf_125, hf_126, hf_127, \
                         hf_128, hf_129, if__125, if__126, if__127, if__128, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = ab_x[k] * hf_125[k]
                   + if__125[k];

        t_186[k] = ab_x[k] * hf_126[k]
                   + if__126[k];

        t_187[k] = ab_x[k] * hf_127[k]
                   + if__127[k];

        t_188[k] = ab_x[k] * hf_128[k]
                   + if__128[k];

        t_189[k] = ab_x[k] * hf_129[k]
                   + if__129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, hf_126, hf_127, \
                         hf_128, hf_129, if__176, if__177, if__178, if__179, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = ab_y[k] * hf_126[k]
                   + if__176[k];

        t_191[k] = ab_y[k] * hf_127[k]
                   + if__177[k];

        t_192[k] = ab_y[k] * hf_128[k]
                   + if__178[k];

        t_193[k] = ab_y[k] * hf_129[k]
                   + if__179[k];

        t_194[k] = ab_z[k] * hf_129[k]
                   + if__189[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, hf_130, hf_131, hf_132, \
                         hf_133, hf_134, if__130, if__131, if__132, if__133, \
                         if__134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = ab_x[k] * hf_130[k]
                   + if__130[k];

        t_196[k] = ab_x[k] * hf_131[k]
                   + if__131[k];

        t_197[k] = ab_x[k] * hf_132[k]
                   + if__132[k];

        t_198[k] = ab_x[k] * hf_133[k]
                   + if__133[k];

        t_199[k] = ab_x[k] * hf_134[k]
                   + if__134[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, hf_135, hf_136, hf_137, \
                         hf_138, hf_139, if__135, if__136, if__137, if__138, \
                         if__139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = ab_x[k] * hf_135[k]
                   + if__135[k];

        t_201[k] = ab_x[k] * hf_136[k]
                   + if__136[k];

        t_202[k] = ab_x[k] * hf_137[k]
                   + if__137[k];

        t_203[k] = ab_x[k] * hf_138[k]
                   + if__138[k];

        t_204[k] = ab_x[k] * hf_139[k]
                   + if__139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, hf_136, hf_137, \
                         hf_138, hf_139, if__186, if__187, if__188, if__189, \
                         if__199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = ab_y[k] * hf_136[k]
                   + if__186[k];

        t_206[k] = ab_y[k] * hf_137[k]
                   + if__187[k];

        t_207[k] = ab_y[k] * hf_138[k]
                   + if__188[k];

        t_208[k] = ab_y[k] * hf_139[k]
                   + if__189[k];

        t_209[k] = ab_z[k] * hf_139[k]
                   + if__199[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, hf_140, hf_141, hf_142, \
                         hf_143, hf_144, if__140, if__141, if__142, if__143, \
                         if__144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = ab_x[k] * hf_140[k]
                   + if__140[k];

        t_211[k] = ab_x[k] * hf_141[k]
                   + if__141[k];

        t_212[k] = ab_x[k] * hf_142[k]
                   + if__142[k];

        t_213[k] = ab_x[k] * hf_143[k]
                   + if__143[k];

        t_214[k] = ab_x[k] * hf_144[k]
                   + if__144[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, hf_145, hf_146, hf_147, \
                         hf_148, hf_149, if__145, if__146, if__147, if__148, \
                         if__149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = ab_x[k] * hf_145[k]
                   + if__145[k];

        t_216[k] = ab_x[k] * hf_146[k]
                   + if__146[k];

        t_217[k] = ab_x[k] * hf_147[k]
                   + if__147[k];

        t_218[k] = ab_x[k] * hf_148[k]
                   + if__148[k];

        t_219[k] = ab_x[k] * hf_149[k]
                   + if__149[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, hf_146, hf_147, \
                         hf_148, hf_149, if__196, if__197, if__198, if__199, \
                         if__209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = ab_y[k] * hf_146[k]
                   + if__196[k];

        t_221[k] = ab_y[k] * hf_147[k]
                   + if__197[k];

        t_222[k] = ab_y[k] * hf_148[k]
                   + if__198[k];

        t_223[k] = ab_y[k] * hf_149[k]
                   + if__199[k];

        t_224[k] = ab_z[k] * hf_149[k]
                   + if__209[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, hf_150, hf_151, hf_152, \
                         hf_153, hf_154, if__150, if__151, if__152, if__153, \
                         if__154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = ab_x[k] * hf_150[k]
                   + if__150[k];

        t_226[k] = ab_x[k] * hf_151[k]
                   + if__151[k];

        t_227[k] = ab_x[k] * hf_152[k]
                   + if__152[k];

        t_228[k] = ab_x[k] * hf_153[k]
                   + if__153[k];

        t_229[k] = ab_x[k] * hf_154[k]
                   + if__154[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, hf_155, hf_156, hf_157, \
                         hf_158, hf_159, if__155, if__156, if__157, if__158, \
                         if__159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = ab_x[k] * hf_155[k]
                   + if__155[k];

        t_231[k] = ab_x[k] * hf_156[k]
                   + if__156[k];

        t_232[k] = ab_x[k] * hf_157[k]
                   + if__157[k];

        t_233[k] = ab_x[k] * hf_158[k]
                   + if__158[k];

        t_234[k] = ab_x[k] * hf_159[k]
                   + if__159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, ab_z, hf_156, hf_157, \
                         hf_158, hf_159, if__216, if__217, if__218, if__219, \
                         if__229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = ab_y[k] * hf_156[k]
                   + if__216[k];

        t_236[k] = ab_y[k] * hf_157[k]
                   + if__217[k];

        t_237[k] = ab_y[k] * hf_158[k]
                   + if__218[k];

        t_238[k] = ab_y[k] * hf_159[k]
                   + if__219[k];

        t_239[k] = ab_z[k] * hf_159[k]
                   + if__229[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, hf_160, hf_161, hf_162, \
                         hf_163, hf_164, if__160, if__161, if__162, if__163, \
                         if__164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = ab_x[k] * hf_160[k]
                   + if__160[k];

        t_241[k] = ab_x[k] * hf_161[k]
                   + if__161[k];

        t_242[k] = ab_x[k] * hf_162[k]
                   + if__162[k];

        t_243[k] = ab_x[k] * hf_163[k]
                   + if__163[k];

        t_244[k] = ab_x[k] * hf_164[k]
                   + if__164[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, hf_165, hf_166, hf_167, \
                         hf_168, hf_169, if__165, if__166, if__167, if__168, \
                         if__169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = ab_x[k] * hf_165[k]
                   + if__165[k];

        t_246[k] = ab_x[k] * hf_166[k]
                   + if__166[k];

        t_247[k] = ab_x[k] * hf_167[k]
                   + if__167[k];

        t_248[k] = ab_x[k] * hf_168[k]
                   + if__168[k];

        t_249[k] = ab_x[k] * hf_169[k]
                   + if__169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, ab_z, hf_166, hf_167, \
                         hf_168, hf_169, if__226, if__227, if__228, if__229, \
                         if__239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = ab_y[k] * hf_166[k]
                   + if__226[k];

        t_251[k] = ab_y[k] * hf_167[k]
                   + if__227[k];

        t_252[k] = ab_y[k] * hf_168[k]
                   + if__228[k];

        t_253[k] = ab_y[k] * hf_169[k]
                   + if__229[k];

        t_254[k] = ab_z[k] * hf_169[k]
                   + if__239[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, hf_170, hf_171, hf_172, \
                         hf_173, hf_174, if__170, if__171, if__172, if__173, \
                         if__174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = ab_x[k] * hf_170[k]
                   + if__170[k];

        t_256[k] = ab_x[k] * hf_171[k]
                   + if__171[k];

        t_257[k] = ab_x[k] * hf_172[k]
                   + if__172[k];

        t_258[k] = ab_x[k] * hf_173[k]
                   + if__173[k];

        t_259[k] = ab_x[k] * hf_174[k]
                   + if__174[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, hf_175, hf_176, hf_177, \
                         hf_178, hf_179, if__175, if__176, if__177, if__178, \
                         if__179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = ab_x[k] * hf_175[k]
                   + if__175[k];

        t_261[k] = ab_x[k] * hf_176[k]
                   + if__176[k];

        t_262[k] = ab_x[k] * hf_177[k]
                   + if__177[k];

        t_263[k] = ab_x[k] * hf_178[k]
                   + if__178[k];

        t_264[k] = ab_x[k] * hf_179[k]
                   + if__179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, ab_z, hf_176, hf_177, \
                         hf_178, hf_179, if__236, if__237, if__238, if__239, \
                         if__249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = ab_y[k] * hf_176[k]
                   + if__236[k];

        t_266[k] = ab_y[k] * hf_177[k]
                   + if__237[k];

        t_267[k] = ab_y[k] * hf_178[k]
                   + if__238[k];

        t_268[k] = ab_y[k] * hf_179[k]
                   + if__239[k];

        t_269[k] = ab_z[k] * hf_179[k]
                   + if__249[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, hf_180, hf_181, hf_182, \
                         hf_183, hf_184, if__180, if__181, if__182, if__183, \
                         if__184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = ab_x[k] * hf_180[k]
                   + if__180[k];

        t_271[k] = ab_x[k] * hf_181[k]
                   + if__181[k];

        t_272[k] = ab_x[k] * hf_182[k]
                   + if__182[k];

        t_273[k] = ab_x[k] * hf_183[k]
                   + if__183[k];

        t_274[k] = ab_x[k] * hf_184[k]
                   + if__184[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, hf_185, hf_186, hf_187, \
                         hf_188, hf_189, if__185, if__186, if__187, if__188, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = ab_x[k] * hf_185[k]
                   + if__185[k];

        t_276[k] = ab_x[k] * hf_186[k]
                   + if__186[k];

        t_277[k] = ab_x[k] * hf_187[k]
                   + if__187[k];

        t_278[k] = ab_x[k] * hf_188[k]
                   + if__188[k];

        t_279[k] = ab_x[k] * hf_189[k]
                   + if__189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, ab_z, hf_186, hf_187, \
                         hf_188, hf_189, if__246, if__247, if__248, if__249, \
                         if__259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_280[k] = ab_y[k] * hf_186[k]
                   + if__246[k];

        t_281[k] = ab_y[k] * hf_187[k]
                   + if__247[k];

        t_282[k] = ab_y[k] * hf_188[k]
                   + if__248[k];

        t_283[k] = ab_y[k] * hf_189[k]
                   + if__249[k];

        t_284[k] = ab_z[k] * hf_189[k]
                   + if__259[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, hf_190, hf_191, hf_192, \
                         hf_193, hf_194, if__190, if__191, if__192, if__193, \
                         if__194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_285[k] = ab_x[k] * hf_190[k]
                   + if__190[k];

        t_286[k] = ab_x[k] * hf_191[k]
                   + if__191[k];

        t_287[k] = ab_x[k] * hf_192[k]
                   + if__192[k];

        t_288[k] = ab_x[k] * hf_193[k]
                   + if__193[k];

        t_289[k] = ab_x[k] * hf_194[k]
                   + if__194[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, hf_195, hf_196, hf_197, \
                         hf_198, hf_199, if__195, if__196, if__197, if__198, \
                         if__199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_290[k] = ab_x[k] * hf_195[k]
                   + if__195[k];

        t_291[k] = ab_x[k] * hf_196[k]
                   + if__196[k];

        t_292[k] = ab_x[k] * hf_197[k]
                   + if__197[k];

        t_293[k] = ab_x[k] * hf_198[k]
                   + if__198[k];

        t_294[k] = ab_x[k] * hf_199[k]
                   + if__199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, ab_z, hf_196, hf_197, \
                         hf_198, hf_199, if__256, if__257, if__258, if__259, \
                         if__269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_295[k] = ab_y[k] * hf_196[k]
                   + if__256[k];

        t_296[k] = ab_y[k] * hf_197[k]
                   + if__257[k];

        t_297[k] = ab_y[k] * hf_198[k]
                   + if__258[k];

        t_298[k] = ab_y[k] * hf_199[k]
                   + if__259[k];

        t_299[k] = ab_z[k] * hf_199[k]
                   + if__269[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, hf_200, hf_201, hf_202, \
                         hf_203, hf_204, if__200, if__201, if__202, if__203, \
                         if__204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_300[k] = ab_x[k] * hf_200[k]
                   + if__200[k];

        t_301[k] = ab_x[k] * hf_201[k]
                   + if__201[k];

        t_302[k] = ab_x[k] * hf_202[k]
                   + if__202[k];

        t_303[k] = ab_x[k] * hf_203[k]
                   + if__203[k];

        t_304[k] = ab_x[k] * hf_204[k]
                   + if__204[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, hf_205, hf_206, hf_207, \
                         hf_208, hf_209, if__205, if__206, if__207, if__208, \
                         if__209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_305[k] = ab_x[k] * hf_205[k]
                   + if__205[k];

        t_306[k] = ab_x[k] * hf_206[k]
                   + if__206[k];

        t_307[k] = ab_x[k] * hf_207[k]
                   + if__207[k];

        t_308[k] = ab_x[k] * hf_208[k]
                   + if__208[k];

        t_309[k] = ab_x[k] * hf_209[k]
                   + if__209[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, ab_z, hf_206, hf_207, \
                         hf_208, hf_209, if__266, if__267, if__268, if__269, \
                         if__279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_310[k] = ab_y[k] * hf_206[k]
                   + if__266[k];

        t_311[k] = ab_y[k] * hf_207[k]
                   + if__267[k];

        t_312[k] = ab_y[k] * hf_208[k]
                   + if__268[k];

        t_313[k] = ab_y[k] * hf_209[k]
                   + if__269[k];

        t_314[k] = ab_z[k] * hf_209[k]
                   + if__279[k];
    }
}

}  // namespace simdtrf
