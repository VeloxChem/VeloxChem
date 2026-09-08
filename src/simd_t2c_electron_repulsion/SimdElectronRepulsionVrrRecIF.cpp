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


#include "SimdElectronRepulsionVrrRecIF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_if_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
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
    const auto *hd_28 = buffer.data(hd + 28);
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
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);

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

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, hd_1, ip0_0, ip1_0, \
                         id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_0 * hd_1[k]
                 + pb_x[k] * id_2[k];

        t_4[k] = pb_y[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, hd_2, ip0_1, ip0_2, ip1_1, \
                         ip1_2, id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * hd_2[k]
                 + pb_x[k] * id_3[k];

        t_6[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_2[k];

        t_7[k] = pb_z[k] * id_2[k];

        t_8[k] = pb_y[k] * id_3[k];

        t_9[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, hd_0, hd_4, \
                         hf_0, id_4, id_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * hf_0[k];

        t_11[k] = f_3 * hd_0[k]
                  + pb_y[k] * id_4[k];

        t_12[k] = pb_z[k] * id_4[k];

        t_13[k] = f_4 * hd_4[k]
                  + pb_x[k] * id_6[k];

        t_14[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, hd_1, hd_2, hf_2, \
                         hf_3, hf_4, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * hf_2[k];

        t_16[k] = f_5 * hd_1[k]
                  + pa_y[k] * hf_3[k];

        t_17[k] = pb_z[k] * id_6[k];

        t_18[k] = f_3 * hd_2[k]
                  + pb_y[k] * id_7[k];

        t_19[k] = pa_y[k] * hf_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, hd_0, hf_0, hf_1, \
                         id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * hf_0[k];

        t_21[k] = pb_y[k] * id_8[k];

        t_22[k] = f_3 * hd_0[k]
                  + pb_z[k] * id_8[k];

        t_23[k] = pa_z[k] * hf_1[k];

        t_24[k] = pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, hd_1, hd_2, \
                         hd_8, hf_3, hf_4, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * hd_8[k]
                  + pb_x[k] * id_11[k];

        t_26[k] = pa_z[k] * hf_3[k];

        t_27[k] = f_3 * hd_1[k]
                  + pb_z[k] * id_10[k];

        t_28[k] = pb_y[k] * id_11[k];

        t_29[k] = f_5 * hd_2[k]
                  + pa_z[k] * hf_4[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, gf0_0, gf1_0, hd_3, \
                         hd_10, hf_5, id_12, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * gf0_0[k]
                  - f_7 * gf1_0[k]
                  + pa_y[k] * hf_5[k];

        t_31[k] = f_8 * hd_3[k]
                  + pb_y[k] * id_12[k];

        t_32[k] = pb_z[k] * id_12[k];

        t_33[k] = f_9 * hd_10[k]
                  + pb_x[k] * id_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_4, gf1_4, hd_11, hf_16, \
                         id_13, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * id_13[k];

        t_35[k] = f_9 * hd_11[k]
                  + pb_x[k] * id_15[k];

        t_36[k] = f_10 * gf0_4[k]
                  - f_11 * gf1_4[k]
                  + pa_x[k] * hf_16[k];

        t_37[k] = pb_z[k] * id_14[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, hd_5, hf_6, \
                         hf_9, hf_10, ip0_3, ip1_3, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * hd_5[k]
                  + pb_y[k] * id_15[k];

        t_39[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_15[k];

        t_40[k] = pa_y[k] * hf_9[k];

        t_41[k] = pa_z[k] * hf_6[k];

        t_42[k] = pa_y[k] * hf_10[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, hd_4, hd_13, \
                         hf_7, hf_8, hf_11, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * hf_7[k];

        t_44[k] = f_9 * hd_13[k]
                  + pb_x[k] * id_17[k];

        t_45[k] = pa_y[k] * hf_11[k];

        t_46[k] = pa_z[k] * hf_8[k];

        t_47[k] = f_3 * hd_4[k]
                  + pb_z[k] * id_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, gf0_0, gf1_0, hd_8, hf_9, \
                         hf_12, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * hd_8[k]
                  + pb_y[k] * id_18[k];

        t_49[k] = pa_y[k] * hf_12[k];

        t_50[k] = f_6 * gf0_0[k]
                  - f_7 * gf1_0[k]
                  + pa_z[k] * hf_9[k];

        t_51[k] = pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, hd_6, hd_16, hd_17, id_19, \
                         id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * hd_6[k]
                  + pb_z[k] * id_19[k];

        t_53[k] = f_9 * hd_16[k]
                  + pb_x[k] * id_21[k];

        t_54[k] = pb_y[k] * id_20[k];

        t_55[k] = f_9 * hd_17[k]
                  + pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, gf0_6, gf1_6, hd_7, hf_22, \
                         ip0_4, ip1_4, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_21[k];

        t_57[k] = f_8 * hd_7[k]
                  + pb_z[k] * id_21[k];

        t_58[k] = pb_y[k] * id_22[k];

        t_59[k] = f_10 * gf0_6[k]
                  - f_11 * gf1_6[k]
                  + pa_x[k] * hf_22[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, gf0_1, gf1_1, hd_9, \
                         hd_19, hf_13, id_23, id_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * gf0_1[k]
                  - f_13 * gf1_1[k]
                  + pa_y[k] * hf_13[k];

        t_61[k] = f_5 * hd_9[k]
                  + pb_y[k] * id_23[k];

        t_62[k] = pb_z[k] * id_23[k];

        t_63[k] = f_5 * hd_19[k]
                  + pb_x[k] * id_25[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, gf0_7, gf1_7, hd_20, hf_26, \
                         id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * id_24[k];

        t_65[k] = f_5 * hd_20[k]
                  + pb_x[k] * id_26[k];

        t_66[k] = f_12 * gf0_7[k]
                  - f_13 * gf1_7[k]
                  + pa_x[k] * hf_26[k];

        t_67[k] = pb_z[k] * id_25[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, hd_9, hd_11, hf_13, \
                         hf_14, ip0_5, ip1_5, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * hd_11[k]
                  + pb_y[k] * id_26[k];

        t_69[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_26[k];

        t_70[k] = pa_z[k] * hf_13[k];

        t_71[k] = pa_z[k] * hf_14[k];

        t_72[k] = f_3 * hd_9[k]
                  + pb_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, hd_10, hd_23, hd_24, \
                         hf_15, hf_16, id_28, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * hf_15[k];

        t_74[k] = f_5 * hd_23[k]
                  + pb_x[k] * id_29[k];

        t_75[k] = f_5 * hd_24[k]
                  + pb_x[k] * id_30[k];

        t_76[k] = pa_z[k] * hf_16[k];

        t_77[k] = f_3 * hd_10[k]
                  + pb_z[k] * id_28[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, hd_11, hd_14, hd_15, \
                         hf_17, hf_18, hf_19, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * hd_14[k]
                  + pb_y[k] * id_30[k];

        t_79[k] = f_5 * hd_11[k]
                  + pa_z[k] * hf_17[k];

        t_80[k] = pa_y[k] * hf_18[k];

        t_81[k] = f_3 * hd_15[k]
                  + pb_y[k] * id_31[k];

        t_82[k] = pa_y[k] * hf_19[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, hd_12, hd_16, hd_26, \
                         hd_27, hf_20, hf_21, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_5 * hd_26[k]
                  + pb_x[k] * id_32[k];

        t_84[k] = f_5 * hd_27[k]
                  + pb_x[k] * id_33[k];

        t_85[k] = pa_y[k] * hf_20[k];

        t_86[k] = f_5 * hd_16[k]
                  + pa_y[k] * hf_21[k];

        t_87[k] = f_8 * hd_12[k]
                  + pb_z[k] * id_32[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, gf0_2, gf1_2, hd_17, hf_18, \
                         hf_22, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * hd_17[k]
                  + pb_y[k] * id_34[k];

        t_89[k] = pa_y[k] * hf_22[k];

        t_90[k] = f_12 * gf0_2[k]
                  - f_13 * gf1_2[k]
                  + pa_z[k] * hf_18[k];

        t_91[k] = pb_y[k] * id_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, hd_15, hd_30, hd_31, id_35, \
                         id_36, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * hd_15[k]
                  + pb_z[k] * id_35[k];

        t_93[k] = f_5 * hd_30[k]
                  + pb_x[k] * id_37[k];

        t_94[k] = pb_y[k] * id_36[k];

        t_95[k] = f_5 * hd_31[k]
                  + pb_x[k] * id_38[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, gf0_8, gf1_8, hd_16, hf_33, \
                         ip0_6, ip1_6, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_37[k];

        t_97[k] = f_5 * hd_16[k]
                  + pb_z[k] * id_37[k];

        t_98[k] = pb_y[k] * id_38[k];

        t_99[k] = f_12 * gf0_8[k]
                  - f_13 * gf1_8[k]
                  + pa_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, pb_z, gf0_3, gf1_3, \
                         hd_18, hd_33, hf_23, id_39, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_10 * gf0_3[k]
                   - f_11 * gf1_3[k]
                   + pa_y[k] * hf_23[k];

        t_101[k] = f_9 * hd_18[k]
                   + pb_y[k] * id_39[k];

        t_102[k] = pb_z[k] * id_39[k];

        t_103[k] = f_8 * hd_33[k]
                   + pb_x[k] * id_41[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, pb_z, gf0_9, gf1_9, hd_34, \
                         hf_37, id_40, id_41, id_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * id_40[k];

        t_105[k] = f_8 * hd_34[k]
                   + pb_x[k] * id_42[k];

        t_106[k] = f_6 * gf0_9[k]
                   - f_7 * gf1_9[k]
                   + pa_x[k] * hf_37[k];

        t_107[k] = pb_z[k] * id_41[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, hd_18, hd_20, \
                         hf_23, hf_24, ip0_7, ip1_7, id_42, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_9 * hd_20[k]
                   + pb_y[k] * id_42[k];

        t_109[k] = f_1 * ip0_7[k]
                   - f_2 * ip1_7[k]
                   + pb_z[k] * id_42[k];

        t_110[k] = pa_z[k] * hf_23[k];

        t_111[k] = pa_z[k] * hf_24[k];

        t_112[k] = f_3 * hd_18[k]
                   + pb_z[k] * id_43[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, hd_19, hd_36, \
                         hd_37, hf_25, hf_26, id_44, id_45, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * hf_25[k];

        t_114[k] = f_8 * hd_36[k]
                   + pb_x[k] * id_45[k];

        t_115[k] = f_8 * hd_37[k]
                   + pb_x[k] * id_46[k];

        t_116[k] = pa_z[k] * hf_26[k];

        t_117[k] = f_3 * hd_19[k]
                   + pb_z[k] * id_44[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pa_z, pb_y, gf0_5, gf1_5, hd_20, \
                         hd_24, hd_25, hf_27, hf_28, id_46, id_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * hd_24[k]
                   + pb_y[k] * id_46[k];

        t_119[k] = f_5 * hd_20[k]
                   + pa_z[k] * hf_27[k];

        t_120[k] = f_6 * gf0_5[k]
                   - f_7 * gf1_5[k]
                   + pa_y[k] * hf_28[k];

        t_121[k] = f_8 * hd_25[k]
                   + pb_y[k] * id_47[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pb_x, pb_z, hd_21, hd_39, hd_40, hd_41, \
                         id_47, id_48, id_49, id_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * hd_21[k]
                   + pb_z[k] * id_47[k];

        t_123[k] = f_8 * hd_39[k]
                   + pb_x[k] * id_48[k];

        t_124[k] = f_8 * hd_40[k]
                   + pb_x[k] * id_49[k];

        t_125[k] = f_8 * hd_41[k]
                   + pb_x[k] * id_50[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, gf0_11, gf1_11, hd_22, hd_28, \
                         hf_38, id_48, id_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_6 * gf0_11[k]
                   - f_7 * gf1_11[k]
                   + pa_x[k] * hf_38[k];

        t_127[k] = f_8 * hd_22[k]
                   + pb_z[k] * id_48[k];

        t_128[k] = f_8 * hd_28[k]
                   + pb_y[k] * id_50[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, gf0_12, gf1_12, hd_29, \
                         hf_29, hf_30, hf_39, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_6 * gf0_12[k]
                   - f_7 * gf1_12[k]
                   + pa_x[k] * hf_39[k];

        t_130[k] = pa_y[k] * hf_29[k];

        t_131[k] = f_3 * hd_29[k]
                   + pb_y[k] * id_51[k];

        t_132[k] = pa_y[k] * hf_30[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pa_y, pb_x, pb_z, hd_26, hd_30, \
                         hd_43, hd_44, hf_31, hf_32, id_52, id_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_8 * hd_43[k]
                   + pb_x[k] * id_52[k];

        t_134[k] = f_8 * hd_44[k]
                   + pb_x[k] * id_53[k];

        t_135[k] = pa_y[k] * hf_31[k];

        t_136[k] = f_5 * hd_30[k]
                   + pa_y[k] * hf_32[k];

        t_137[k] = f_5 * hd_26[k]
                   + pb_z[k] * id_52[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pa_z, pb_y, gf0_5, gf1_5, hd_31, \
                         hf_29, hf_33, id_54, id_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * hd_31[k]
                   + pb_y[k] * id_54[k];

        t_139[k] = pa_y[k] * hf_33[k];

        t_140[k] = f_10 * gf0_5[k]
                   - f_11 * gf1_5[k]
                   + pa_z[k] * hf_29[k];

        t_141[k] = pb_y[k] * id_55[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, hd_29, hd_46, hd_47, \
                         id_55, id_56, id_57, id_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_9 * hd_29[k]
                   + pb_z[k] * id_55[k];

        t_143[k] = f_8 * hd_46[k]
                   + pb_x[k] * id_57[k];

        t_144[k] = pb_y[k] * id_56[k];

        t_145[k] = f_8 * hd_47[k]
                   + pb_x[k] * id_58[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, pb_z, gf0_14, gf1_14, hd_30, \
                         hf_43, ip0_8, ip1_8, id_57, id_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * ip0_8[k]
                   - f_2 * ip1_8[k]
                   + pb_y[k] * id_57[k];

        t_147[k] = f_9 * hd_30[k]
                   + pb_z[k] * id_57[k];

        t_148[k] = pb_y[k] * id_58[k];

        t_149[k] = f_6 * gf0_14[k]
                   - f_7 * gf1_14[k]
                   + pa_x[k] * hf_43[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_x, pb_x, pb_y, pb_z, hd_32, \
                         hd_48, hd_49, hf_44, id_59, id_60, id_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * hd_48[k]
                   + pa_x[k] * hf_44[k];

        t_151[k] = f_4 * hd_32[k]
                   + pb_y[k] * id_59[k];

        t_152[k] = pb_z[k] * id_59[k];

        t_153[k] = f_3 * hd_49[k]
                   + pb_x[k] * id_61[k];

        t_154[k] = pb_z[k] * id_60[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_x, pb_x, pb_z, hd_50, hf_46, \
                         hf_47, hf_48, id_61, id_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_3 * hd_50[k]
                   + pb_x[k] * id_62[k];

        t_156[k] = pa_x[k] * hf_46[k];

        t_157[k] = pb_z[k] * id_61[k];

        t_158[k] = pa_x[k] * hf_47[k];

        t_159[k] = pa_x[k] * hf_48[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_z, pb_x, pb_z, hd_32, hd_53, \
                         hf_34, hf_35, hf_36, id_63, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_z[k] * hf_34[k];

        t_161[k] = pa_z[k] * hf_35[k];

        t_162[k] = f_3 * hd_32[k]
                   + pb_z[k] * id_63[k];

        t_163[k] = pa_z[k] * hf_36[k];

        t_164[k] = f_3 * hd_53[k]
                   + pb_x[k] * id_64[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, pa_x, pb_x, hd_54, hd_55, \
                         hf_49, hf_50, hf_51, hf_52, hf_53, id_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_3 * hd_54[k]
                   + pb_x[k] * id_65[k];

        t_166[k] = pa_x[k] * hf_49[k];

        t_167[k] = pa_x[k] * hf_50[k];

        t_168[k] = pa_x[k] * hf_51[k];

        t_169[k] = pa_x[k] * hf_52[k];

        t_170[k] = f_5 * hd_55[k]
                   + pa_x[k] * hf_53[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pb_x, pb_y, pb_z, hd_35, hd_38, hd_56, \
                         hd_57, id_66, id_67, id_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_5 * hd_38[k]
                   + pb_y[k] * id_66[k];

        t_172[k] = f_8 * hd_35[k]
                   + pb_z[k] * id_66[k];

        t_173[k] = f_3 * hd_56[k]
                   + pb_x[k] * id_67[k];

        t_174[k] = f_3 * hd_57[k]
                   + pb_x[k] * id_68[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, pa_x, pb_x, hd_58, hd_59, \
                         hf_54, hf_55, hf_56, hf_57, hf_58, id_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * hd_58[k]
                   + pb_x[k] * id_69[k];

        t_176[k] = pa_x[k] * hf_54[k];

        t_177[k] = pa_x[k] * hf_55[k];

        t_178[k] = pa_x[k] * hf_56[k];

        t_179[k] = pa_x[k] * hf_57[k];

        t_180[k] = f_5 * hd_59[k]
                   + pa_x[k] * hf_58[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_x, pb_y, pb_z, hd_38, hd_42, hd_60, \
                         hd_61, id_70, id_71, id_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * hd_42[k]
                   + pb_y[k] * id_70[k];

        t_182[k] = f_5 * hd_38[k]
                   + pb_z[k] * id_70[k];

        t_183[k] = f_3 * hd_60[k]
                   + pb_x[k] * id_71[k];

        t_184[k] = f_3 * hd_61[k]
                   + pb_x[k] * id_72[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, pa_x, pa_y, pb_x, hd_62, \
                         hf_40, hf_59, hf_60, hf_61, hf_62, id_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_3 * hd_62[k]
                   + pb_x[k] * id_73[k];

        t_186[k] = pa_x[k] * hf_59[k];

        t_187[k] = pa_x[k] * hf_60[k];

        t_188[k] = pa_x[k] * hf_61[k];

        t_189[k] = pa_x[k] * hf_62[k];

        t_190[k] = pa_y[k] * hf_40[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pa_y, pb_x, pb_y, hd_45, hd_64, \
                         hd_65, hf_41, hf_42, id_74, id_75, id_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_3 * hd_45[k]
                   + pb_y[k] * id_74[k];

        t_192[k] = pa_y[k] * hf_41[k];

        t_193[k] = f_3 * hd_64[k]
                   + pb_x[k] * id_75[k];

        t_194[k] = f_3 * hd_65[k]
                   + pb_x[k] * id_76[k];

        t_195[k] = pa_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, t_201, pa_x, pb_y, hd_67, hf_63, \
                         hf_64, hf_65, hf_66, hf_67, id_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * hf_63[k];

        t_197[k] = pa_x[k] * hf_64[k];

        t_198[k] = pa_x[k] * hf_65[k];

        t_199[k] = pa_x[k] * hf_66[k];

        t_200[k] = f_5 * hd_67[k]
                   + pa_x[k] * hf_67[k];

        t_201[k] = pb_y[k] * id_77[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, hd_45, hd_68, hd_69, \
                         id_77, id_78, id_79, id_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_4 * hd_45[k]
                   + pb_z[k] * id_77[k];

        t_203[k] = f_3 * hd_68[k]
                   + pb_x[k] * id_79[k];

        t_204[k] = pb_y[k] * id_78[k];

        t_205[k] = f_3 * hd_69[k]
                   + pb_x[k] * id_80[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, pa_x, pb_x, pb_y, hf_69, hf_70, \
                         hf_71, ip0_9, ip1_9, id_80, id_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_x[k] * hf_69[k];

        t_207[k] = pa_x[k] * hf_70[k];

        t_208[k] = pb_y[k] * id_80[k];

        t_209[k] = pa_x[k] * hf_71[k];

        t_210[k] = f_1 * ip0_9[k]
                   - f_2 * ip1_9[k]
                   + pb_x[k] * id_81[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pb_x, pb_y, pb_z, hd_48, id_81, \
                         id_82, id_83, id_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_0 * hd_48[k]
                   + pb_y[k] * id_81[k];

        t_212[k] = pb_z[k] * id_81[k];

        t_213[k] = pb_x[k] * id_82[k];

        t_214[k] = pb_x[k] * id_83[k];

        t_215[k] = pb_x[k] * id_84[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_y, pb_z, hd_49, hd_50, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_82, id_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_0 * hd_49[k]
                   + f_1 * ip0_10[k]
                   - f_2 * ip1_10[k]
                   + pb_y[k] * id_82[k];

        t_217[k] = pb_z[k] * id_82[k];

        t_218[k] = f_0 * hd_50[k]
                   + pb_y[k] * id_84[k];

        t_219[k] = f_1 * ip0_11[k]
                   - f_2 * ip1_11[k]
                   + pb_z[k] * id_84[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, t_225, pa_z, pb_x, pb_z, hd_48, \
                         hf_44, hf_45, id_85, id_86, id_87, id_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_z[k] * hf_44[k];

        t_221[k] = pa_z[k] * hf_45[k];

        t_222[k] = f_3 * hd_48[k]
                   + pb_z[k] * id_85[k];

        t_223[k] = pb_x[k] * id_86[k];

        t_224[k] = pb_x[k] * id_87[k];

        t_225[k] = pb_x[k] * id_88[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_z, pb_y, pb_z, hd_49, hd_50, hd_54, \
                         hf_46, hf_48, id_86, id_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * hf_46[k];

        t_227[k] = f_3 * hd_49[k]
                   + pb_z[k] * id_86[k];

        t_228[k] = f_4 * hd_54[k]
                   + pb_y[k] * id_88[k];

        t_229[k] = f_5 * hd_50[k]
                   + pa_z[k] * hf_48[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, pb_x, pb_y, pb_z, hd_51, hd_55, \
                         ip0_12, ip1_12, id_89, id_90, id_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * ip0_12[k]
                   - f_2 * ip1_12[k]
                   + pb_x[k] * id_89[k];

        t_231[k] = f_9 * hd_55[k]
                   + pb_y[k] * id_89[k];

        t_232[k] = f_8 * hd_51[k]
                   + pb_z[k] * id_89[k];

        t_233[k] = pb_x[k] * id_90[k];

        t_234[k] = pb_x[k] * id_91[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pa_z, pb_x, pb_y, pb_z, gf0_9, gf1_9, \
                         hd_52, hd_58, hf_49, id_90, id_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = pb_x[k] * id_92[k];

        t_236[k] = f_6 * gf0_9[k]
                   - f_7 * gf1_9[k]
                   + pa_z[k] * hf_49[k];

        t_237[k] = f_8 * hd_52[k]
                   + pb_z[k] * id_90[k];

        t_238[k] = f_9 * hd_58[k]
                   + pb_y[k] * id_92[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_y, pb_x, pb_y, pb_z, gf0_12, gf1_12, \
                         hd_55, hd_59, hf_57, ip0_13, ip1_13, id_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * gf0_12[k]
                   - f_11 * gf1_12[k]
                   + pa_y[k] * hf_57[k];

        t_240[k] = f_1 * ip0_13[k]
                   - f_2 * ip1_13[k]
                   + pb_x[k] * id_93[k];

        t_241[k] = f_5 * hd_59[k]
                   + pb_y[k] * id_93[k];

        t_242[k] = f_5 * hd_55[k]
                   + pb_z[k] * id_93[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, pa_z, pb_x, pb_z, gf0_10, gf1_10, \
                         hd_56, hf_54, id_94, id_95, id_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pb_x[k] * id_94[k];

        t_244[k] = pb_x[k] * id_95[k];

        t_245[k] = pb_x[k] * id_96[k];

        t_246[k] = f_12 * gf0_10[k]
                   - f_13 * gf1_10[k]
                   + pa_z[k] * hf_54[k];

        t_247[k] = f_5 * hd_56[k]
                   + pb_z[k] * id_94[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_y, pb_x, pb_y, gf0_13, gf1_13, hd_62, \
                         hd_63, hf_62, ip0_14, ip1_14, id_96, id_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_5 * hd_62[k]
                   + pb_y[k] * id_96[k];

        t_249[k] = f_12 * gf0_13[k]
                   - f_13 * gf1_13[k]
                   + pa_y[k] * hf_62[k];

        t_250[k] = f_1 * ip0_14[k]
                   - f_2 * ip1_14[k]
                   + pb_x[k] * id_97[k];

        t_251[k] = f_8 * hd_63[k]
                   + pb_y[k] * id_97[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, pa_z, pb_x, pb_z, gf0_11, gf1_11, \
                         hd_59, hf_59, id_97, id_98, id_99, id_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_9 * hd_59[k]
                   + pb_z[k] * id_97[k];

        t_253[k] = pb_x[k] * id_98[k];

        t_254[k] = pb_x[k] * id_99[k];

        t_255[k] = pb_x[k] * id_100[k];

        t_256[k] = f_10 * gf0_11[k]
                   - f_11 * gf1_11[k]
                   + pa_z[k] * hf_59[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_y, pb_y, pb_z, gf0_14, gf1_14, hd_60, \
                         hd_66, hf_66, hf_67, id_98, id_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_9 * hd_60[k]
                   + pb_z[k] * id_98[k];

        t_258[k] = f_8 * hd_66[k]
                   + pb_y[k] * id_100[k];

        t_259[k] = f_6 * gf0_14[k]
                   - f_7 * gf1_14[k]
                   + pa_y[k] * hf_66[k];

        t_260[k] = pa_y[k] * hf_67[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, pa_y, pb_x, pb_y, hd_67, hf_68, \
                         id_101, id_102, id_103, id_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * hd_67[k]
                   + pb_y[k] * id_101[k];

        t_262[k] = pa_y[k] * hf_68[k];

        t_263[k] = pb_x[k] * id_102[k];

        t_264[k] = pb_x[k] * id_103[k];

        t_265[k] = pb_x[k] * id_104[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_y, pb_y, pb_z, hd_64, hd_68, hd_69, \
                         hf_69, hf_71, id_102, id_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_5 * hd_68[k]
                   + pa_y[k] * hf_69[k];

        t_267[k] = f_4 * hd_64[k]
                   + pb_z[k] * id_102[k];

        t_268[k] = f_3 * hd_69[k]
                   + pb_y[k] * id_104[k];

        t_269[k] = pa_y[k] * hf_71[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, hd_67, \
                         ip0_15, ip1_15, id_105, id_106, id_107, \
                         id_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * ip0_15[k]
                   - f_2 * ip1_15[k]
                   + pb_x[k] * id_105[k];

        t_271[k] = pb_y[k] * id_105[k];

        t_272[k] = f_0 * hd_67[k]
                   + pb_z[k] * id_105[k];

        t_273[k] = pb_x[k] * id_106[k];

        t_274[k] = pb_x[k] * id_107[k];

        t_275[k] = pb_x[k] * id_108[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_y, pb_z, hd_68, hd_69, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_106, id_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_106[k];

        t_277[k] = f_0 * hd_68[k]
                   + pb_z[k] * id_106[k];

        t_278[k] = pb_y[k] * id_108[k];

        t_279[k] = f_0 * hd_69[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_108[k];
    }
}

auto
compute_prim_if_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_42 = buffer.data(gf1 + 42);
    const auto *gf1_46 = buffer.data(gf1 + 46);
    const auto *gf1_54 = buffer.data(gf1 + 54);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
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
    const auto *hd_28 = buffer.data(hd + 28);
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
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
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
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_81 = buffer.data(hf + 81);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, hd_1, hd_2, ip0_0, \
                         ip1_0, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_0 * hd_1[k]
                 + pb_x[k] * id_1[k];

        t_4[k] = f_0 * hd_2[k]
                 + pb_x[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, hf_0, ip0_1, ip0_2, ip1_1, \
                         ip1_2, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_6[k] = pb_y[k] * id_2[k];

        t_7[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_8[k] = pa_y[k] * hf_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, hd_0, hd_1, hd_2, hd_4, \
                         hf_3, id_3, id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * hd_0[k]
                 + pb_y[k] * id_3[k];

        t_10[k] = f_4 * hd_4[k]
                  + pb_x[k] * id_4[k];

        t_11[k] = f_5 * hd_1[k]
                  + pa_y[k] * hf_3[k];

        t_12[k] = f_3 * hd_2[k]
                  + pb_y[k] * id_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_z, hd_0, hd_8, \
                         hf_0, hf_3, hf_5, id_6, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hf_5[k];

        t_14[k] = pa_z[k] * hf_0[k];

        t_15[k] = f_3 * hd_0[k]
                  + pb_z[k] * id_6[k];

        t_16[k] = f_4 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_17[k] = pa_z[k] * hf_3[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_z, gf0_0, gf1_0, hd_1, hd_2, hf_5, \
                         hf_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * hd_1[k]
                  + pb_z[k] * id_7[k];

        t_19[k] = f_5 * hd_2[k]
                  + pa_z[k] * hf_5[k];

        t_20[k] = f_6 * gf0_0[k]
                  - f_7 * gf1_0[k]
                  + pa_y[k] * hf_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pb_x, pb_y, pb_z, gf0_5, gf1_13, \
                         hd_3, hd_10, hf_14, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * hd_3[k]
                  + pb_y[k] * id_9[k];

        t_22[k] = pb_z[k] * id_9[k];

        t_23[k] = f_9 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_24[k] = f_10 * gf0_5[k]
                  - f_11 * gf1_13[k]
                  + pa_x[k] * hf_14[k];

        t_25[k] = pb_z[k] * id_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, pb_z, hd_5, hf_7, hf_9, \
                         ip0_3, ip1_3, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hd_5[k]
                  + pb_y[k] * id_11[k];

        t_27[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_11[k];

        t_28[k] = pa_y[k] * hf_9[k];

        t_29[k] = pa_z[k] * hf_7[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_y, pb_z, gf0_0, gf1_0, hd_4, \
                         hd_8, hf_8, hf_10, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hd_4[k]
                  + pb_z[k] * id_12[k];

        t_31[k] = f_3 * hd_8[k]
                  + pb_y[k] * id_13[k];

        t_32[k] = pa_y[k] * hf_10[k];

        t_33[k] = f_6 * gf0_0[k]
                  - f_7 * gf1_0[k]
                  + pa_z[k] * hf_8[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, hd_6, hd_7, hd_16, \
                         ip0_4, ip1_4, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * id_14[k];

        t_35[k] = f_8 * hd_6[k]
                  + pb_z[k] * id_14[k];

        t_36[k] = f_9 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_37[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_15[k];

        t_38[k] = f_8 * hd_7[k]
                  + pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pa_y, pb_y, gf0_1, gf0_8, gf1_6, \
                         gf1_19, hd_9, hf_11, hf_23, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_y[k] * id_16[k];

        t_40[k] = f_10 * gf0_8[k]
                  - f_11 * gf1_19[k]
                  + pa_x[k] * hf_23[k];

        t_41[k] = f_12 * gf0_1[k]
                  - f_13 * gf1_6[k]
                  + pa_y[k] * hf_11[k];

        t_42[k] = f_5 * hd_9[k]
                  + pb_y[k] * id_17[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_z, gf0_9, gf1_22, hd_18, \
                         hf_27, id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_z[k] * id_17[k];

        t_44[k] = f_5 * hd_18[k]
                  + pb_x[k] * id_18[k];

        t_45[k] = f_12 * gf0_9[k]
                  - f_13 * gf1_22[k]
                  + pa_x[k] * hf_27[k];

        t_46[k] = pb_z[k] * id_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, hd_9, hd_11, hf_11, \
                         hf_14, ip0_5, ip1_5, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * hd_11[k]
                  + pb_y[k] * id_19[k];

        t_48[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_19[k];

        t_49[k] = pa_z[k] * hf_11[k];

        t_50[k] = f_3 * hd_9[k]
                  + pb_z[k] * id_20[k];

        t_51[k] = pa_z[k] * hf_14[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pa_z, pb_y, pb_z, hd_10, hd_11, hd_13, \
                         hf_16, hf_17, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * hd_10[k]
                  + pb_z[k] * id_21[k];

        t_53[k] = f_8 * hd_13[k]
                  + pb_y[k] * id_22[k];

        t_54[k] = f_5 * hd_11[k]
                  + pa_z[k] * hf_16[k];

        t_55[k] = pa_y[k] * hf_17[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_y, pb_z, hd_12, hd_15, hd_16, \
                         hf_19, hf_21, hf_23, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * hf_19[k];

        t_57[k] = f_5 * hd_15[k]
                  + pa_y[k] * hf_21[k];

        t_58[k] = f_8 * hd_12[k]
                  + pb_z[k] * id_23[k];

        t_59[k] = f_3 * hd_16[k]
                  + pb_y[k] * id_24[k];

        t_60[k] = pa_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_z, pb_x, pb_y, pb_z, gf0_2, gf1_8, hd_14, \
                         hd_27, hf_17, id_25, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_12 * gf0_2[k]
                  - f_13 * gf1_8[k]
                  + pa_z[k] * hf_17[k];

        t_62[k] = pb_y[k] * id_25[k];

        t_63[k] = f_5 * hd_14[k]
                  + pb_z[k] * id_25[k];

        t_64[k] = f_5 * hd_27[k]
                  + pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_y, pb_z, gf0_10, gf1_26, hd_15, \
                         hf_37, ip0_6, ip1_6, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_26[k];

        t_66[k] = f_5 * hd_15[k]
                  + pb_z[k] * id_26[k];

        t_67[k] = pb_y[k] * id_27[k];

        t_68[k] = f_12 * gf0_10[k]
                  - f_13 * gf1_26[k]
                  + pa_x[k] * hf_37[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_x, pb_y, pb_z, gf0_3, gf1_11, hd_17, \
                         hd_29, hf_24, id_28, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * gf0_3[k]
                  - f_11 * gf1_11[k]
                  + pa_y[k] * hf_24[k];

        t_70[k] = f_9 * hd_17[k]
                  + pb_y[k] * id_28[k];

        t_71[k] = pb_z[k] * id_28[k];

        t_72[k] = f_8 * hd_29[k]
                  + pb_x[k] * id_29[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_y, pb_z, gf0_11, gf1_30, hd_19, \
                         hf_39, ip0_7, ip1_7, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * gf0_11[k]
                  - f_7 * gf1_30[k]
                  + pa_x[k] * hf_39[k];

        t_74[k] = pb_z[k] * id_29[k];

        t_75[k] = f_9 * hd_19[k]
                  + pb_y[k] * id_30[k];

        t_76[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_30[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_z, pb_y, pb_z, hd_17, hd_18, hd_22, \
                         hf_24, hf_27, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * hf_24[k];

        t_78[k] = f_3 * hd_17[k]
                  + pb_z[k] * id_31[k];

        t_79[k] = pa_z[k] * hf_27[k];

        t_80[k] = f_3 * hd_18[k]
                  + pb_z[k] * id_32[k];

        t_81[k] = f_5 * hd_22[k]
                  + pb_y[k] * id_33[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pa_z, pb_z, gf0_6, gf1_15, hd_19, hd_20, \
                         hf_29, hf_30, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_5 * hd_19[k]
                  + pa_z[k] * hf_29[k];

        t_83[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_15[k]
                  + pa_y[k] * hf_30[k];

        t_84[k] = f_8 * hd_20[k]
                  + pb_z[k] * id_34[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pa_x, pb_y, pb_z, gf0_13, gf1_39, hd_21, hd_24, \
                         hf_40, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_39[k]
                  + pa_x[k] * hf_40[k];

        t_86[k] = f_8 * hd_21[k]
                  + pb_z[k] * id_35[k];

        t_87[k] = f_8 * hd_24[k]
                  + pb_y[k] * id_36[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, gf0_15, gf1_42, hd_26, hf_31, \
                         hf_33, hf_35, hf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_42[k]
                  + pa_x[k] * hf_41[k];

        t_89[k] = pa_y[k] * hf_31[k];

        t_90[k] = pa_y[k] * hf_33[k];

        t_91[k] = f_5 * hd_26[k]
                  + pa_y[k] * hf_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_y, pb_z, gf0_6, gf1_15, hd_23, \
                         hd_27, hf_31, hf_37, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * hd_23[k]
                  + pb_z[k] * id_37[k];

        t_93[k] = f_3 * hd_27[k]
                  + pb_y[k] * id_38[k];

        t_94[k] = pa_y[k] * hf_37[k];

        t_95[k] = f_10 * gf0_6[k]
                  - f_11 * gf1_15[k]
                  + pa_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, hd_25, hd_26, hd_33, \
                         ip0_8, ip1_8, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_y[k] * id_39[k];

        t_97[k] = f_9 * hd_25[k]
                  + pb_z[k] * id_39[k];

        t_98[k] = f_8 * hd_33[k]
                  + pb_x[k] * id_41[k];

        t_99[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_40[k];

        t_100[k] = f_9 * hd_26[k]
                   + pb_z[k] * id_40[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, gf0_17, gf1_54, hd_28, hd_34, \
                         hf_44, hf_45, id_41, id_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_y[k] * id_41[k];

        t_102[k] = f_6 * gf0_17[k]
                   - f_7 * gf1_54[k]
                   + pa_x[k] * hf_44[k];

        t_103[k] = f_5 * hd_34[k]
                   + pa_x[k] * hf_45[k];

        t_104[k] = f_4 * hd_28[k]
                   + pb_y[k] * id_42[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_x, pa_z, pb_x, hd_35, hf_38, \
                         hf_48, hf_50, hf_51, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * hd_35[k]
                   + pb_x[k] * id_43[k];

        t_106[k] = pa_x[k] * hf_48[k];

        t_107[k] = pa_x[k] * hf_50[k];

        t_108[k] = pa_x[k] * hf_51[k];

        t_109[k] = pa_z[k] * hf_38[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_x, pb_z, hd_28, hd_40, hf_53, \
                         hf_54, hf_55, hf_56, id_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * hd_28[k]
                   + pb_z[k] * id_44[k];

        t_111[k] = pa_x[k] * hf_53[k];

        t_112[k] = pa_x[k] * hf_54[k];

        t_113[k] = pa_x[k] * hf_55[k];

        t_114[k] = f_5 * hd_40[k]
                   + pa_x[k] * hf_56[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, pa_x, pb_z, hd_30, hd_43, \
                         hf_59, hf_60, hf_61, hf_62, hf_63, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_8 * hd_30[k]
                   + pb_z[k] * id_45[k];

        t_116[k] = pa_x[k] * hf_59[k];

        t_117[k] = pa_x[k] * hf_60[k];

        t_118[k] = pa_x[k] * hf_61[k];

        t_119[k] = pa_x[k] * hf_62[k];

        t_120[k] = f_5 * hd_43[k]
                   + pa_x[k] * hf_63[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, t_126, pa_x, pa_y, pb_z, hd_31, \
                         hf_42, hf_66, hf_67, hf_68, hf_69, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * hd_31[k]
                   + pb_z[k] * id_46[k];

        t_122[k] = pa_x[k] * hf_66[k];

        t_123[k] = pa_x[k] * hf_67[k];

        t_124[k] = pa_x[k] * hf_68[k];

        t_125[k] = pa_x[k] * hf_69[k];

        t_126[k] = pa_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_x, pa_y, hd_48, hf_43, hf_70, \
                         hf_71, hf_72, hf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_y[k] * hf_43[k];

        t_128[k] = pa_x[k] * hf_70[k];

        t_129[k] = pa_x[k] * hf_71[k];

        t_130[k] = pa_x[k] * hf_72[k];

        t_131[k] = f_5 * hd_48[k]
                   + pa_x[k] * hf_74[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pa_x, pb_x, pb_z, hd_32, hd_50, \
                         hf_78, hf_79, hf_81, id_47, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * hd_32[k]
                   + pb_z[k] * id_47[k];

        t_133[k] = f_3 * hd_50[k]
                   + pb_x[k] * id_48[k];

        t_134[k] = pa_x[k] * hf_78[k];

        t_135[k] = pa_x[k] * hf_79[k];

        t_136[k] = pa_x[k] * hf_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pb_x, pb_y, hd_34, hd_35, ip0_9, \
                         ip0_10, ip1_9, ip1_10, id_49, id_50, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_1 * ip0_9[k]
                   - f_2 * ip1_9[k]
                   + pb_x[k] * id_49[k];

        t_138[k] = f_0 * hd_34[k]
                   + pb_y[k] * id_49[k];

        t_139[k] = pb_x[k] * id_50[k];

        t_140[k] = pb_x[k] * id_51[k];

        t_141[k] = f_0 * hd_35[k]
                   + f_1 * ip0_10[k]
                   - f_2 * ip1_10[k]
                   + pb_y[k] * id_50[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, pa_z, pb_y, pb_z, hd_34, hd_36, \
                         hf_45, ip0_11, ip1_11, id_50, id_51, id_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = pb_z[k] * id_50[k];

        t_143[k] = f_0 * hd_36[k]
                   + pb_y[k] * id_51[k];

        t_144[k] = f_1 * ip0_11[k]
                   - f_2 * ip1_11[k]
                   + pb_z[k] * id_51[k];

        t_145[k] = pa_z[k] * hf_45[k];

        t_146[k] = f_3 * hd_34[k]
                   + pb_z[k] * id_52[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pb_y, pb_z, hd_35, hd_36, hd_39, \
                         hf_48, hf_51, id_53, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pa_z[k] * hf_48[k];

        t_148[k] = f_3 * hd_35[k]
                   + pb_z[k] * id_53[k];

        t_149[k] = f_4 * hd_39[k]
                   + pb_y[k] * id_54[k];

        t_150[k] = f_5 * hd_36[k]
                   + pa_z[k] * hf_51[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pb_x, pb_z, hd_37, ip0_12, ip1_12, id_55, \
                         id_56, id_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_1 * ip0_12[k]
                   - f_2 * ip1_12[k]
                   + pb_x[k] * id_55[k];

        t_152[k] = f_8 * hd_37[k]
                   + pb_z[k] * id_55[k];

        t_153[k] = pb_x[k] * id_56[k];

        t_154[k] = pb_x[k] * id_57[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pa_z, pb_y, pb_z, gf0_11, gf1_30, hd_38, hd_42, \
                         hf_52, id_56, id_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_6 * gf0_11[k]
                   - f_7 * gf1_30[k]
                   + pa_z[k] * hf_52[k];

        t_156[k] = f_8 * hd_38[k]
                   + pb_z[k] * id_56[k];

        t_157[k] = f_9 * hd_42[k]
                   + pb_y[k] * id_57[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_y, pb_x, pb_z, gf0_15, gf1_42, hd_40, \
                         hf_62, ip0_13, ip1_13, id_58, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_10 * gf0_15[k]
                   - f_11 * gf1_42[k]
                   + pa_y[k] * hf_62[k];

        t_159[k] = f_1 * ip0_13[k]
                   - f_2 * ip1_13[k]
                   + pb_x[k] * id_58[k];

        t_160[k] = f_5 * hd_40[k]
                   + pb_z[k] * id_58[k];

        t_161[k] = pb_x[k] * id_59[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_z, pb_x, pb_y, pb_z, gf0_12, gf1_34, \
                         hd_41, hd_45, hf_59, id_59, id_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pb_x[k] * id_60[k];

        t_163[k] = f_12 * gf0_12[k]
                   - f_13 * gf1_34[k]
                   + pa_z[k] * hf_59[k];

        t_164[k] = f_5 * hd_41[k]
                   + pb_z[k] * id_59[k];

        t_165[k] = f_5 * hd_45[k]
                   + pb_y[k] * id_60[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_y, pb_x, pb_z, gf0_16, gf1_46, hd_43, \
                         hf_69, ip0_14, ip1_14, id_61, id_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * gf0_16[k]
                   - f_13 * gf1_46[k]
                   + pa_y[k] * hf_69[k];

        t_167[k] = f_1 * ip0_14[k]
                   - f_2 * ip1_14[k]
                   + pb_x[k] * id_61[k];

        t_168[k] = f_9 * hd_43[k]
                   + pb_z[k] * id_61[k];

        t_169[k] = pb_x[k] * id_62[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_z, pb_x, pb_y, pb_z, gf0_13, gf1_39, \
                         hd_44, hd_47, hf_66, id_62, id_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_x[k] * id_63[k];

        t_171[k] = f_10 * gf0_13[k]
                   - f_11 * gf1_39[k]
                   + pa_z[k] * hf_66[k];

        t_172[k] = f_9 * hd_44[k]
                   + pb_z[k] * id_62[k];

        t_173[k] = f_8 * hd_47[k]
                   + pb_y[k] * id_63[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pa_y, pb_z, gf0_17, gf1_54, hd_46, \
                         hd_49, hf_73, hf_74, hf_75, hf_78, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * gf0_17[k]
                   - f_7 * gf1_54[k]
                   + pa_y[k] * hf_73[k];

        t_175[k] = pa_y[k] * hf_74[k];

        t_176[k] = pa_y[k] * hf_75[k];

        t_177[k] = f_5 * hd_49[k]
                   + pa_y[k] * hf_78[k];

        t_178[k] = f_4 * hd_46[k]
                   + pb_z[k] * id_64[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pa_y, pb_x, pb_y, pb_z, hd_48, hd_50, \
                         hf_81, ip0_15, ip1_15, id_65, id_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_3 * hd_50[k]
                   + pb_y[k] * id_65[k];

        t_180[k] = pa_y[k] * hf_81[k];

        t_181[k] = f_1 * ip0_15[k]
                   - f_2 * ip1_15[k]
                   + pb_x[k] * id_66[k];

        t_182[k] = f_0 * hd_48[k]
                   + pb_z[k] * id_66[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pb_x, pb_y, pb_z, hd_49, ip0_16, \
                         ip1_16, id_67, id_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pb_x[k] * id_67[k];

        t_184[k] = pb_x[k] * id_68[k];

        t_185[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_67[k];

        t_186[k] = f_0 * hd_49[k]
                   + pb_z[k] * id_67[k];

        t_187[k] = pb_y[k] * id_68[k];
    }

#pragma omp simd aligned(t_188, pb_z, hd_50, ip0_17, ip1_17, id_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_0 * hd_50[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_68[k];
    }
}

auto
compute_prim_if_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_29 = buffer.data(hd + 29);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_4, hf_6, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_4[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_5, gf1_5, hf_11, ip0_3, ip1_3, id_4, \
                         id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_5[k]
                 + pa_x[k] * hf_11[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_8, hf_7, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_1, \
                         gf1_8, hf_8, hf_19, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_19[k];

        t_18[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_8[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_9, gf1_9, hd_10, hf_23, \
                         ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_9[k]
                  + pa_x[k] * hf_23[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_14, hf_14, \
                         ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_14[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_14[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, gf1_3, \
                         gf1_10, hf_20, hf_31, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_31[k];

        t_30[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_20[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_11, gf1_11, hd_15, \
                         hf_32, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_15[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_x[k] * hf_32[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_16, hf_26, \
                         ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_26[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_16[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_17, gf1_17, hf_33, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_33[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_18, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_18[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_11, gf1_11, hf_40, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_15, gf1_15, hd_22, \
                         hf_46, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_22[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_15[k]
                  + pa_y[k] * hf_46[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_12, gf0_16, \
                         gf1_12, gf1_16, hd_25, hf_44, hf_52, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_z[k] * hf_44[k];

        t_58[k] = f_10 * hd_25[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_16[k]
                  + pa_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_13, gf1_13, hf_50, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_13[k]
                  + pa_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_17, gf1_17, hd_26, \
                         hf_53, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_26[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_y[k] * hf_53[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_29, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_29[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_39 = buffer.data(gf1 + 39);
    const auto *gf1_44 = buffer.data(gf1 + 44);
    const auto *gf1_49 = buffer.data(gf1 + 49);
    const auto *gf1_51 = buffer.data(gf1 + 51);
    const auto *gf1_55 = buffer.data(gf1 + 55);
    const auto *gf1_62 = buffer.data(gf1 + 62);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_79 = buffer.data(hf + 79);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_7, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_7[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_5, gf1_17, hf_17, ip0_3, ip1_3, \
                         id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_17[k]
                 + pa_x[k] * hf_17[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_10, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_10[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_7, \
                         gf1_27, hf_14, hf_27, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_27[k]
                  + pa_x[k] * hf_27[k];

        t_18[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_7[k]
                  + pa_y[k] * hf_14[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_9, gf1_30, hd_12, \
                         hf_31, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_30[k]
                  + pa_x[k] * hf_31[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_2, gf1_10, hd_16, \
                         hf_22, ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_10[k]
                  + pa_z[k] * hf_22[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, \
                         gf1_14, gf1_34, hf_28, hf_44, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_34[k]
                  + pa_x[k] * hf_44[k];

        t_30[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_14[k]
                  + pa_y[k] * hf_28[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_11, gf1_39, hd_17, \
                         hf_47, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_39[k]
                  + pa_x[k] * hf_47[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_6, gf1_22, hd_18, \
                         hf_39, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_22[k]
                  + pa_z[k] * hf_39[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_17, gf1_62, hf_52, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_62[k]
                  + pa_x[k] * hf_52[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_11, gf1_39, hf_62, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_39[k]
                  + pa_z[k] * hf_62[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_15, gf1_51, hd_25, \
                         hf_69, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_51[k]
                  + pa_y[k] * hf_69[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_12, gf0_16, \
                         gf1_44, gf1_55, hd_28, hf_67, hf_75, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_44[k]
                  + pa_z[k] * hf_67[k];

        t_58[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_55[k]
                  + pa_y[k] * hf_75[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_13, gf1_49, hf_73, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_49[k]
                  + pa_z[k] * hf_73[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_17, gf1_62, hd_29, \
                         hf_79, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_62[k]
                  + pa_y[k] * hf_79[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
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

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_39 = buffer.data(gf0 + 39);
    const auto *gf0_44 = buffer.data(gf0 + 44);
    const auto *gf0_49 = buffer.data(gf0 + 49);
    const auto *gf0_51 = buffer.data(gf0 + 51);
    const auto *gf0_55 = buffer.data(gf0 + 55);
    const auto *gf0_62 = buffer.data(gf0 + 62);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_41 = buffer.data(gf1 + 41);
    const auto *gf1_43 = buffer.data(gf1 + 43);
    const auto *gf1_46 = buffer.data(gf1 + 46);
    const auto *gf1_53 = buffer.data(gf1 + 53);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_68 = buffer.data(hf + 68);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pa_z, pb_z, hd_1, hf_0, hf_3, \
                         hf_5, ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_8[k] = pa_y[k] * hf_5[k];

        t_9[k] = pa_z[k] * hf_0[k];

        t_10[k] = pa_z[k] * hf_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_y, pa_z, pb_x, pb_z, gf0_0, gf1_0, hd_2, \
                         hd_7, hf_5, hf_6, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_5[k];

        t_12[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_y[k] * hf_6[k];

        t_13[k] = pb_z[k] * id_3[k];

        t_14[k] = f_6 * hd_7[k]
                  + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_z, pb_z, gf0_17, gf1_14, hf_7, \
                         hf_13, ip0_3, ip1_3, id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * gf0_17[k]
                  - f_8 * gf1_14[k]
                  + pa_x[k] * hf_13[k];

        t_16[k] = pb_z[k] * id_4[k];

        t_17[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];

        t_18[k] = pa_z[k] * hf_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_11, \
                         hf_8, hf_9, id_6, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * hf_9[k];

        t_20[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_8[k];

        t_21[k] = pb_y[k] * id_6[k];

        t_22[k] = f_6 * hd_11[k]
                  + pb_x[k] * id_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_y, gf0_27, gf1_22, hf_21, ip0_4, ip1_4, \
                         id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];

        t_24[k] = pb_y[k] * id_8[k];

        t_25[k] = f_7 * gf0_27[k]
                  - f_8 * gf1_22[k]
                  + pa_x[k] * hf_21[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, pb_z, gf0_7, gf1_7, hd_13, hf_10, id_9, \
                         id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_9 * gf0_7[k]
                  - f_10 * gf1_7[k]
                  + pa_y[k] * hf_10[k];

        t_27[k] = pb_z[k] * id_9[k];

        t_28[k] = f_3 * hd_13[k]
                  + pb_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, gf0_30, gf1_25, hf_10, \
                         hf_25, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * gf0_30[k]
                  - f_10 * gf1_25[k]
                  + pa_x[k] * hf_25[k];

        t_30[k] = pb_z[k] * id_10[k];

        t_31[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];

        t_32[k] = pa_z[k] * hf_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, gf0_10, gf1_9, hd_8, hd_10, \
                         hf_13, hf_15, hf_16, hf_19, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * hf_13[k];

        t_34[k] = f_3 * hd_8[k]
                  + pa_z[k] * hf_15[k];

        t_35[k] = f_3 * hd_10[k]
                  + pa_y[k] * hf_19[k];

        t_36[k] = pa_y[k] * hf_21[k];

        t_37[k] = f_9 * gf0_10[k]
                  - f_10 * gf1_9[k]
                  + pa_z[k] * hf_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, hd_17, ip0_6, ip1_6, id_12, \
                         id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * id_12[k];

        t_39[k] = f_3 * hd_17[k]
                  + pb_x[k] * id_14[k];

        t_40[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];

        t_41[k] = pb_y[k] * id_14[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pa_y, pb_z, gf0_14, gf0_34, gf1_11, gf1_28, \
                         hf_22, hf_34, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_9 * gf0_34[k]
                  - f_10 * gf1_28[k]
                  + pa_x[k] * hf_34[k];

        t_43[k] = f_7 * gf0_14[k]
                  - f_8 * gf1_11[k]
                  + pa_y[k] * hf_22[k];

        t_44[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, gf0_39, gf1_33, hd_18, \
                         hf_36, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_16[k];

        t_46[k] = f_4 * gf0_39[k]
                  - f_5 * gf1_33[k]
                  + pa_x[k] * hf_36[k];

        t_47[k] = pb_z[k] * id_16[k];

        t_48[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pa_z, gf0_22, gf1_17, hd_14, hf_22, \
                         hf_25, hf_27, hf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_z[k] * hf_22[k];

        t_50[k] = pa_z[k] * hf_25[k];

        t_51[k] = f_3 * hd_14[k]
                  + pa_z[k] * hf_27[k];

        t_52[k] = f_4 * gf0_22[k]
                  - f_5 * gf1_17[k]
                  + pa_y[k] * hf_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pa_y, gf0_49, gf0_51, gf1_41, gf1_43, \
                         hd_16, hf_32, hf_34, hf_37, hf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * gf0_49[k]
                  - f_5 * gf1_41[k]
                  + pa_x[k] * hf_37[k];

        t_54[k] = f_4 * gf0_51[k]
                  - f_5 * gf1_43[k]
                  + pa_x[k] * hf_38[k];

        t_55[k] = f_3 * hd_16[k]
                  + pa_y[k] * hf_32[k];

        t_56[k] = pa_y[k] * hf_34[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_z, pb_x, pb_y, gf0_22, gf1_17, hd_19, \
                         hf_29, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_7 * gf0_22[k]
                  - f_8 * gf1_17[k]
                  + pa_z[k] * hf_29[k];

        t_58[k] = pb_y[k] * id_18[k];

        t_59[k] = f_11 * hd_19[k]
                  + pb_x[k] * id_20[k];

        t_60[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_x, pa_z, pb_y, gf0_62, gf1_53, \
                         hd_20, hf_35, hf_40, hf_41, hf_44, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_y[k] * id_20[k];

        t_62[k] = f_4 * gf0_62[k]
                  - f_5 * gf1_53[k]
                  + pa_x[k] * hf_40[k];

        t_63[k] = f_3 * hd_20[k]
                  + pa_x[k] * hf_41[k];

        t_64[k] = pa_x[k] * hf_44[k];

        t_65[k] = pa_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, hd_25, hd_28, hd_33, hf_49, hf_55, \
                         hf_63, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * hd_25[k]
                  + pa_x[k] * hf_49[k];

        t_67[k] = f_3 * hd_28[k]
                  + pa_x[k] * hf_55[k];

        t_68[k] = f_3 * hd_33[k]
                  + pa_x[k] * hf_63[k];

        t_69[k] = pa_x[k] * hf_68[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pb_x, pb_y, pb_z, hd_21, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_21, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_71[k] = pb_x[k] * id_22[k];

        t_72[k] = pb_x[k] * id_23[k];

        t_73[k] = f_0 * hd_21[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_74[k] = pb_z[k] * id_22[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_z, hd_22, hf_41, hf_44, hf_46, \
                         ip0_11, ip1_11, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];

        t_76[k] = pa_z[k] * hf_41[k];

        t_77[k] = pa_z[k] * hf_44[k];

        t_78[k] = f_3 * hd_22[k]
                  + pa_z[k] * hf_46[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_z, pb_x, gf0_39, gf1_33, hf_47, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_80[k] = pb_x[k] * id_25[k];

        t_81[k] = pb_x[k] * id_26[k];

        t_82[k] = f_4 * gf0_39[k]
                  - f_5 * gf1_33[k]
                  + pa_z[k] * hf_47[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_x, pb_y, gf0_51, gf1_43, hd_27, \
                         hf_54, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_6 * hd_27[k]
                  + pb_y[k] * id_26[k];

        t_84[k] = f_7 * gf0_51[k]
                  - f_8 * gf1_43[k]
                  + pa_y[k] * hf_54[k];

        t_85[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_86[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_x, pb_y, gf0_44, gf0_55, \
                         gf1_36, gf1_46, hd_30, hf_52, hf_60, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * id_29[k];

        t_88[k] = f_9 * gf0_44[k]
                  - f_10 * gf1_36[k]
                  + pa_z[k] * hf_52[k];

        t_89[k] = f_3 * hd_30[k]
                  + pb_y[k] * id_29[k];

        t_90[k] = f_9 * gf0_55[k]
                  - f_10 * gf1_46[k]
                  + pa_y[k] * hf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_z, pb_x, gf0_49, gf1_41, hf_58, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_92[k] = pb_x[k] * id_31[k];

        t_93[k] = pb_x[k] * id_32[k];

        t_94[k] = f_7 * gf0_49[k]
                  - f_8 * gf1_41[k]
                  + pa_z[k] * hf_58[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_y, pb_y, gf0_62, gf1_53, hd_32, hd_34, \
                         hf_62, hf_66, hf_68, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * hd_32[k]
                  + pb_y[k] * id_32[k];

        t_96[k] = f_4 * gf0_62[k]
                  - f_5 * gf1_53[k]
                  + pa_y[k] * hf_62[k];

        t_97[k] = f_3 * hd_34[k]
                  + pa_y[k] * hf_66[k];

        t_98[k] = pa_y[k] * hf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_33, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_100[k] = pb_x[k] * id_34[k];

        t_101[k] = pb_x[k] * id_35[k];

        t_102[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_34[k];

        t_103[k] = pb_y[k] * id_35[k];
    }

#pragma omp simd aligned(t_104, pb_z, hd_35, ip0_17, ip1_17, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_0 * hd_35[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_5, gf1_11, hf_11, ip0_3, ip1_3, \
                         id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_11[k]
                 + pa_x[k] * hf_11[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_6, \
                         gf1_19, hf_8, hf_19, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_19[k]
                  + pa_x[k] * hf_19[k];

        t_18[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_9, gf1_21, hd_12, \
                         hf_23, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_21[k]
                  + pa_x[k] * hf_23[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_2, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, gf1_8, \
                         gf1_23, hf_20, hf_31, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_23[k]
                  + pa_x[k] * hf_31[k];

        t_30[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_11, gf1_27, hd_17, \
                         hf_33, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_27[k]
                  + pa_x[k] * hf_33[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_6, gf1_14, hd_18, \
                         hf_26, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_14[k]
                  + pa_z[k] * hf_26[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_17, gf1_44, hf_35, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_44[k]
                  + pa_x[k] * hf_35[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_11, gf1_27, hf_42, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_27[k]
                  + pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_15, gf1_36, hd_25, \
                         hf_48, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_36[k]
                  + pa_y[k] * hf_48[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_12, gf0_16, \
                         gf1_30, gf1_38, hd_28, hf_46, hf_54, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_30[k]
                  + pa_z[k] * hf_46[k];

        t_58[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_38[k]
                  + pa_y[k] * hf_54[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_13, gf1_34, hf_52, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_34[k]
                  + pa_z[k] * hf_52[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_17, gf1_44, hd_29, \
                         hf_56, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_44[k]
                  + pa_y[k] * hf_56[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_44 = buffer.data(gf0 + 44);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_40 = buffer.data(gf1 + 40);
    const auto *gf1_43 = buffer.data(gf1 + 43);
    const auto *gf1_50 = buffer.data(gf1 + 50);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_59 = buffer.data(hf + 59);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_11, gf1_13, hf_12, ip0_3, ip1_3, \
                         id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_11[k]
                 - f_7 * gf1_13[k]
                 + pa_x[k] * hf_12[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_6, gf0_19, gf1_7, \
                         gf1_21, hf_9, hf_20, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_19[k]
                  - f_7 * gf1_21[k]
                  + pa_x[k] * hf_20[k];

        t_18[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_7[k]
                  + pa_y[k] * hf_9[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_21, gf1_23, hd_12, \
                         hf_24, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_21[k]
                  - f_9 * gf1_23[k]
                  + pa_x[k] * hf_24[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_7, gf1_8, hd_16, hf_15, \
                         ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_7[k]
                  - f_9 * gf1_8[k]
                  + pa_z[k] * hf_15[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_8, gf0_23, \
                         gf1_10, gf1_25, hf_21, hf_32, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_23[k]
                  - f_9 * gf1_25[k]
                  + pa_x[k] * hf_32[k];

        t_30[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_10[k]
                  + pa_y[k] * hf_21[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_27, gf1_30, hd_17, \
                         hf_34, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_30[k]
                  + pa_x[k] * hf_34[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_14, gf1_16, hd_18, \
                         hf_27, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_14[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_27[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_44, gf1_50, hf_36, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_44[k]
                  - f_4 * gf1_50[k]
                  + pa_x[k] * hf_36[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_27, gf1_30, hf_43, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_30[k]
                  + pa_z[k] * hf_43[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_36, gf1_40, hd_25, \
                         hf_50, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_36[k]
                  - f_7 * gf1_40[k]
                  + pa_y[k] * hf_50[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_30, gf0_38, \
                         gf1_33, gf1_43, hd_28, hf_48, hf_56, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_30[k]
                  - f_9 * gf1_33[k]
                  + pa_z[k] * hf_48[k];

        t_58[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_38[k]
                  - f_9 * gf1_43[k]
                  + pa_y[k] * hf_56[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_34, gf1_38, hf_54, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_34[k]
                  - f_7 * gf1_38[k]
                  + pa_z[k] * hf_54[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_44, gf1_50, hd_29, \
                         hf_59, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_44[k]
                  - f_4 * gf1_50[k]
                  + pa_y[k] * hf_59[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_40 = buffer.data(gf0 + 40);
    const auto *gf0_43 = buffer.data(gf0 + 43);
    const auto *gf0_50 = buffer.data(gf0 + 50);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_13, gf1_11, hf_11, ip0_3, ip1_3, \
                         id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_13[k]
                 - f_7 * gf1_11[k]
                 + pa_x[k] * hf_11[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_7, gf0_21, gf1_6, \
                         gf1_19, hf_8, hf_19, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_21[k]
                  - f_7 * gf1_19[k]
                  + pa_x[k] * hf_19[k];

        t_18[k] = f_8 * gf0_7[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_23, gf1_21, hd_12, \
                         hf_23, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_23[k]
                  - f_9 * gf1_21[k]
                  + pa_x[k] * hf_23[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_8, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_8[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_10, gf0_25, \
                         gf1_8, gf1_23, hf_20, hf_31, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_25[k]
                  - f_9 * gf1_23[k]
                  + pa_x[k] * hf_31[k];

        t_30[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_30, gf1_27, hd_17, \
                         hf_32, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_30[k]
                  - f_4 * gf1_27[k]
                  + pa_x[k] * hf_32[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_16, gf1_14, hd_18, \
                         hf_26, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_16[k]
                  - f_7 * gf1_14[k]
                  + pa_z[k] * hf_26[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_50, gf1_44, hf_33, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_50[k]
                  - f_4 * gf1_44[k]
                  + pa_x[k] * hf_33[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_30, gf1_27, hf_40, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_30[k]
                  - f_4 * gf1_27[k]
                  + pa_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_40, gf1_36, hd_25, \
                         hf_46, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_40[k]
                  - f_7 * gf1_36[k]
                  + pa_y[k] * hf_46[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_33, gf0_43, \
                         gf1_30, gf1_38, hd_28, hf_44, hf_52, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_33[k]
                  - f_9 * gf1_30[k]
                  + pa_z[k] * hf_44[k];

        t_58[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_43[k]
                  - f_9 * gf1_38[k]
                  + pa_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_38, gf1_34, hf_50, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_38[k]
                  - f_7 * gf1_34[k]
                  + pa_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_50, gf1_44, hd_29, \
                         hf_53, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_50[k]
                  - f_4 * gf1_44[k]
                  + pa_y[k] * hf_53[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_44 = buffer.data(gf0 + 44);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_11, gf1_11, hf_11, ip0_3, ip1_3, \
                         id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_11[k]
                 - f_7 * gf1_11[k]
                 + pa_x[k] * hf_11[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_6, gf0_19, gf1_6, \
                         gf1_19, hf_8, hf_19, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_19[k]
                  - f_7 * gf1_19[k]
                  + pa_x[k] * hf_19[k];

        t_18[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_21, gf1_21, hd_12, \
                         hf_23, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_21[k]
                  - f_9 * gf1_21[k]
                  + pa_x[k] * hf_23[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_7, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_7[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_8, gf0_23, gf1_8, \
                         gf1_23, hf_20, hf_31, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_23[k]
                  - f_9 * gf1_23[k]
                  + pa_x[k] * hf_31[k];

        t_30[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_27, gf1_27, hd_17, \
                         hf_33, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_27[k]
                  + pa_x[k] * hf_33[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_14, gf1_14, hd_18, \
                         hf_26, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_14[k]
                  - f_7 * gf1_14[k]
                  + pa_z[k] * hf_26[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_44, gf1_44, hf_35, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_44[k]
                  - f_4 * gf1_44[k]
                  + pa_x[k] * hf_35[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_27, gf1_27, hf_42, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_27[k]
                  + pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_36, gf1_36, hd_25, \
                         hf_48, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_36[k]
                  - f_7 * gf1_36[k]
                  + pa_y[k] * hf_48[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_30, gf0_38, \
                         gf1_30, gf1_38, hd_28, hf_46, hf_54, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_30[k]
                  - f_9 * gf1_30[k]
                  + pa_z[k] * hf_46[k];

        t_58[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_38[k]
                  - f_9 * gf1_38[k]
                  + pa_y[k] * hf_54[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_34, gf1_34, hf_52, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_34[k]
                  - f_7 * gf1_34[k]
                  + pa_z[k] * hf_52[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_44, gf1_44, hd_29, \
                         hf_56, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_44[k]
                  - f_4 * gf1_44[k]
                  + pa_y[k] * hf_56[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gf0, const size_t gf1,
                                     const size_t hd, const size_t hf, const size_t ip0,
                                     const size_t ip1, const size_t id, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_21 = buffer.data(gf0 + 21);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_36 = buffer.data(gf0 + 36);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_44 = buffer.data(gf0 + 44);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_36 = buffer.data(gf1 + 36);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, \
                         ip0_2, ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_11, gf1_11, hf_11, ip0_3, ip1_3, \
                         id_4, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_11[k]
                 - f_7 * gf1_11[k]
                 + pa_x[k] * hf_11[k];

        t_10[k] = pb_z[k] * id_4[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_6[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_8[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_6, gf0_19, gf1_6, \
                         gf1_19, hf_8, hf_19, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * gf0_19[k]
                  - f_7 * gf1_19[k]
                  + pa_x[k] * hf_19[k];

        t_18[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_19[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_21, gf1_21, hd_12, \
                         hf_23, ip0_5, ip1_5, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_10[k];

        t_21[k] = f_8 * gf0_21[k]
                  - f_9 * gf1_21[k]
                  + pa_x[k] * hf_23[k];

        t_22[k] = pb_z[k] * id_10[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_7, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_7[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_25[k] = pb_y[k] * id_12[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_14[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_8, gf0_23, gf1_8, \
                         gf1_23, hf_20, hf_31, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_8 * gf0_23[k]
                  - f_9 * gf1_23[k]
                  + pa_x[k] * hf_31[k];

        t_30[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_31[k] = pb_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_27, gf1_27, hd_17, \
                         hf_32, ip0_7, ip1_7, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_27[k]
                  + pa_x[k] * hf_32[k];

        t_34[k] = pb_z[k] * id_16[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_14, gf1_14, hd_18, \
                         hf_26, ip0_8, ip1_8, id_18, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_14[k]
                  - f_7 * gf1_14[k]
                  + pa_z[k] * hf_26[k];

        t_37[k] = pb_y[k] * id_18[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_20[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_44, gf1_44, hf_33, \
                         ip0_9, ip1_9, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_3 * gf0_44[k]
                  - f_4 * gf1_44[k]
                  + pa_x[k] * hf_33[k];

        t_42[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_21[k];

        t_43[k] = pb_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * id_23[k];

        t_45[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_22[k];

        t_46[k] = pb_z[k] * id_22[k];

        t_47[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, gf0_27, gf1_27, hf_40, ip0_12, \
                         ip1_12, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = f_3 * gf0_27[k]
                  - f_4 * gf1_27[k]
                  + pa_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, gf0_36, gf1_36, hd_25, \
                         hf_46, ip0_13, ip1_13, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = f_6 * gf0_36[k]
                  - f_7 * gf1_36[k]
                  + pa_y[k] * hf_46[k];

        t_54[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_27[k];

        t_55[k] = pb_x[k] * id_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pa_z, pb_x, pb_y, gf0_30, gf0_38, \
                         gf1_30, gf1_38, hd_28, hf_44, hf_52, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * id_29[k];

        t_57[k] = f_8 * gf0_30[k]
                  - f_9 * gf1_30[k]
                  + pa_z[k] * hf_44[k];

        t_58[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_29[k];

        t_59[k] = f_8 * gf0_38[k]
                  - f_9 * gf1_38[k]
                  + pa_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, gf0_34, gf1_34, hf_50, ip0_14, \
                         ip1_14, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];

        t_62[k] = pb_x[k] * id_32[k];

        t_63[k] = f_6 * gf0_34[k]
                  - f_7 * gf1_34[k]
                  + pa_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, gf0_44, gf1_44, hd_29, \
                         hf_53, ip0_15, ip1_15, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_32[k];

        t_65[k] = f_3 * gf0_44[k]
                  - f_4 * gf1_44[k]
                  + pa_y[k] * hf_53[k];

        t_66[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_35[k];

        t_69[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = pb_y[k] * id_35[k];

        t_71[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);

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

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_17 = buffer.data(ip1 + 17);
    const auto *ip1_21 = buffer.data(ip1 + 21);
    const auto *ip1_25 = buffer.data(ip1 + 25);
    const auto *ip1_26 = buffer.data(ip1 + 26);
    const auto *ip1_27 = buffer.data(ip1 + 27);
    const auto *ip1_29 = buffer.data(ip1 + 29);
    const auto *ip1_31 = buffer.data(ip1 + 31);
    const auto *ip1_33 = buffer.data(ip1 + 33);
    const auto *ip1_36 = buffer.data(ip1 + 36);
    const auto *ip1_37 = buffer.data(ip1 + 37);
    const auto *ip1_38 = buffer.data(ip1 + 38);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, hd_0, hd_1, hd_2, ip0_0, ip0_1, \
                         ip1_0, ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = f_0 * hd_1[k]
                 + pb_x[k] * id_1[k];

        t_2[k] = f_0 * hd_2[k]
                 + pb_x[k] * id_2[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, hd_0, hd_4, hf_0, ip0_2, \
                         ip1_2, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_5[k] = pa_y[k] * hf_0[k];

        t_6[k] = f_3 * hd_0[k]
                 + pb_y[k] * id_3[k];

        t_7[k] = f_4 * hd_4[k]
                 + pb_x[k] * id_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pa_z, pb_x, pb_z, hd_0, hd_1, hd_6, hf_0, \
                         hf_1, id_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * hd_1[k]
                 + pa_y[k] * hf_1[k];

        t_9[k] = pa_z[k] * hf_0[k];

        t_10[k] = f_3 * hd_0[k]
                  + pb_z[k] * id_5[k];

        t_11[k] = f_4 * hd_6[k]
                  + pb_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pb_y, gf0_0, gf1_0, hd_2, hd_3, hf_2, \
                         hf_3, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * hd_2[k]
                  + pa_z[k] * hf_2[k];

        t_13[k] = f_6 * gf0_0[k]
                  - f_7 * gf1_0[k]
                  + pa_y[k] * hf_3[k];

        t_14[k] = f_8 * hd_3[k]
                  + pb_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, pb_z, gf0_4, gf1_4, hd_8, hf_6, ip0_3, \
                         ip1_6, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_16[k] = f_10 * gf0_4[k]
                  - f_11 * gf1_4[k]
                  + pa_x[k] * hf_6[k];

        t_17[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_6[k]
                  + pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_x, pb_z, gf0_0, gf1_0, hd_5, hd_12, hf_4, \
                         id_10, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_6 * gf0_0[k]
                  - f_7 * gf1_0[k]
                  + pa_z[k] * hf_4[k];

        t_19[k] = f_8 * hd_5[k]
                  + pb_z[k] * id_10[k];

        t_20[k] = f_9 * hd_12[k]
                  + pb_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, gf0_1, gf0_6, gf1_1, gf1_6, hf_5, \
                         hf_8, ip0_4, ip1_8, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_11[k];

        t_22[k] = f_10 * gf0_6[k]
                  - f_11 * gf1_6[k]
                  + pa_x[k] * hf_8[k];

        t_23[k] = f_12 * gf0_1[k]
                  - f_13 * gf1_1[k]
                  + pa_y[k] * hf_5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_y, gf0_7, gf1_7, hd_7, hd_14, hf_10, \
                         id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd_7[k]
                  + pb_y[k] * id_13[k];

        t_25[k] = f_5 * hd_14[k]
                  + pb_x[k] * id_14[k];

        t_26[k] = f_12 * gf0_7[k]
                  - f_13 * gf1_7[k]
                  + pa_x[k] * hf_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_z, gf0_2, gf1_2, hd_10, hf_7, \
                         ip0_5, ip1_11, id_15, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_15[k];

        t_28[k] = pa_y[k] * hf_7[k];

        t_29[k] = f_12 * gf0_2[k]
                  - f_13 * gf1_2[k]
                  + pa_z[k] * hf_7[k];

        t_30[k] = f_5 * hd_10[k]
                  + pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pb_x, pb_y, gf0_8, gf1_8, hd_19, hf_13, \
                         ip0_6, ip1_14, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * hd_19[k]
                  + pb_x[k] * id_19[k];

        t_32[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_14[k]
                  + pb_y[k] * id_18[k];

        t_33[k] = f_12 * gf0_8[k]
                  - f_13 * gf1_8[k]
                  + pa_x[k] * hf_13[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_x, pb_y, gf0_3, gf1_3, hd_13, hd_21, hf_9, \
                         id_20, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_10 * gf0_3[k]
                  - f_11 * gf1_3[k]
                  + pa_y[k] * hf_9[k];

        t_35[k] = f_9 * hd_13[k]
                  + pb_y[k] * id_20[k];

        t_36[k] = f_8 * hd_21[k]
                  + pb_x[k] * id_21[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pa_y, pb_z, gf0_5, gf0_9, gf1_5, gf1_9, \
                         hf_11, hf_14, ip0_7, ip1_17, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_6 * gf0_9[k]
                  - f_7 * gf1_9[k]
                  + pa_x[k] * hf_14[k];

        t_38[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_22[k];

        t_39[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_y[k] * hf_11[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pa_z, gf0_5, gf0_11, gf0_12, \
                         gf1_5, gf1_11, gf1_12, hf_12, hf_15, hf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_6 * gf0_11[k]
                  - f_7 * gf1_11[k]
                  + pa_x[k] * hf_15[k];

        t_41[k] = f_6 * gf0_12[k]
                  - f_7 * gf1_12[k]
                  + pa_x[k] * hf_16[k];

        t_42[k] = pa_y[k] * hf_12[k];

        t_43[k] = f_10 * gf0_5[k]
                  - f_11 * gf1_5[k]
                  + pa_z[k] * hf_12[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_y, pb_z, hd_17, hd_25, ip0_8, ip1_21, \
                         id_27, id_28, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * hd_17[k]
                  + pb_z[k] * id_27[k];

        t_45[k] = f_8 * hd_25[k]
                  + pb_x[k] * id_29[k];

        t_46[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_21[k]
                  + pb_y[k] * id_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, pb_y, gf0_14, gf1_14, hd_20, \
                         hd_26, hd_27, hf_17, hf_18, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * gf0_14[k]
                  - f_7 * gf1_14[k]
                  + pa_x[k] * hf_17[k];

        t_48[k] = f_5 * hd_26[k]
                  + pa_x[k] * hf_18[k];

        t_49[k] = f_4 * hd_20[k]
                  + pb_y[k] * id_30[k];

        t_50[k] = f_3 * hd_27[k]
                  + pb_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pa_x, hd_43, hf_19, hf_22, hf_23, \
                         hf_24, hf_25, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_x[k] * hf_19[k];

        t_52[k] = pa_x[k] * hf_22[k];

        t_53[k] = pa_x[k] * hf_23[k];

        t_54[k] = pa_x[k] * hf_24[k];

        t_55[k] = pa_x[k] * hf_25[k];

        t_56[k] = f_5 * hd_43[k]
                  + pa_x[k] * hf_27[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_x, pb_z, hd_24, hd_45, hf_29, ip0_9, \
                         ip1_25, id_36, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_4 * hd_24[k]
                  + pb_z[k] * id_36[k];

        t_58[k] = f_3 * hd_45[k]
                  + pb_x[k] * id_37[k];

        t_59[k] = pa_x[k] * hf_29[k];

        t_60[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_25[k]
                  + pb_x[k] * id_38[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pb_z, hd_26, hd_27, hd_28, ip0_10, \
                         ip0_11, ip1_26, ip1_27, id_38, id_39, id_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * hd_26[k]
                  + pb_y[k] * id_38[k];

        t_62[k] = f_0 * hd_27[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_26[k]
                  + pb_y[k] * id_39[k];

        t_63[k] = f_0 * hd_28[k]
                  + pb_y[k] * id_40[k];

        t_64[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_27[k]
                  + pb_z[k] * id_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_y, pb_z, hd_27, hd_28, hd_31, hf_19, \
                         hf_20, id_41, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_z[k] * hf_19[k];

        t_66[k] = f_3 * hd_27[k]
                  + pb_z[k] * id_41[k];

        t_67[k] = f_4 * hd_31[k]
                  + pb_y[k] * id_43[k];

        t_68[k] = f_5 * hd_28[k]
                  + pa_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_z, pb_x, pb_z, gf0_9, gf1_9, hd_29, hf_21, \
                         ip0_12, ip1_29, id_44, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_29[k]
                  + pb_x[k] * id_44[k];

        t_70[k] = f_6 * gf0_9[k]
                  - f_7 * gf1_9[k]
                  + pa_z[k] * hf_21[k];

        t_71[k] = f_8 * hd_29[k]
                  + pb_z[k] * id_45[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pb_x, pb_y, gf0_12, gf1_12, hd_35, hf_23, \
                         ip0_13, ip1_31, id_47, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * hd_35[k]
                  + pb_y[k] * id_47[k];

        t_73[k] = f_10 * gf0_12[k]
                  - f_11 * gf1_12[k]
                  + pa_y[k] * hf_23[k];

        t_74[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_31[k]
                  + pb_x[k] * id_48[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_y, pb_z, gf0_10, gf1_10, hd_33, hd_39, \
                         hf_22, id_49, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_12 * gf0_10[k]
                  - f_13 * gf1_10[k]
                  + pa_z[k] * hf_22[k];

        t_76[k] = f_5 * hd_33[k]
                  + pb_z[k] * id_49[k];

        t_77[k] = f_5 * hd_39[k]
                  + pb_y[k] * id_51[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_y, pa_z, pb_x, gf0_11, gf0_13, gf1_11, gf1_13, \
                         hf_24, hf_25, ip0_14, ip1_33, id_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_12 * gf0_13[k]
                  - f_13 * gf1_13[k]
                  + pa_y[k] * hf_25[k];

        t_79[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_33[k]
                  + pb_x[k] * id_52[k];

        t_80[k] = f_10 * gf0_11[k]
                  - f_11 * gf1_11[k]
                  + pa_z[k] * hf_24[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pb_y, pb_z, gf0_14, gf1_14, hd_37, \
                         hd_42, hd_44, hf_26, hf_28, id_53, id_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_9 * hd_37[k]
                  + pb_z[k] * id_53[k];

        t_82[k] = f_8 * hd_42[k]
                  + pb_y[k] * id_55[k];

        t_83[k] = f_6 * gf0_14[k]
                  - f_7 * gf1_14[k]
                  + pa_y[k] * hf_26[k];

        t_84[k] = f_5 * hd_44[k]
                  + pa_y[k] * hf_28[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pb_x, pb_y, pb_z, hd_40, hd_45, hf_29, \
                         ip0_15, ip1_36, id_56, id_58, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_4 * hd_40[k]
                  + pb_z[k] * id_56[k];

        t_86[k] = f_3 * hd_45[k]
                  + pb_y[k] * id_58[k];

        t_87[k] = pa_y[k] * hf_29[k];

        t_88[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_36[k]
                  + pb_x[k] * id_59[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, hd_43, hd_44, hd_45, ip0_16, \
                         ip0_17, ip1_37, ip1_38, id_59, id_60, id_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * hd_43[k]
                  + pb_z[k] * id_59[k];

        t_90[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_37[k]
                  + pb_y[k] * id_60[k];

        t_91[k] = f_0 * hd_44[k]
                  + pb_z[k] * id_60[k];

        t_92[k] = f_0 * hd_45[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_38[k]
                  + pb_z[k] * id_61[k];
    }
}

auto
compute_prim_if_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.5 / p;

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

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_4 = buffer.data(hd + 4);
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
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
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
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
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
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_71 = buffer.data(hf + 71);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_z[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, hd_1, hd_2, hf_0, hf_3, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * id_2[k];

        t_6[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_7[k] = pa_y[k] * hf_0[k];

        t_8[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_9[k] = f_4 * hd_2[k]
                 + pb_y[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_z, hd_0, hd_1, hf_0, \
                         hf_3, hf_5, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * hf_5[k];

        t_11[k] = pa_z[k] * hf_0[k];

        t_12[k] = f_4 * hd_0[k]
                  + pb_z[k] * id_6[k];

        t_13[k] = pa_z[k] * hf_3[k];

        t_14[k] = f_4 * hd_1[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pa_z, pb_y, pb_z, gf0_0, gf1_0, hd_2, \
                         hf_5, hf_6, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_5[k];

        t_17[k] = f_5 * gf0_0[k]
                  - f_6 * gf1_0[k]
                  + pa_y[k] * hf_6[k];

        t_18[k] = pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_y, pb_z, gf0_5, gf1_5, hd_5, \
                         hd_10, hf_13, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_20[k] = f_8 * gf0_5[k]
                  - f_9 * gf1_5[k]
                  + pa_x[k] * hf_13[k];

        t_21[k] = pb_z[k] * id_10[k];

        t_22[k] = f_10 * hd_5[k]
                  + pb_y[k] * id_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, pb_z, hd_4, hf_7, hf_9, ip0_3, \
                         ip1_3, id_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_11[k];

        t_24[k] = pa_y[k] * hf_9[k];

        t_25[k] = pa_z[k] * hf_7[k];

        t_26[k] = f_4 * hd_4[k]
                  + pb_z[k] * id_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, gf0_0, gf1_0, hd_8, hf_8, \
                         hf_10, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * hd_8[k]
                  + pb_y[k] * id_13[k];

        t_28[k] = pa_y[k] * hf_10[k];

        t_29[k] = f_5 * gf0_0[k]
                  - f_6 * gf1_0[k]
                  + pa_z[k] * hf_8[k];

        t_30[k] = pb_y[k] * id_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, hd_6, hd_7, hd_16, \
                         ip0_4, ip1_4, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * hd_6[k]
                  + pb_z[k] * id_14[k];

        t_32[k] = f_7 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_33[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_15[k];

        t_34[k] = f_10 * hd_7[k]
                  + pb_z[k] * id_15[k];

        t_35[k] = pb_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pa_y, pb_z, gf0_1, gf0_8, gf1_1, gf1_8, \
                         hf_11, hf_19, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_8 * gf0_8[k]
                  - f_9 * gf1_8[k]
                  + pa_x[k] * hf_19[k];

        t_37[k] = f_11 * gf0_1[k]
                  - f_12 * gf1_1[k]
                  + pa_y[k] * hf_11[k];

        t_38[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_x, pb_y, pb_z, gf0_10, gf1_10, \
                         hd_11, hd_18, hf_22, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * hd_18[k]
                  + pb_x[k] * id_18[k];

        t_40[k] = f_11 * gf0_10[k]
                  - f_12 * gf1_10[k]
                  + pa_x[k] * hf_22[k];

        t_41[k] = pb_z[k] * id_18[k];

        t_42[k] = f_3 * hd_11[k]
                  + pb_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_z, pb_z, hd_9, hd_10, hf_11, hf_13, \
                         ip0_5, ip1_5, id_19, id_20, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_19[k];

        t_44[k] = pa_z[k] * hf_11[k];

        t_45[k] = f_4 * hd_9[k]
                  + pb_z[k] * id_20[k];

        t_46[k] = pa_z[k] * hf_13[k];

        t_47[k] = f_4 * hd_10[k]
                  + pb_z[k] * id_21[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, hd_11, hd_13, hd_15, \
                         hf_14, hf_15, hf_16, hf_18, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * hd_13[k]
                  + pb_y[k] * id_22[k];

        t_49[k] = f_3 * hd_11[k]
                  + pa_z[k] * hf_14[k];

        t_50[k] = pa_y[k] * hf_15[k];

        t_51[k] = pa_y[k] * hf_16[k];

        t_52[k] = f_3 * hd_15[k]
                  + pa_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, gf0_2, gf1_2, hd_12, \
                         hd_16, hf_15, hf_19, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * hd_12[k]
                  + pb_z[k] * id_23[k];

        t_54[k] = f_4 * hd_16[k]
                  + pb_y[k] * id_24[k];

        t_55[k] = pa_y[k] * hf_19[k];

        t_56[k] = f_11 * gf0_2[k]
                  - f_12 * gf1_2[k]
                  + pa_z[k] * hf_15[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pb_x, pb_y, pb_z, hd_14, hd_15, hd_28, \
                         ip0_6, ip1_6, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * id_25[k];

        t_58[k] = f_3 * hd_14[k]
                  + pb_z[k] * id_25[k];

        t_59[k] = f_3 * hd_28[k]
                  + pb_x[k] * id_27[k];

        t_60[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_26[k];

        t_61[k] = f_3 * hd_15[k]
                  + pb_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_3, \
                         gf1_12, hf_20, hf_29, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * id_27[k];

        t_63[k] = f_11 * gf0_12[k]
                  - f_12 * gf1_12[k]
                  + pa_x[k] * hf_29[k];

        t_64[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_3[k]
                  + pa_y[k] * hf_20[k];

        t_65[k] = pb_z[k] * id_28[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pb_x, pb_y, pb_z, gf0_13, gf1_13, \
                         hd_19, hd_30, hf_32, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_10 * hd_30[k]
                  + pb_x[k] * id_29[k];

        t_67[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_13[k]
                  + pa_x[k] * hf_32[k];

        t_68[k] = pb_z[k] * id_29[k];

        t_69[k] = f_7 * hd_19[k]
                  + pb_y[k] * id_30[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_z, hd_17, hd_18, hf_20, hf_22, \
                         ip0_7, ip1_7, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_30[k];

        t_71[k] = pa_z[k] * hf_20[k];

        t_72[k] = f_4 * hd_17[k]
                  + pb_z[k] * id_31[k];

        t_73[k] = pa_z[k] * hf_22[k];

        t_74[k] = f_4 * hd_18[k]
                  + pb_z[k] * id_32[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pa_z, pb_y, gf0_6, gf1_6, hd_19, hd_22, \
                         hf_23, hf_24, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * hd_22[k]
                  + pb_y[k] * id_33[k];

        t_76[k] = f_3 * hd_19[k]
                  + pa_z[k] * hf_23[k];

        t_77[k] = f_5 * gf0_6[k]
                  - f_6 * gf1_6[k]
                  + pa_y[k] * hf_24[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_y, pb_z, gf0_15, gf1_15, hd_20, \
                         hd_21, hd_25, hf_33, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * hd_20[k]
                  + pb_z[k] * id_34[k];

        t_79[k] = f_5 * gf0_15[k]
                  - f_6 * gf1_15[k]
                  + pa_x[k] * hf_33[k];

        t_80[k] = f_10 * hd_21[k]
                  + pb_z[k] * id_35[k];

        t_81[k] = f_10 * hd_25[k]
                  + pb_y[k] * id_36[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, gf0_17, gf1_17, hd_27, hf_25, \
                         hf_26, hf_28, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_5 * gf0_17[k]
                  - f_6 * gf1_17[k]
                  + pa_x[k] * hf_34[k];

        t_83[k] = pa_y[k] * hf_25[k];

        t_84[k] = pa_y[k] * hf_26[k];

        t_85[k] = f_3 * hd_27[k]
                  + pa_y[k] * hf_28[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pa_z, pb_y, pb_z, gf0_6, gf1_6, hd_24, \
                         hd_28, hf_25, hf_29, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * hd_24[k]
                  + pb_z[k] * id_37[k];

        t_87[k] = f_4 * hd_28[k]
                  + pb_y[k] * id_38[k];

        t_88[k] = pa_y[k] * hf_29[k];

        t_89[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_6[k]
                  + pa_z[k] * hf_25[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pb_x, pb_y, pb_z, hd_26, hd_27, hd_36, \
                         ip0_8, ip1_8, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_y[k] * id_39[k];

        t_91[k] = f_7 * hd_26[k]
                  + pb_z[k] * id_39[k];

        t_92[k] = f_10 * hd_36[k]
                  + pb_x[k] * id_41[k];

        t_93[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_40[k];

        t_94[k] = f_7 * hd_27[k]
                  + pb_z[k] * id_40[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_x, pb_x, pb_y, gf0_20, gf1_20, hd_37, \
                         hd_38, hf_38, hf_39, id_41, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_y[k] * id_41[k];

        t_96[k] = f_5 * gf0_20[k]
                  - f_6 * gf1_20[k]
                  + pa_x[k] * hf_38[k];

        t_97[k] = f_3 * hd_37[k]
                  + pa_x[k] * hf_39[k];

        t_98[k] = f_4 * hd_38[k]
                  + pb_x[k] * id_43[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, pa_x, pa_z, pb_z, hd_29, \
                         hf_30, hf_42, hf_44, hf_45, hf_47, id_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * hf_42[k];

        t_100[k] = pa_x[k] * hf_44[k];

        t_101[k] = pa_x[k] * hf_45[k];

        t_102[k] = pa_z[k] * hf_30[k];

        t_103[k] = f_4 * hd_29[k]
                   + pb_z[k] * id_44[k];

        t_104[k] = pa_x[k] * hf_47[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, t_110, pa_x, pb_z, hd_31, hd_43, \
                         hf_48, hf_49, hf_50, hf_51, hf_52, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_x[k] * hf_48[k];

        t_106[k] = pa_x[k] * hf_49[k];

        t_107[k] = f_3 * hd_43[k]
                   + pa_x[k] * hf_50[k];

        t_108[k] = f_10 * hd_31[k]
                   + pb_z[k] * id_45[k];

        t_109[k] = pa_x[k] * hf_51[k];

        t_110[k] = pa_x[k] * hf_52[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, t_116, pa_x, pb_z, hd_32, hd_46, \
                         hf_53, hf_54, hf_55, hf_56, hf_57, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_x[k] * hf_53[k];

        t_112[k] = pa_x[k] * hf_54[k];

        t_113[k] = f_3 * hd_46[k]
                   + pa_x[k] * hf_55[k];

        t_114[k] = f_3 * hd_32[k]
                   + pb_z[k] * id_46[k];

        t_115[k] = pa_x[k] * hf_56[k];

        t_116[k] = pa_x[k] * hf_57[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, pa_x, pa_y, hf_35, \
                         hf_36, hf_58, hf_59, hf_60, hf_61, hf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_x[k] * hf_58[k];

        t_118[k] = pa_x[k] * hf_59[k];

        t_119[k] = pa_y[k] * hf_35[k];

        t_120[k] = pa_y[k] * hf_36[k];

        t_121[k] = pa_x[k] * hf_60[k];

        t_122[k] = pa_x[k] * hf_61[k];

        t_123[k] = pa_x[k] * hf_62[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_x, pb_x, pb_z, hd_35, hd_51, \
                         hd_53, hf_64, hf_68, hf_69, id_47, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_3 * hd_51[k]
                   + pa_x[k] * hf_64[k];

        t_125[k] = f_13 * hd_35[k]
                   + pb_z[k] * id_47[k];

        t_126[k] = f_4 * hd_53[k]
                   + pb_x[k] * id_48[k];

        t_127[k] = pa_x[k] * hf_68[k];

        t_128[k] = pa_x[k] * hf_69[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pa_x, pb_x, pb_z, hf_71, ip0_9, \
                         ip1_9, id_49, id_50, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pa_x[k] * hf_71[k];

        t_130[k] = f_1 * ip0_9[k]
                   - f_2 * ip1_9[k]
                   + pb_x[k] * id_49[k];

        t_131[k] = pb_z[k] * id_49[k];

        t_132[k] = pb_x[k] * id_50[k];

        t_133[k] = pb_x[k] * id_51[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pb_y, pb_z, hd_38, hd_39, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_50, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * hd_38[k]
                   + f_1 * ip0_10[k]
                   - f_2 * ip1_10[k]
                   + pb_y[k] * id_50[k];

        t_135[k] = pb_z[k] * id_50[k];

        t_136[k] = f_0 * hd_39[k]
                   + pb_y[k] * id_51[k];

        t_137[k] = f_1 * ip0_11[k]
                   - f_2 * ip1_11[k]
                   + pb_z[k] * id_51[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pa_z, pb_x, pb_z, hd_37, hd_38, \
                         hf_39, hf_42, id_52, id_53, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pa_z[k] * hf_39[k];

        t_139[k] = f_4 * hd_37[k]
                   + pb_z[k] * id_52[k];

        t_140[k] = pb_x[k] * id_54[k];

        t_141[k] = pa_z[k] * hf_42[k];

        t_142[k] = f_4 * hd_38[k]
                   + pb_z[k] * id_53[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_z, pb_x, pb_y, pb_z, hd_39, hd_40, \
                         hd_42, hf_45, ip0_12, ip1_12, id_54, id_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_13 * hd_42[k]
                   + pb_y[k] * id_54[k];

        t_144[k] = f_3 * hd_39[k]
                   + pa_z[k] * hf_45[k];

        t_145[k] = f_1 * ip0_12[k]
                   - f_2 * ip1_12[k]
                   + pb_x[k] * id_55[k];

        t_146[k] = f_10 * hd_40[k]
                   + pb_z[k] * id_55[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pa_z, pb_x, pb_y, pb_z, gf0_13, \
                         gf1_13, hd_41, hd_45, hf_46, id_56, id_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pb_x[k] * id_56[k];

        t_148[k] = pb_x[k] * id_57[k];

        t_149[k] = f_5 * gf0_13[k]
                   - f_6 * gf1_13[k]
                   + pa_z[k] * hf_46[k];

        t_150[k] = f_10 * hd_41[k]
                   + pb_z[k] * id_56[k];

        t_151[k] = f_7 * hd_45[k]
                   + pb_y[k] * id_57[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_y, pb_x, pb_z, gf0_17, gf1_17, hd_43, \
                         hf_54, ip0_13, ip1_13, id_58, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_8 * gf0_17[k]
                   - f_9 * gf1_17[k]
                   + pa_y[k] * hf_54[k];

        t_153[k] = f_1 * ip0_13[k]
                   - f_2 * ip1_13[k]
                   + pb_x[k] * id_58[k];

        t_154[k] = f_3 * hd_43[k]
                   + pb_z[k] * id_58[k];

        t_155[k] = pb_x[k] * id_59[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pa_z, pb_x, pb_y, pb_z, gf0_14, gf1_14, \
                         hd_44, hd_48, hf_51, id_59, id_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_x[k] * id_60[k];

        t_157[k] = f_11 * gf0_14[k]
                   - f_12 * gf1_14[k]
                   + pa_z[k] * hf_51[k];

        t_158[k] = f_3 * hd_44[k]
                   + pb_z[k] * id_59[k];

        t_159[k] = f_3 * hd_48[k]
                   + pb_y[k] * id_60[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_y, pb_x, pb_z, gf0_19, gf1_19, hd_46, \
                         hf_59, ip0_14, ip1_14, id_61, id_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_11 * gf0_19[k]
                   - f_12 * gf1_19[k]
                   + pa_y[k] * hf_59[k];

        t_161[k] = f_1 * ip0_14[k]
                   - f_2 * ip1_14[k]
                   + pb_x[k] * id_61[k];

        t_162[k] = f_7 * hd_46[k]
                   + pb_z[k] * id_61[k];

        t_163[k] = pb_x[k] * id_62[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pb_x, pb_y, pb_z, gf0_15, gf1_15, \
                         hd_47, hd_50, hf_56, id_62, id_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pb_x[k] * id_63[k];

        t_165[k] = f_8 * gf0_15[k]
                   - f_9 * gf1_15[k]
                   + pa_z[k] * hf_56[k];

        t_166[k] = f_7 * hd_47[k]
                   + pb_z[k] * id_62[k];

        t_167[k] = f_10 * hd_50[k]
                   + pb_y[k] * id_63[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pa_y, pb_x, gf0_20, gf1_20, hd_52, \
                         hf_63, hf_64, hf_66, hf_68, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_5 * gf0_20[k]
                   - f_6 * gf1_20[k]
                   + pa_y[k] * hf_63[k];

        t_169[k] = pa_y[k] * hf_64[k];

        t_170[k] = pa_y[k] * hf_66[k];

        t_171[k] = pb_x[k] * id_64[k];

        t_172[k] = f_3 * hd_52[k]
                   + pa_y[k] * hf_68[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_y, pb_x, pb_y, pb_z, hd_49, hd_53, \
                         hf_71, ip0_15, ip1_15, id_64, id_65, id_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_13 * hd_49[k]
                   + pb_z[k] * id_64[k];

        t_174[k] = f_4 * hd_53[k]
                   + pb_y[k] * id_65[k];

        t_175[k] = pa_y[k] * hf_71[k];

        t_176[k] = f_1 * ip0_15[k]
                   - f_2 * ip1_15[k]
                   + pb_x[k] * id_66[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, pb_x, pb_y, pb_z, hd_51, \
                         hd_52, ip0_16, ip1_16, id_66, id_67, id_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pb_y[k] * id_66[k];

        t_178[k] = f_0 * hd_51[k]
                   + pb_z[k] * id_66[k];

        t_179[k] = pb_x[k] * id_67[k];

        t_180[k] = pb_x[k] * id_68[k];

        t_181[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_67[k];

        t_182[k] = f_0 * hd_52[k]
                   + pb_z[k] * id_67[k];
    }

#pragma omp simd aligned(t_183, t_184, pb_y, pb_z, hd_53, ip0_17, ip1_17, \
                         id_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pb_y[k] * id_68[k];

        t_184[k] = f_0 * hd_53[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_68[k];
    }
}

auto
compute_prim_if_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.5 / p;

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

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_26 = buffer.data(gf1 + 26);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hd_0, hd_1, hf_0, hf_3, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_8[k] = pa_z[k] * hf_0[k];

        t_9[k] = f_4 * hd_0[k]
                 + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, pb_z, gf0_0, gf1_0, hd_2, \
                         hd_8, hf_5, hf_6, id_7, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_5[k];

        t_11[k] = f_5 * gf0_0[k]
                  - f_6 * gf1_0[k]
                  + pa_y[k] * hf_6[k];

        t_12[k] = pb_z[k] * id_7[k];

        t_13[k] = f_7 * hd_8[k]
                  + pb_x[k] * id_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_z, gf0_5, gf1_7, hf_10, ip0_3, ip1_3, \
                         id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * gf0_5[k]
                  - f_9 * gf1_7[k]
                  + pa_x[k] * hf_10[k];

        t_15[k] = pb_z[k] * id_8[k];

        t_16[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, pb_z, gf0_0, gf1_0, hd_5, \
                         hd_12, hf_7, id_10, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * gf0_0[k]
                  - f_6 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_10 * hd_5[k]
                  + pb_z[k] * id_10[k];

        t_20[k] = f_7 * hd_12[k]
                  + pb_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pb_y, gf0_8, gf1_10, hf_13, ip0_4, ip1_4, \
                         id_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_11[k];

        t_22[k] = pb_y[k] * id_12[k];

        t_23[k] = f_8 * gf0_8[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_x, pb_z, gf0_1, gf1_3, hd_14, hf_8, id_13, \
                         id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_11 * gf0_1[k]
                  - f_12 * gf1_3[k]
                  + pa_y[k] * hf_8[k];

        t_25[k] = pb_z[k] * id_13[k];

        t_26[k] = f_3 * hd_14[k]
                  + pb_x[k] * id_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pb_z, gf0_10, gf1_12, hf_11, \
                         hf_16, ip0_5, ip1_5, id_14, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * gf0_10[k]
                  - f_12 * gf1_12[k]
                  + pa_x[k] * hf_16[k];

        t_28[k] = pb_z[k] * id_14[k];

        t_29[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_15[k];

        t_30[k] = pa_y[k] * hf_11[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_z, pb_x, pb_y, pb_z, gf0_2, gf1_4, hd_10, \
                         hd_19, hf_11, id_17, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_11 * gf0_2[k]
                  - f_12 * gf1_4[k]
                  + pa_z[k] * hf_11[k];

        t_32[k] = pb_y[k] * id_17[k];

        t_33[k] = f_3 * hd_10[k]
                  + pb_z[k] * id_17[k];

        t_34[k] = f_3 * hd_19[k]
                  + pb_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_y, gf0_12, gf1_14, hf_20, ip0_6, ip1_6, \
                         id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_18[k];

        t_36[k] = pb_y[k] * id_19[k];

        t_37[k] = f_11 * gf0_12[k]
                  - f_12 * gf1_14[k]
                  + pa_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_x, pb_z, gf0_3, gf1_5, hd_21, hf_14, \
                         id_20, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_5[k]
                  + pa_y[k] * hf_14[k];

        t_39[k] = pb_z[k] * id_20[k];

        t_40[k] = f_10 * hd_21[k]
                  + pb_x[k] * id_21[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_x, pb_z, gf0_13, gf1_16, hf_21, ip0_7, ip1_7, \
                         id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_16[k]
                  + pa_x[k] * hf_21[k];

        t_42[k] = pb_z[k] * id_21[k];

        t_43[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, gf0_6, gf0_15, gf0_17, gf1_8, \
                         gf1_19, gf1_21, hf_17, hf_18, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * gf0_6[k]
                  - f_6 * gf1_8[k]
                  + pa_y[k] * hf_17[k];

        t_45[k] = f_5 * gf0_15[k]
                  - f_6 * gf1_19[k]
                  + pa_x[k] * hf_22[k];

        t_46[k] = f_5 * gf0_17[k]
                  - f_6 * gf1_21[k]
                  + pa_x[k] * hf_23[k];

        t_47[k] = pa_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, pb_z, gf0_6, gf1_8, hd_17, \
                         hd_25, hf_18, id_27, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_8[k]
                  + pa_z[k] * hf_18[k];

        t_49[k] = pb_y[k] * id_27[k];

        t_50[k] = f_7 * hd_17[k]
                  + pb_z[k] * id_27[k];

        t_51[k] = f_10 * hd_25[k]
                  + pb_x[k] * id_29[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_y, gf0_20, gf1_26, hd_26, hf_24, \
                         hf_25, ip0_8, ip1_8, id_28, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_28[k];

        t_53[k] = pb_y[k] * id_29[k];

        t_54[k] = f_5 * gf0_20[k]
                  - f_6 * gf1_26[k]
                  + pa_x[k] * hf_24[k];

        t_55[k] = f_3 * hd_26[k]
                  + pa_x[k] * hf_25[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, pa_x, hd_39, hf_28, hf_32, hf_34, \
                         hf_35, hf_37, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * hf_28[k];

        t_57[k] = pa_x[k] * hf_32[k];

        t_58[k] = pa_x[k] * hf_34[k];

        t_59[k] = pa_x[k] * hf_35[k];

        t_60[k] = pa_x[k] * hf_37[k];

        t_61[k] = f_3 * hd_39[k]
                  + pa_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_x, pb_x, pb_z, hd_24, hf_44, ip0_9, \
                         ip1_9, id_36, id_38, id_39, id_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_13 * hd_24[k]
                  + pb_z[k] * id_36[k];

        t_63[k] = pa_x[k] * hf_44[k];

        t_64[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_38[k];

        t_65[k] = pb_x[k] * id_39[k];

        t_66[k] = pb_x[k] * id_40[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_y, pb_z, hd_27, hd_28, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_39, id_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * hd_27[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_39[k];

        t_68[k] = pb_z[k] * id_39[k];

        t_69[k] = f_0 * hd_28[k]
                  + pb_y[k] * id_40[k];

        t_70[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_40[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_z, pb_y, pb_z, hd_27, hd_28, hd_30, hf_28, \
                         hf_30, id_41, id_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_z[k] * hf_28[k];

        t_72[k] = f_4 * hd_27[k]
                  + pb_z[k] * id_41[k];

        t_73[k] = f_13 * hd_30[k]
                  + pb_y[k] * id_42[k];

        t_74[k] = f_3 * hd_28[k]
                  + pa_z[k] * hf_30[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_x, gf0_13, gf1_16, hf_31, ip0_12, \
                         ip1_12, id_43, id_44, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_43[k];

        t_76[k] = pb_x[k] * id_44[k];

        t_77[k] = pb_x[k] * id_45[k];

        t_78[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_16[k]
                  + pa_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pa_y, pb_y, pb_z, gf0_17, gf1_21, hd_29, hd_33, \
                         hf_34, id_44, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_10 * hd_29[k]
                  + pb_z[k] * id_44[k];

        t_80[k] = f_7 * hd_33[k]
                  + pb_y[k] * id_45[k];

        t_81[k] = f_8 * gf0_17[k]
                  - f_9 * gf1_21[k]
                  + pa_y[k] * hf_34[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_z, pb_x, gf0_14, gf1_18, hf_32, ip0_13, \
                         ip1_13, id_46, id_47, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_46[k];

        t_83[k] = pb_x[k] * id_47[k];

        t_84[k] = pb_x[k] * id_48[k];

        t_85[k] = f_11 * gf0_14[k]
                  - f_12 * gf1_18[k]
                  + pa_z[k] * hf_32[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_y, pb_y, pb_z, gf0_19, gf1_23, hd_32, hd_36, \
                         hf_37, id_47, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * hd_32[k]
                  + pb_z[k] * id_47[k];

        t_87[k] = f_3 * hd_36[k]
                  + pb_y[k] * id_48[k];

        t_88[k] = f_11 * gf0_19[k]
                  - f_12 * gf1_23[k]
                  + pa_y[k] * hf_37[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_z, pb_x, gf0_15, gf1_19, hf_35, ip0_14, \
                         ip1_14, id_49, id_50, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_49[k];

        t_90[k] = pb_x[k] * id_50[k];

        t_91[k] = pb_x[k] * id_51[k];

        t_92[k] = f_8 * gf0_15[k]
                  - f_9 * gf1_19[k]
                  + pa_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_y, pb_z, gf0_20, gf1_26, hd_35, \
                         hd_38, hd_40, hf_38, hf_42, id_50, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_7 * hd_35[k]
                  + pb_z[k] * id_50[k];

        t_94[k] = f_10 * hd_38[k]
                  + pb_y[k] * id_51[k];

        t_95[k] = f_5 * gf0_20[k]
                  - f_6 * gf1_26[k]
                  + pa_y[k] * hf_38[k];

        t_96[k] = f_3 * hd_40[k]
                  + pa_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pb_x, pb_y, pb_z, hd_37, hd_41, hf_44, \
                         ip0_15, ip1_15, id_52, id_53, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * hd_37[k]
                  + pb_z[k] * id_52[k];

        t_98[k] = f_4 * hd_41[k]
                  + pb_y[k] * id_53[k];

        t_99[k] = pa_y[k] * hf_44[k];

        t_100[k] = f_1 * ip0_15[k]
                   - f_2 * ip1_15[k]
                   + pb_x[k] * id_54[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, pb_x, pb_y, pb_z, hd_39, \
                         hd_40, ip0_16, ip1_16, id_54, id_55, id_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * hd_39[k]
                   + pb_z[k] * id_54[k];

        t_102[k] = pb_x[k] * id_55[k];

        t_103[k] = pb_x[k] * id_56[k];

        t_104[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_55[k];

        t_105[k] = f_0 * hd_40[k]
                   + pb_z[k] * id_55[k];

        t_106[k] = pb_y[k] * id_56[k];
    }

#pragma omp simd aligned(t_107, pb_z, hd_41, ip0_17, ip1_17, id_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * hd_41[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_56[k];
    }
}

auto
compute_prim_if_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_89 = buffer.data(hf + 89);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_7, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_7[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_5, hd_6, hf_17, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_17[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_10, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_10[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_1, \
                         gf1_8, hf_14, hf_27, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_27[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_14[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_9, gf1_10, hd_12, \
                         hf_31, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_31[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_16, hf_22, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_22[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, gf1_3, \
                         gf1_12, hf_28, hf_45, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_45[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_28[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_11, gf1_13, hd_17, \
                         hf_48, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pa_x[k] * hf_48[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_18, hf_40, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_40[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_y, gf0_17, gf1_20, hf_55, hf_60, \
                         hf_89, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_20[k]
                  + pa_x[k] * hf_55[k];

        t_44[k] = pa_x[k] * hf_60[k];

        t_45[k] = pa_x[k] * hf_89[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, hd_20, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_47[k] = pb_x[k] * id_26[k];

        t_48[k] = pb_x[k] * id_27[k];

        t_49[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_50[k] = pb_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_x, pb_z, hf_60, ip0_11, ip0_12, \
                         ip1_11, ip1_12, id_27, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_52[k] = pa_z[k] * hf_60[k];

        t_53[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_54[k] = pb_x[k] * id_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_x, pb_y, gf0_11, gf0_15, \
                         gf1_13, gf1_17, hd_25, hf_65, hf_72, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * id_31[k];

        t_56[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pa_z[k] * hf_65[k];

        t_57[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_58[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_17[k]
                  + pa_y[k] * hf_72[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, gf0_12, gf1_14, hf_70, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_60[k] = pb_x[k] * id_33[k];

        t_61[k] = pb_x[k] * id_34[k];

        t_62[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_14[k]
                  + pa_z[k] * hf_70[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pb_x, pb_y, gf0_16, gf1_19, hd_28, \
                         hf_78, ip0_14, ip1_14, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_64[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_19[k]
                  + pa_y[k] * hf_78[k];

        t_65[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_66[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, pb_y, gf0_13, gf0_17, \
                         gf1_15, gf1_20, hd_29, hf_76, hf_82, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * id_37[k];

        t_68[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_15[k]
                  + pa_z[k] * hf_76[k];

        t_69[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_70[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_20[k]
                  + pa_y[k] * hf_82[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pb_x, pb_y, hf_89, ip0_15, \
                         ip0_16, ip1_15, ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_y[k] * hf_89[k];

        t_72[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_73[k] = pb_x[k] * id_40[k];

        t_74[k] = pb_x[k] * id_41[k];

        t_75[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];
    }

#pragma omp simd aligned(t_76, t_77, pb_y, pb_z, hd_32, ip0_17, ip1_17, \
                         id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * id_41[k];

        t_77[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_33 = buffer.data(gf1 + 33);
    const auto *gf1_35 = buffer.data(gf1 + 35);
    const auto *gf1_38 = buffer.data(gf1 + 38);
    const auto *gf1_44 = buffer.data(gf1 + 44);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_74 = buffer.data(hf + 74);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_z[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, hd_1, hf_0, hf_3, hf_6, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * id_2[k];

        t_6[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_7[k] = pa_y[k] * hf_0[k];

        t_8[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_9[k] = pa_y[k] * hf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_y, gf0_0, gf1_0, hd_2, \
                         hf_0, hf_3, hf_6, hf_7, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * hf_0[k];

        t_11[k] = pa_z[k] * hf_3[k];

        t_12[k] = pb_y[k] * id_5[k];

        t_13[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_6[k];

        t_14[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_y[k] * hf_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, gf0_5, gf1_12, hd_7, hf_14, \
                         id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * id_6[k];

        t_16[k] = f_6 * hd_7[k]
                  + pb_x[k] * id_7[k];

        t_17[k] = f_7 * gf0_5[k]
                  - f_8 * gf1_12[k]
                  + pa_x[k] * hf_14[k];

        t_18[k] = pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_8, hf_9, \
                         hf_10, ip0_3, ip1_3, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_8[k];

        t_20[k] = pa_z[k] * hf_8[k];

        t_21[k] = pa_y[k] * hf_10[k];

        t_22[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, hd_11, ip0_4, ip1_4, id_9, id_10, \
                         id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_y[k] * id_9[k];

        t_24[k] = f_6 * hd_11[k]
                  + pb_x[k] * id_11[k];

        t_25[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_10[k];

        t_26[k] = pb_y[k] * id_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, gf0_1, gf0_8, gf1_6, gf1_17, \
                         hf_11, hf_22, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_17[k]
                  + pa_x[k] * hf_22[k];

        t_28[k] = f_9 * gf0_1[k]
                  - f_10 * gf1_6[k]
                  + pa_y[k] * hf_11[k];

        t_29[k] = pb_z[k] * id_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, gf0_10, gf1_20, hd_13, \
                         hf_26, ip0_5, ip1_5, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hd_13[k]
                  + pb_x[k] * id_13[k];

        t_31[k] = f_9 * gf0_10[k]
                  - f_10 * gf1_20[k]
                  + pa_x[k] * hf_26[k];

        t_32[k] = pb_z[k] * id_13[k];

        t_33[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pa_z, hd_8, hd_10, hf_11, hf_14, \
                         hf_16, hf_20, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * hf_11[k];

        t_35[k] = pa_z[k] * hf_14[k];

        t_36[k] = f_3 * hd_8[k]
                  + pa_z[k] * hf_16[k];

        t_37[k] = f_3 * hd_10[k]
                  + pa_y[k] * hf_20[k];

        t_38[k] = pa_y[k] * hf_22[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_y, gf0_2, gf1_8, hd_17, hf_17, \
                         ip0_6, ip1_6, id_15, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * gf0_2[k]
                  - f_10 * gf1_8[k]
                  + pa_z[k] * hf_17[k];

        t_40[k] = pb_y[k] * id_15[k];

        t_41[k] = f_3 * hd_17[k]
                  + pb_x[k] * id_17[k];

        t_42[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, \
                         gf1_10, gf1_23, hf_23, hf_35, id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * id_17[k];

        t_44[k] = f_9 * gf0_12[k]
                  - f_10 * gf1_23[k]
                  + pa_x[k] * hf_35[k];

        t_45[k] = f_7 * gf0_3[k]
                  - f_8 * gf1_10[k]
                  + pa_y[k] * hf_23[k];

        t_46[k] = pb_z[k] * id_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, pb_z, gf0_13, gf1_27, hd_18, \
                         hf_38, ip0_7, ip1_7, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_19[k];

        t_48[k] = f_4 * gf0_13[k]
                  - f_5 * gf1_27[k]
                  + pa_x[k] * hf_38[k];

        t_49[k] = pb_z[k] * id_19[k];

        t_50[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_20[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pa_z, gf0_6, gf1_14, hd_14, hf_23, \
                         hf_26, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * hf_23[k];

        t_52[k] = pa_z[k] * hf_26[k];

        t_53[k] = f_3 * hd_14[k]
                  + pa_z[k] * hf_28[k];

        t_54[k] = f_4 * gf0_6[k]
                  - f_5 * gf1_14[k]
                  + pa_y[k] * hf_29[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pa_y, gf0_15, gf0_17, gf1_33, gf1_35, \
                         hd_16, hf_33, hf_35, hf_39, hf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * gf0_15[k]
                  - f_5 * gf1_33[k]
                  + pa_x[k] * hf_39[k];

        t_56[k] = f_4 * gf0_17[k]
                  - f_5 * gf1_35[k]
                  + pa_x[k] * hf_40[k];

        t_57[k] = f_3 * hd_16[k]
                  + pa_y[k] * hf_33[k];

        t_58[k] = pa_y[k] * hf_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, pb_y, gf0_6, gf1_14, hd_19, \
                         hf_30, ip0_8, ip1_8, id_21, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * gf0_6[k]
                  - f_8 * gf1_14[k]
                  + pa_z[k] * hf_30[k];

        t_60[k] = pb_y[k] * id_21[k];

        t_61[k] = f_11 * hd_19[k]
                  + pb_x[k] * id_23[k];

        t_62[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_22[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_x, pb_x, pb_y, gf0_20, gf1_44, hd_20, \
                         hd_21, hf_43, hf_44, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_y[k] * id_23[k];

        t_64[k] = f_4 * gf0_20[k]
                  - f_5 * gf1_44[k]
                  + pa_x[k] * hf_43[k];

        t_65[k] = f_3 * hd_20[k]
                  + pa_x[k] * hf_44[k];

        t_66[k] = f_12 * hd_21[k]
                  + pb_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pa_x, pa_z, hd_25, hd_28, hd_33, hf_36, \
                         hf_48, hf_53, hf_59, hf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pa_x[k] * hf_48[k];

        t_68[k] = pa_z[k] * hf_36[k];

        t_69[k] = f_3 * hd_25[k]
                  + pa_x[k] * hf_53[k];

        t_70[k] = f_3 * hd_28[k]
                  + pa_x[k] * hf_59[k];

        t_71[k] = f_3 * hd_33[k]
                  + pa_x[k] * hf_68[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_x, pb_x, pb_z, hd_35, hf_74, ip0_9, \
                         ip1_9, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_12 * hd_35[k]
                  + pb_x[k] * id_25[k];

        t_73[k] = pa_x[k] * hf_74[k];

        t_74[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_26[k];

        t_75[k] = pb_z[k] * id_26[k];

        t_76[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pb_x, pb_y, pb_z, hd_21, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_x[k] * id_28[k];

        t_78[k] = f_0 * hd_21[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_27[k];

        t_79[k] = pb_z[k] * id_27[k];

        t_80[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_28[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_z, pb_x, hd_22, hf_44, hf_48, hf_50, \
                         ip0_12, ip1_12, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_z[k] * hf_44[k];

        t_82[k] = pb_x[k] * id_30[k];

        t_83[k] = pa_z[k] * hf_48[k];

        t_84[k] = f_3 * hd_22[k]
                  + pa_z[k] * hf_50[k];

        t_85[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_z, pb_x, pb_y, gf0_13, gf1_27, hd_27, \
                         hf_51, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * id_32[k];

        t_87[k] = pb_x[k] * id_33[k];

        t_88[k] = f_4 * gf0_13[k]
                  - f_5 * gf1_27[k]
                  + pa_z[k] * hf_51[k];

        t_89[k] = f_6 * hd_27[k]
                  + pb_y[k] * id_33[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_y, pb_x, gf0_17, gf1_35, hf_58, ip0_13, \
                         ip1_13, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_7 * gf0_17[k]
                  - f_8 * gf1_35[k]
                  + pa_y[k] * hf_58[k];

        t_91[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_34[k];

        t_92[k] = pb_x[k] * id_35[k];

        t_93[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pa_z, pb_y, gf0_14, gf0_19, gf1_30, gf1_38, \
                         hd_30, hf_56, hf_64, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * gf0_14[k]
                  - f_10 * gf1_30[k]
                  + pa_z[k] * hf_56[k];

        t_95[k] = f_3 * hd_30[k]
                  + pb_y[k] * id_36[k];

        t_96[k] = f_9 * gf0_19[k]
                  - f_10 * gf1_38[k]
                  + pa_y[k] * hf_64[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_z, pb_x, gf0_15, gf1_33, hf_62, ip0_14, \
                         ip1_14, id_37, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_37[k];

        t_98[k] = pb_x[k] * id_38[k];

        t_99[k] = pb_x[k] * id_39[k];

        t_100[k] = f_7 * gf0_15[k]
                   - f_8 * gf1_33[k]
                   + pa_z[k] * hf_62[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pb_x, pb_y, gf0_20, gf1_44, hd_32, \
                         hd_34, hf_67, hf_72, id_39, id_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_11 * hd_32[k]
                   + pb_y[k] * id_39[k];

        t_102[k] = f_4 * gf0_20[k]
                   - f_5 * gf1_44[k]
                   + pa_y[k] * hf_67[k];

        t_103[k] = pb_x[k] * id_40[k];

        t_104[k] = f_3 * hd_34[k]
                   + pa_y[k] * hf_72[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_y, pb_x, pb_y, hd_35, hf_74, \
                         ip0_15, ip1_15, id_41, id_42, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * hd_35[k]
                   + pb_y[k] * id_41[k];

        t_106[k] = pa_y[k] * hf_74[k];

        t_107[k] = f_1 * ip0_15[k]
                   - f_2 * ip1_15[k]
                   + pb_x[k] * id_42[k];

        t_108[k] = pb_y[k] * id_42[k];

        t_109[k] = pb_x[k] * id_43[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_x, pb_y, pb_z, hd_35, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_43, id_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * id_44[k];

        t_111[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_43[k];

        t_112[k] = pb_y[k] * id_44[k];

        t_113[k] = f_0 * hd_35[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_44[k];
    }
}

auto
compute_prim_if_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_27 = buffer.data(gf0 + 27);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_33 = buffer.data(gf0 + 33);
    const auto *gf0_35 = buffer.data(gf0 + 35);
    const auto *gf0_38 = buffer.data(gf0 + 38);
    const auto *gf0_44 = buffer.data(gf0 + 44);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_35 = buffer.data(gf1 + 35);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hd_1, hd_2, hf_0, hf_3, \
                         hf_5, ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_8[k] = pa_z[k] * hf_0[k];

        t_9[k] = f_3 * hd_2[k]
                 + pa_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, id_5, \
                         id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_y[k] * hf_6[k];

        t_11[k] = pb_z[k] * id_5[k];

        t_12[k] = f_6 * hd_6[k]
                  + pb_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_z, gf0_12, gf1_10, hf_10, ip0_3, ip1_3, \
                         id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * gf0_12[k]
                  - f_8 * gf1_10[k]
                  + pa_x[k] * hf_10[k];

        t_14[k] = pb_z[k] * id_6[k];

        t_15[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_17[k] = pb_y[k] * id_8[k];

        t_18[k] = f_6 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_19[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, pb_z, gf0_6, gf0_17, gf1_6, \
                         gf1_13, hf_8, hf_13, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_y[k] * id_10[k];

        t_21[k] = f_7 * gf0_17[k]
                  - f_8 * gf1_13[k]
                  + pa_x[k] * hf_13[k];

        t_22[k] = f_9 * gf0_6[k]
                  - f_10 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_23[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, gf0_20, gf1_15, hd_12, \
                         hf_16, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_25[k] = f_9 * gf0_20[k]
                  - f_10 * gf1_15[k]
                  + pa_x[k] * hf_16[k];

        t_26[k] = pb_z[k] * id_12[k];

        t_27[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pa_z, pb_x, pb_y, gf0_8, gf1_7, hd_16, \
                         hf_11, id_14, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * hf_11[k];

        t_29[k] = f_9 * gf0_8[k]
                  - f_10 * gf1_7[k]
                  + pa_z[k] * hf_11[k];

        t_30[k] = pb_y[k] * id_14[k];

        t_31[k] = f_3 * hd_16[k]
                  + pb_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pb_y, gf0_23, gf1_17, hf_20, ip0_6, ip1_6, \
                         id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];

        t_33[k] = pb_y[k] * id_16[k];

        t_34[k] = f_9 * gf0_23[k]
                  - f_10 * gf1_17[k]
                  + pa_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_x, pb_z, gf0_10, gf1_8, hd_17, hf_14, \
                         id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * gf0_10[k]
                  - f_8 * gf1_8[k]
                  + pa_y[k] * hf_14[k];

        t_36[k] = pb_z[k] * id_17[k];

        t_37[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_z, gf0_27, gf1_21, hf_21, ip0_7, ip1_7, \
                         id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_4 * gf0_27[k]
                  - f_5 * gf1_21[k]
                  + pa_x[k] * hf_21[k];

        t_39[k] = pb_z[k] * id_18[k];

        t_40[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pa_y, gf0_14, gf0_33, gf0_35, gf1_11, \
                         gf1_25, gf1_27, hf_17, hf_18, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_11[k]
                  + pa_y[k] * hf_17[k];

        t_42[k] = f_4 * gf0_33[k]
                  - f_5 * gf1_25[k]
                  + pa_x[k] * hf_22[k];

        t_43[k] = f_4 * gf0_35[k]
                  - f_5 * gf1_27[k]
                  + pa_x[k] * hf_23[k];

        t_44[k] = pa_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_z, pb_x, pb_y, gf0_14, gf1_11, hd_18, \
                         hf_18, ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_7 * gf0_14[k]
                  - f_8 * gf1_11[k]
                  + pa_z[k] * hf_18[k];

        t_46[k] = pb_y[k] * id_20[k];

        t_47[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_48[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_x, pb_y, gf0_44, gf1_35, hd_19, \
                         hf_24, hf_25, hf_28, hf_32, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_y[k] * id_22[k];

        t_50[k] = f_4 * gf0_44[k]
                  - f_5 * gf1_35[k]
                  + pa_x[k] * hf_24[k];

        t_51[k] = f_3 * hd_19[k]
                  + pa_x[k] * hf_25[k];

        t_52[k] = pa_x[k] * hf_28[k];

        t_53[k] = pa_x[k] * hf_32[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_x, hd_30, hf_34, hf_35, hf_37, \
                         hf_39, hf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * hf_34[k];

        t_55[k] = pa_x[k] * hf_35[k];

        t_56[k] = pa_x[k] * hf_37[k];

        t_57[k] = f_3 * hd_30[k]
                  + pa_x[k] * hf_39[k];

        t_58[k] = pa_x[k] * hf_44[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pb_x, pb_y, pb_z, hd_20, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_60[k] = pb_x[k] * id_26[k];

        t_61[k] = pb_x[k] * id_27[k];

        t_62[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_63[k] = pb_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_x, pb_z, hd_21, hf_28, hf_30, \
                         ip0_11, ip0_12, ip1_11, ip1_12, id_27, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_65[k] = pa_z[k] * hf_28[k];

        t_66[k] = f_3 * hd_21[k]
                  + pa_z[k] * hf_30[k];

        t_67[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_x, pb_y, gf0_27, gf1_21, hd_25, \
                         hf_31, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * id_30[k];

        t_69[k] = pb_x[k] * id_31[k];

        t_70[k] = f_4 * gf0_27[k]
                  - f_5 * gf1_21[k]
                  + pa_z[k] * hf_31[k];

        t_71[k] = f_6 * hd_25[k]
                  + pb_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pb_x, gf0_35, gf1_27, hf_34, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * gf0_35[k]
                  - f_8 * gf1_27[k]
                  + pa_y[k] * hf_34[k];

        t_73[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_74[k] = pb_x[k] * id_33[k];

        t_75[k] = pb_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_y, pa_z, pb_y, gf0_30, gf0_38, gf1_24, gf1_29, \
                         hd_28, hf_32, hf_37, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_9 * gf0_30[k]
                  - f_10 * gf1_24[k]
                  + pa_z[k] * hf_32[k];

        t_77[k] = f_3 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_78[k] = f_9 * gf0_38[k]
                  - f_10 * gf1_29[k]
                  + pa_y[k] * hf_37[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_z, pb_x, gf0_33, gf1_25, hf_35, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_80[k] = pb_x[k] * id_36[k];

        t_81[k] = pb_x[k] * id_37[k];

        t_82[k] = f_7 * gf0_33[k]
                  - f_8 * gf1_25[k]
                  + pa_z[k] * hf_35[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_y, gf0_44, gf1_35, hd_29, hd_31, \
                         hf_38, hf_42, hf_44, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_84[k] = f_4 * gf0_44[k]
                  - f_5 * gf1_35[k]
                  + pa_y[k] * hf_38[k];

        t_85[k] = f_3 * hd_31[k]
                  + pa_y[k] * hf_42[k];

        t_86[k] = pa_y[k] * hf_44[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_88[k] = pb_x[k] * id_40[k];

        t_89[k] = pb_x[k] * id_41[k];

        t_90[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_91[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_92, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_62 = buffer.data(hf + 62);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, \
                         ip0_2, ip1_2, id_2, id_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_7[k] = pb_z[k] * id_5[k];

        t_8[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, gf0_5, gf1_5, hf_11, ip0_3, ip1_3, id_6, \
                         id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_5[k]
                 + pa_x[k] * hf_11[k];

        t_10[k] = pb_z[k] * id_6[k];

        t_11[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_13[k] = pb_y[k] * id_8[k];

        t_14[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_15[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_1, \
                         gf1_8, hf_8, hf_19, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * id_10[k];

        t_17[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_19[k];

        t_18[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_8[k];

        t_19[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, gf0_10, gf1_10, hd_12, \
                         hf_23, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_21[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_23[k];

        t_22[k] = pb_z[k] * id_12[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_16, hf_14, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_14[k];

        t_25[k] = pb_y[k] * id_14[k];

        t_26[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_27[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_3, \
                         gf1_12, hf_20, hf_31, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * id_16[k];

        t_29[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_31[k];

        t_30[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_20[k];

        t_31[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, gf0_13, gf1_13, hd_17, \
                         hf_33, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_33[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_13[k]
                  + pa_x[k] * hf_33[k];

        t_34[k] = pb_z[k] * id_18[k];

        t_35[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_18, hf_26, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_26[k];

        t_37[k] = pb_y[k] * id_20[k];

        t_38[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_39[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_x, pb_y, gf0_20, gf1_20, hd_20, \
                         hf_35, hf_39, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * id_22[k];

        t_41[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_20[k]
                  + pa_x[k] * hf_35[k];

        t_42[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_43[k] = pa_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_x, pb_x, hd_32, hf_62, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_45[k] = pa_x[k] * hf_62[k];

        t_46[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_47[k] = pb_x[k] * id_26[k];

        t_48[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_y, pb_z, hd_20, ip0_10, ip0_11, ip1_10, ip1_11, \
                         id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_50[k] = pb_z[k] * id_26[k];

        t_51[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_z, pb_x, gf0_13, gf1_13, hf_42, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_53[k] = pb_x[k] * id_30[k];

        t_54[k] = pb_x[k] * id_31[k];

        t_55[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_13[k]
                  + pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pb_x, pb_y, gf0_17, gf1_17, hd_25, \
                         hf_48, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_57[k] = f_6 * gf0_17[k]
                  - f_7 * gf1_17[k]
                  + pa_y[k] * hf_48[k];

        t_58[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_59[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, pb_x, pb_y, gf0_14, gf0_19, \
                         gf1_14, gf1_19, hd_28, hf_46, hf_54, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pb_x[k] * id_34[k];

        t_61[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_14[k]
                  + pa_z[k] * hf_46[k];

        t_62[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_63[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_19[k]
                  + pa_y[k] * hf_54[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_x, gf0_15, gf1_15, hf_52, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_65[k] = pb_x[k] * id_36[k];

        t_66[k] = pb_x[k] * id_37[k];

        t_67[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_15[k]
                  + pa_z[k] * hf_52[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pb_y, gf0_20, gf1_20, hd_29, hd_32, \
                         hf_56, hf_62, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_69[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_20[k]
                  + pa_y[k] * hf_56[k];

        t_70[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_71[k] = pa_y[k] * hf_62[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_73[k] = pb_x[k] * id_40[k];

        t_74[k] = pb_x[k] * id_41[k];

        t_75[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_76[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_77, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_34 = buffer.data(gf1 + 34);
    const auto *gf1_41 = buffer.data(gf1 + 41);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_68 = buffer.data(hf + 68);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_7, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_7[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_11, hd_6, hf_13, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_11[k]
                  + pa_x[k] * hf_13[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_8, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_8[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_7, \
                         gf1_16, hf_10, hf_21, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_16[k]
                  + pa_x[k] * hf_21[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_7[k]
                  + pa_y[k] * hf_10[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_10, gf1_18, hd_12, \
                         hf_25, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_18[k]
                  + pa_x[k] * hf_25[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_8, hd_16, hf_16, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_8[k]
                  + pa_z[k] * hf_16[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_9, \
                         gf1_20, hf_22, hf_33, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_20[k]
                  + pa_x[k] * hf_33[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_9[k]
                  + pa_y[k] * hf_22[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_13, gf1_25, hd_17, \
                         hf_35, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_25[k]
                  + pa_x[k] * hf_35[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_13, hd_18, \
                         hf_28, ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_13[k]
                  + pa_z[k] * hf_28[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_20, gf1_41, hd_20, \
                         hf_37, hf_42, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_41[k]
                  + pa_x[k] * hf_37[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_42[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_68, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_68[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_42, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_13, gf1_25, hf_45, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_25[k]
                  + pa_z[k] * hf_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_17, gf1_32, hd_25, \
                         hf_52, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_17[k]
                  - f_7 * gf1_32[k]
                  + pa_y[k] * hf_52[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_14, gf0_19, \
                         gf1_28, gf1_34, hd_28, hf_50, hf_58, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_28[k]
                  + pa_z[k] * hf_50[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_34[k]
                  + pa_y[k] * hf_58[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_15, gf1_30, hf_56, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_30[k]
                  + pa_z[k] * hf_56[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_20, gf1_41, hd_29, hd_32, \
                         hf_61, hf_68, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_41[k]
                  + pa_y[k] * hf_61[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_68[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_25 = buffer.data(gf0 + 25);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_34 = buffer.data(gf0 + 34);
    const auto *gf0_41 = buffer.data(gf0 + 41);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_38 = buffer.data(gf1 + 38);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pa_z, pb_z, hd_2, hf_0, hf_5, ip0_2, ip1_2, \
                         id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * hd_2[k]
                 + pa_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, id_5, \
                         id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * gf0_0[k]
                 - f_5 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_10[k] = pb_z[k] * id_5[k];

        t_11[k] = f_6 * hd_6[k]
                  + pb_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_z, gf0_11, gf1_10, hf_11, ip0_3, ip1_3, \
                         id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gf0_11[k]
                  - f_8 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_13[k] = pb_z[k] * id_6[k];

        t_14[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_16[k] = pb_y[k] * id_8[k];

        t_17[k] = f_6 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_18[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_y, pb_y, pb_z, gf0_7, gf0_16, gf1_6, \
                         gf1_15, hf_8, hf_19, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * id_10[k];

        t_20[k] = f_7 * gf0_16[k]
                  - f_8 * gf1_15[k]
                  + pa_x[k] * hf_19[k];

        t_21[k] = f_9 * gf0_7[k]
                  - f_10 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_22[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, pb_z, gf0_18, gf1_17, hd_12, \
                         hf_23, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_24[k] = f_9 * gf0_18[k]
                  - f_10 * gf1_17[k]
                  + pa_x[k] * hf_23[k];

        t_25[k] = pb_z[k] * id_12[k];

        t_26[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_x, pb_y, gf0_8, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_9 * gf0_8[k]
                  - f_10 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_28[k] = pb_y[k] * id_14[k];

        t_29[k] = f_3 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_30[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_y, pb_y, pb_z, gf0_9, gf0_20, gf1_8, \
                         gf1_19, hf_20, hf_31, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * id_16[k];

        t_32[k] = f_9 * gf0_20[k]
                  - f_10 * gf1_19[k]
                  + pa_x[k] * hf_31[k];

        t_33[k] = f_7 * gf0_9[k]
                  - f_8 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_34[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_x, pb_z, gf0_25, gf1_23, hd_17, \
                         hf_33, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_36[k] = f_4 * gf0_25[k]
                  - f_5 * gf1_23[k]
                  + pa_x[k] * hf_33[k];

        t_37[k] = pb_z[k] * id_18[k];

        t_38[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_y, gf0_13, gf1_12, hd_18, \
                         hf_26, ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * gf0_13[k]
                  - f_8 * gf1_12[k]
                  + pa_z[k] * hf_26[k];

        t_40[k] = pb_y[k] * id_20[k];

        t_41[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_42[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, gf0_41, gf1_38, hd_20, \
                         hf_35, hf_39, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * id_22[k];

        t_44[k] = f_4 * gf0_41[k]
                  - f_5 * gf1_38[k]
                  + pa_x[k] * hf_35[k];

        t_45[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_46[k] = pa_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pb_x, hd_32, hf_62, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_48[k] = pa_x[k] * hf_62[k];

        t_49[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_50[k] = pb_x[k] * id_26[k];

        t_51[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_z, pb_y, pb_z, hd_20, hf_39, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_53[k] = pb_z[k] * id_26[k];

        t_54[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_55[k] = pa_z[k] * hf_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, hd_21, hf_41, ip0_12, ip1_12, \
                         id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * hd_21[k]
                  + pa_z[k] * hf_41[k];

        t_57[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_58[k] = pb_x[k] * id_30[k];

        t_59[k] = pb_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_y, pa_z, pb_y, gf0_25, gf0_32, gf1_23, gf1_30, \
                         hd_25, hf_42, hf_48, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * gf0_25[k]
                  - f_5 * gf1_23[k]
                  + pa_z[k] * hf_42[k];

        t_61[k] = f_6 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_62[k] = f_7 * gf0_32[k]
                  - f_8 * gf1_30[k]
                  + pa_y[k] * hf_48[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_x, gf0_28, gf1_26, hf_46, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_64[k] = pb_x[k] * id_33[k];

        t_65[k] = pb_x[k] * id_34[k];

        t_66[k] = f_9 * gf0_28[k]
                  - f_10 * gf1_26[k]
                  + pa_z[k] * hf_46[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pb_x, pb_y, gf0_34, gf1_32, hd_28, \
                         hf_54, ip0_14, ip1_14, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_68[k] = f_9 * gf0_34[k]
                  - f_10 * gf1_32[k]
                  + pa_y[k] * hf_54[k];

        t_69[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_70[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pa_z, pb_x, pb_y, gf0_30, gf0_41, \
                         gf1_28, gf1_38, hd_29, hf_52, hf_56, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_x[k] * id_37[k];

        t_72[k] = f_7 * gf0_30[k]
                  - f_8 * gf1_28[k]
                  + pa_z[k] * hf_52[k];

        t_73[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_74[k] = f_4 * gf0_41[k]
                  - f_5 * gf1_38[k]
                  + pa_y[k] * hf_56[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_x, pb_y, hd_31, hd_32, hf_60, hf_62, \
                         ip0_15, ip1_15, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * hd_31[k]
                  + pa_y[k] * hf_60[k];

        t_76[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_77[k] = pa_y[k] * hf_62[k];

        t_78[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, pb_z, hd_32, ip0_16, \
                         ip0_17, ip1_16, ip1_17, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_x[k] * id_40[k];

        t_80[k] = pb_x[k] * id_41[k];

        t_81[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_82[k] = pb_y[k] * id_41[k];

        t_83[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_38 = buffer.data(gf0 + 38);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_21 = buffer.data(gf1 + 21);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_25 = buffer.data(gf1 + 25);
    const auto *gf1_27 = buffer.data(gf1 + 27);
    const auto *gf1_29 = buffer.data(gf1 + 29);
    const auto *gf1_35 = buffer.data(gf1 + 35);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_41 = buffer.data(hf + 41);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_10, gf1_10, hd_6, \
                         hf_10, ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_10[k]
                  + pa_x[k] * hf_10[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_6, gf0_15, gf1_6, \
                         gf1_13, hf_8, hf_13, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_13[k]
                  + pa_x[k] * hf_13[k];

        t_20[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_17, gf1_15, hd_12, \
                         hf_16, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_17[k]
                  - f_9 * gf1_15[k]
                  + pa_x[k] * hf_16[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_7, gf1_7, hd_16, hf_11, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_7[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_11[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_8, gf0_19, gf1_8, \
                         gf1_17, hf_14, hf_19, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_17[k]
                  + pa_x[k] * hf_19[k];

        t_32[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_14[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_23, gf1_21, hd_17, \
                         hf_20, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_21[k]
                  + pa_x[k] * hf_20[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_12, gf1_11, hd_18, \
                         hf_17, ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_12[k]
                  - f_7 * gf1_11[k]
                  + pa_z[k] * hf_17[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_y, gf0_38, gf1_35, hf_21, hf_25, \
                         hf_41, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_38[k]
                  - f_4 * gf1_35[k]
                  + pa_x[k] * hf_21[k];

        t_44[k] = pa_x[k] * hf_25[k];

        t_45[k] = pa_x[k] * hf_41[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, hd_20, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_47[k] = pb_x[k] * id_26[k];

        t_48[k] = pb_x[k] * id_27[k];

        t_49[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_50[k] = pb_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_x, pb_z, hf_25, ip0_11, ip0_12, \
                         ip1_11, ip1_12, id_27, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_52[k] = pa_z[k] * hf_25[k];

        t_53[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_54[k] = pb_x[k] * id_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_x, pb_y, gf0_23, gf0_30, \
                         gf1_21, gf1_27, hd_25, hf_28, hf_31, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * id_31[k];

        t_56[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_21[k]
                  + pa_z[k] * hf_28[k];

        t_57[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_58[k] = f_6 * gf0_30[k]
                  - f_7 * gf1_27[k]
                  + pa_y[k] * hf_31[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, gf0_26, gf1_24, hf_29, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_60[k] = pb_x[k] * id_33[k];

        t_61[k] = pb_x[k] * id_34[k];

        t_62[k] = f_8 * gf0_26[k]
                  - f_9 * gf1_24[k]
                  + pa_z[k] * hf_29[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pb_x, pb_y, gf0_32, gf1_29, hd_28, \
                         hf_34, ip0_14, ip1_14, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_64[k] = f_8 * gf0_32[k]
                  - f_9 * gf1_29[k]
                  + pa_y[k] * hf_34[k];

        t_65[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_66[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, pb_y, gf0_28, gf0_38, \
                         gf1_25, gf1_35, hd_29, hf_32, hf_35, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * id_37[k];

        t_68[k] = f_6 * gf0_28[k]
                  - f_7 * gf1_25[k]
                  + pa_z[k] * hf_32[k];

        t_69[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_70[k] = f_3 * gf0_38[k]
                  - f_4 * gf1_35[k]
                  + pa_y[k] * hf_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pb_x, pb_y, hf_41, ip0_15, \
                         ip0_16, ip1_15, ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_y[k] * hf_41[k];

        t_72[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_73[k] = pb_x[k] * id_40[k];

        t_74[k] = pb_x[k] * id_41[k];

        t_75[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];
    }

#pragma omp simd aligned(t_76, t_77, pb_y, pb_z, hd_32, ip0_17, ip1_17, \
                         id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * id_41[k];

        t_77[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_62 = buffer.data(hf + 62);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_5, hd_6, hf_11, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_1, \
                         gf1_8, hf_8, hf_19, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_19[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_9, gf1_10, hd_12, \
                         hf_23, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_23[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_16, hf_14, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_14[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, gf1_3, \
                         gf1_12, hf_20, hf_31, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_31[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_20[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_11, gf1_13, hd_17, \
                         hf_33, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pa_x[k] * hf_33[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_18, hf_26, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_26[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_y, gf0_17, gf1_20, hf_35, hf_39, \
                         hf_62, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_20[k]
                  + pa_x[k] * hf_35[k];

        t_44[k] = pa_x[k] * hf_39[k];

        t_45[k] = pa_x[k] * hf_62[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, hd_20, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_47[k] = pb_x[k] * id_26[k];

        t_48[k] = pb_x[k] * id_27[k];

        t_49[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_50[k] = pb_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_x, pb_z, hf_39, ip0_11, ip0_12, \
                         ip1_11, ip1_12, id_27, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_52[k] = pa_z[k] * hf_39[k];

        t_53[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_54[k] = pb_x[k] * id_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_x, pb_y, gf0_11, gf0_15, \
                         gf1_13, gf1_17, hd_25, hf_42, hf_48, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * id_31[k];

        t_56[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pa_z[k] * hf_42[k];

        t_57[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_58[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_17[k]
                  + pa_y[k] * hf_48[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, gf0_12, gf1_14, hf_46, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_60[k] = pb_x[k] * id_33[k];

        t_61[k] = pb_x[k] * id_34[k];

        t_62[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_14[k]
                  + pa_z[k] * hf_46[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pb_x, pb_y, gf0_16, gf1_19, hd_28, \
                         hf_54, ip0_14, ip1_14, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_64[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_19[k]
                  + pa_y[k] * hf_54[k];

        t_65[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_66[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, pb_y, gf0_13, gf0_17, \
                         gf1_15, gf1_20, hd_29, hf_52, hf_56, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * id_37[k];

        t_68[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_15[k]
                  + pa_z[k] * hf_52[k];

        t_69[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_70[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_20[k]
                  + pa_y[k] * hf_56[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pb_x, pb_y, hf_62, ip0_15, \
                         ip0_16, ip1_15, ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_y[k] * hf_62[k];

        t_72[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_73[k] = pb_x[k] * id_40[k];

        t_74[k] = pb_x[k] * id_41[k];

        t_75[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];
    }

#pragma omp simd aligned(t_76, t_77, pb_y, pb_z, hd_32, ip0_17, ip1_17, \
                         id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * id_41[k];

        t_77[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_38 = buffer.data(gf1 + 38);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_62 = buffer.data(hf + 62);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_10, hd_6, hf_11, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_6, \
                         gf1_15, hf_8, hf_19, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_15[k]
                  + pa_x[k] * hf_19[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_10, gf1_17, hd_12, \
                         hf_23, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_17[k]
                  + pa_x[k] * hf_23[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_8, \
                         gf1_19, hf_20, hf_31, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_19[k]
                  + pa_x[k] * hf_31[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_13, gf1_23, hd_17, \
                         hf_33, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_33[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_12, hd_18, \
                         hf_26, ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_12[k]
                  + pa_z[k] * hf_26[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_20, gf1_38, hd_20, \
                         hf_35, hf_39, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_38[k]
                  + pa_x[k] * hf_35[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_62, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_62[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_39, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_13, gf1_23, hf_42, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_23[k]
                  + pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_17, gf1_30, hd_25, \
                         hf_48, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_17[k]
                  - f_7 * gf1_30[k]
                  + pa_y[k] * hf_48[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_14, gf0_19, \
                         gf1_26, gf1_32, hd_28, hf_46, hf_54, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_26[k]
                  + pa_z[k] * hf_46[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_32[k]
                  + pa_y[k] * hf_54[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_15, gf1_28, hf_52, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_28[k]
                  + pa_z[k] * hf_52[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_20, gf1_38, hd_29, hd_32, \
                         hf_56, hf_62, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_38[k]
                  + pa_y[k] * hf_56[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_62[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_23 = buffer.data(gf0 + 23);
    const auto *gf0_26 = buffer.data(gf0 + 26);
    const auto *gf0_28 = buffer.data(gf0 + 28);
    const auto *gf0_30 = buffer.data(gf0 + 30);
    const auto *gf0_32 = buffer.data(gf0 + 32);
    const auto *gf0_38 = buffer.data(gf0 + 38);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_23 = buffer.data(gf1 + 23);
    const auto *gf1_26 = buffer.data(gf1 + 26);
    const auto *gf1_28 = buffer.data(gf1 + 28);
    const auto *gf1_30 = buffer.data(gf1 + 30);
    const auto *gf1_32 = buffer.data(gf1 + 32);
    const auto *gf1_38 = buffer.data(gf1 + 38);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_62 = buffer.data(hf + 62);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_10, gf1_10, hd_6, \
                         hf_11, ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_6, gf0_15, gf1_6, \
                         gf1_15, hf_8, hf_19, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_15[k]
                  + pa_x[k] * hf_19[k];

        t_20[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_6[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_17, gf1_17, hd_12, \
                         hf_23, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_17[k]
                  - f_9 * gf1_17[k]
                  + pa_x[k] * hf_23[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_7, gf1_7, hd_16, hf_14, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_7[k]
                  - f_9 * gf1_7[k]
                  + pa_z[k] * hf_14[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_8, gf0_19, gf1_8, \
                         gf1_19, hf_20, hf_31, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_19[k]
                  + pa_x[k] * hf_31[k];

        t_32[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_y[k] * hf_20[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_23, gf1_23, hd_17, \
                         hf_33, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_33[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_12, gf1_12, hd_18, \
                         hf_26, ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_12[k]
                  - f_7 * gf1_12[k]
                  + pa_z[k] * hf_26[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_38, gf1_38, hd_20, \
                         hf_35, hf_39, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_38[k]
                  - f_4 * gf1_38[k]
                  + pa_x[k] * hf_35[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_62, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_62[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_39, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_23, gf1_23, hf_42, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_23[k]
                  + pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_30, gf1_30, hd_25, \
                         hf_48, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_30[k]
                  - f_7 * gf1_30[k]
                  + pa_y[k] * hf_48[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_26, gf0_32, \
                         gf1_26, gf1_32, hd_28, hf_46, hf_54, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_26[k]
                  - f_9 * gf1_26[k]
                  + pa_z[k] * hf_46[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_32[k]
                  - f_9 * gf1_32[k]
                  + pa_y[k] * hf_54[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_28, gf1_28, hf_52, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_28[k]
                  - f_7 * gf1_28[k]
                  + pa_z[k] * hf_52[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_38, gf1_38, hd_29, hd_32, \
                         hf_56, hf_62, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_38[k]
                  - f_4 * gf1_38[k]
                  + pa_y[k] * hf_56[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_62[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_23 = buffer.data(hd + 23);

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

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_35 = buffer.data(id + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_0, hf_0, hf_1, \
                         ip0_0, ip1_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pa_y[k] * hf_0[k];

        t_2[k] = pa_z[k] * hf_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gf0_0, gf0_4, gf0_6, gf1_0, gf1_4, gf1_6, \
                         hf_2, hf_4, hf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gf0_4[k]
                 - f_6 * gf1_4[k]
                 + pa_x[k] * hf_4[k];

        t_5[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];

        t_6[k] = f_5 * gf0_6[k]
                 - f_6 * gf1_6[k]
                 + pa_x[k] * hf_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, gf0_1, gf0_2, gf0_7, gf1_1, \
                         gf1_2, gf1_7, hf_3, hf_5, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gf0_1[k]
                 - f_8 * gf1_1[k]
                 + pa_y[k] * hf_3[k];

        t_8[k] = f_7 * gf0_7[k]
                 - f_8 * gf1_7[k]
                 + pa_x[k] * hf_8[k];

        t_9[k] = pa_y[k] * hf_5[k];

        t_10[k] = f_7 * gf0_2[k]
                  - f_8 * gf1_2[k]
                  + pa_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gf0_3, gf0_8, gf0_9, gf1_3, gf1_8, \
                         gf1_9, hf_7, hf_11, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_8[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = f_5 * gf0_3[k]
                  - f_6 * gf1_3[k]
                  + pa_y[k] * hf_7[k];

        t_13[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pa_x[k] * hf_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_y, gf0_5, gf0_11, gf0_12, gf1_5, \
                         gf1_11, gf1_12, hf_9, hf_10, hf_13, hf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pa_y[k] * hf_9[k];

        t_15[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_x[k] * hf_13[k];

        t_16[k] = f_3 * gf0_12[k]
                  - f_4 * gf1_12[k]
                  + pa_x[k] * hf_14[k];

        t_17[k] = pa_y[k] * hf_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, gf0_5, gf0_14, gf1_5, \
                         gf1_14, hf_10, hf_15, hf_16, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gf0_5[k]
                  - f_6 * gf1_5[k]
                  + pa_z[k] * hf_10[k];

        t_19[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_15[k];

        t_20[k] = pa_x[k] * hf_16[k];

        t_21[k] = pa_x[k] * hf_18[k];

        t_22[k] = pa_x[k] * hf_19[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pa_z, pb_y, hd_16, hf_16, hf_20, \
                         hf_21, hf_23, ip0_1, ip1_1, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * hf_20[k];

        t_24[k] = pa_x[k] * hf_21[k];

        t_25[k] = pa_x[k] * hf_23[k];

        t_26[k] = f_0 * hd_16[k]
                  + f_1 * ip0_1[k]
                  - f_2 * ip1_1[k]
                  + pb_y[k] * id_26[k];

        t_27[k] = pa_z[k] * hf_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, gf0_9, gf0_10, gf0_12, gf1_9, gf1_10, \
                         gf1_12, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pa_z[k] * hf_17[k];

        t_29[k] = f_5 * gf0_12[k]
                  - f_6 * gf1_12[k]
                  + pa_y[k] * hf_19[k];

        t_30[k] = f_7 * gf0_10[k]
                  - f_8 * gf1_10[k]
                  + pa_z[k] * hf_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gf0_11, gf0_13, gf0_14, gf1_11, \
                         gf1_13, gf1_14, hf_20, hf_21, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * gf0_13[k]
                  - f_8 * gf1_13[k]
                  + pa_y[k] * hf_21[k];

        t_32[k] = f_5 * gf0_11[k]
                  - f_6 * gf1_11[k]
                  + pa_z[k] * hf_20[k];

        t_33[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_y[k] * hf_22[k];

        t_34[k] = pa_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_35, pb_z, hd_23, ip0_2, ip1_2, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hd_23[k]
                  + f_1 * ip0_2[k]
                  - f_2 * ip1_2[k]
                  + pb_z[k] * id_35[k];
    }
}

auto
compute_prim_if_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_23 = buffer.data(hd + 23);

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

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_20 = buffer.data(ip1 + 20);
    const auto *ip1_32 = buffer.data(ip1 + 32);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_61 = buffer.data(id + 61);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_0, hf_0, hf_1, \
                         ip0_0, ip1_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pa_y[k] * hf_0[k];

        t_2[k] = pa_z[k] * hf_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gf0_0, gf0_4, gf0_6, gf1_0, gf1_4, gf1_6, \
                         hf_2, hf_4, hf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gf0_4[k]
                 - f_6 * gf1_4[k]
                 + pa_x[k] * hf_4[k];

        t_5[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];

        t_6[k] = f_5 * gf0_6[k]
                 - f_6 * gf1_6[k]
                 + pa_x[k] * hf_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, gf0_1, gf0_2, gf0_7, gf1_1, \
                         gf1_2, gf1_7, hf_3, hf_5, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gf0_1[k]
                 - f_8 * gf1_1[k]
                 + pa_y[k] * hf_3[k];

        t_8[k] = f_7 * gf0_7[k]
                 - f_8 * gf1_7[k]
                 + pa_x[k] * hf_8[k];

        t_9[k] = pa_y[k] * hf_5[k];

        t_10[k] = f_7 * gf0_2[k]
                  - f_8 * gf1_2[k]
                  + pa_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gf0_3, gf0_8, gf0_9, gf1_3, gf1_8, \
                         gf1_9, hf_7, hf_11, hf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_8[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = f_5 * gf0_3[k]
                  - f_6 * gf1_3[k]
                  + pa_y[k] * hf_7[k];

        t_13[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pa_x[k] * hf_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_y, gf0_5, gf0_11, gf0_12, gf1_5, \
                         gf1_11, gf1_12, hf_9, hf_10, hf_13, hf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pa_y[k] * hf_9[k];

        t_15[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_x[k] * hf_13[k];

        t_16[k] = f_3 * gf0_12[k]
                  - f_4 * gf1_12[k]
                  + pa_x[k] * hf_14[k];

        t_17[k] = pa_y[k] * hf_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, gf0_5, gf0_14, gf1_5, \
                         gf1_14, hf_10, hf_15, hf_16, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gf0_5[k]
                  - f_6 * gf1_5[k]
                  + pa_z[k] * hf_10[k];

        t_19[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_15[k];

        t_20[k] = pa_x[k] * hf_16[k];

        t_21[k] = pa_x[k] * hf_18[k];

        t_22[k] = pa_x[k] * hf_19[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pa_z, pb_y, hd_16, hf_16, hf_20, \
                         hf_21, hf_23, ip0_1, ip1_20, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * hf_20[k];

        t_24[k] = pa_x[k] * hf_21[k];

        t_25[k] = pa_x[k] * hf_23[k];

        t_26[k] = f_0 * hd_16[k]
                  + f_1 * ip0_1[k]
                  - f_2 * ip1_20[k]
                  + pb_y[k] * id_39[k];

        t_27[k] = pa_z[k] * hf_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, gf0_9, gf0_10, gf0_12, gf1_9, gf1_10, \
                         gf1_12, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pa_z[k] * hf_17[k];

        t_29[k] = f_5 * gf0_12[k]
                  - f_6 * gf1_12[k]
                  + pa_y[k] * hf_19[k];

        t_30[k] = f_7 * gf0_10[k]
                  - f_8 * gf1_10[k]
                  + pa_z[k] * hf_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gf0_11, gf0_13, gf0_14, gf1_11, \
                         gf1_13, gf1_14, hf_20, hf_21, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * gf0_13[k]
                  - f_8 * gf1_13[k]
                  + pa_y[k] * hf_21[k];

        t_32[k] = f_5 * gf0_11[k]
                  - f_6 * gf1_11[k]
                  + pa_z[k] * hf_20[k];

        t_33[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_y[k] * hf_22[k];

        t_34[k] = pa_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_35, pb_z, hd_23, ip0_2, ip1_32, id_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hd_23[k]
                  + f_1 * ip0_2[k]
                  - f_2 * ip1_32[k]
                  + pb_z[k] * id_61[k];
    }
}

auto
compute_prim_if_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_29 = buffer.data(hd + 29);

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

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_19 = buffer.data(ip0 + 19);
    const auto *ip0_20 = buffer.data(ip0 + 20);
    const auto *ip0_21 = buffer.data(ip0 + 21);
    const auto *ip0_30 = buffer.data(ip0 + 30);
    const auto *ip0_31 = buffer.data(ip0 + 31);
    const auto *ip0_32 = buffer.data(ip0 + 32);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_24 = buffer.data(ip1 + 24);
    const auto *ip1_25 = buffer.data(ip1 + 25);
    const auto *ip1_26 = buffer.data(ip1 + 26);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip0_2, ip1_0, \
                         ip1_1, ip1_2, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_2[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gf0_0, gf0_4, gf1_0, \
                         gf1_4, hf_0, hf_1, hf_2, hf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * hf_0[k];

        t_4[k] = pa_z[k] * hf_0[k];

        t_5[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];

        t_6[k] = f_5 * gf0_4[k]
                 - f_6 * gf1_4[k]
                 + pa_x[k] * hf_4[k];

        t_7[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, gf0_1, gf0_6, gf0_7, gf1_1, gf1_6, \
                         gf1_7, hf_3, hf_5, hf_6, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gf0_6[k]
                 - f_6 * gf1_6[k]
                 + pa_x[k] * hf_6[k];

        t_9[k] = f_7 * gf0_1[k]
                 - f_8 * gf1_1[k]
                 + pa_y[k] * hf_3[k];

        t_10[k] = f_7 * gf0_7[k]
                  - f_8 * gf1_7[k]
                  + pa_x[k] * hf_8[k];

        t_11[k] = pa_y[k] * hf_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pa_z, gf0_2, gf0_3, gf0_8, gf1_2, \
                         gf1_3, gf1_8, hf_5, hf_7, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gf0_2[k]
                  - f_8 * gf1_2[k]
                  + pa_z[k] * hf_5[k];

        t_13[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_8[k]
                  + pa_x[k] * hf_11[k];

        t_14[k] = f_5 * gf0_3[k]
                  - f_6 * gf1_3[k]
                  + pa_y[k] * hf_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, gf0_5, gf0_9, gf0_11, gf1_5, gf1_9, \
                         gf1_11, hf_9, hf_12, hf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pa_x[k] * hf_12[k];

        t_16[k] = f_3 * gf0_5[k]
                  - f_4 * gf1_5[k]
                  + pa_y[k] * hf_9[k];

        t_17[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_x[k] * hf_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pa_z, gf0_5, gf0_12, gf0_14, \
                         gf1_5, gf1_12, gf1_14, hf_10, hf_14, hf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gf0_12[k]
                  - f_4 * gf1_12[k]
                  + pa_x[k] * hf_14[k];

        t_19[k] = pa_y[k] * hf_10[k];

        t_20[k] = f_5 * gf0_5[k]
                  - f_6 * gf1_5[k]
                  + pa_z[k] * hf_10[k];

        t_21[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, hf_16, hf_18, hf_19, hf_20, \
                         hf_21, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * hf_16[k];

        t_23[k] = pa_x[k] * hf_18[k];

        t_24[k] = pa_x[k] * hf_19[k];

        t_25[k] = pa_x[k] * hf_20[k];

        t_26[k] = pa_x[k] * hf_21[k];

        t_27[k] = pa_x[k] * hf_23[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, pb_y, pb_z, hd_19, ip0_19, ip0_20, ip0_21, \
                         ip1_13, ip1_14, ip1_15, id_28, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * ip0_19[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_28[k];

        t_29[k] = f_0 * hd_19[k]
                  + f_1 * ip0_20[k]
                  - f_2 * ip1_14[k]
                  + pb_y[k] * id_29[k];

        t_30[k] = f_1 * ip0_21[k]
                  - f_2 * ip1_15[k]
                  + pb_z[k] * id_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gf0_9, gf0_10, gf0_12, gf1_9, \
                         gf1_10, gf1_12, hf_16, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_z[k] * hf_16[k];

        t_32[k] = f_3 * gf0_9[k]
                  - f_4 * gf1_9[k]
                  + pa_z[k] * hf_17[k];

        t_33[k] = f_5 * gf0_12[k]
                  - f_6 * gf1_12[k]
                  + pa_y[k] * hf_19[k];

        t_34[k] = f_7 * gf0_10[k]
                  - f_8 * gf1_10[k]
                  + pa_z[k] * hf_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, gf0_11, gf0_13, gf0_14, gf1_11, \
                         gf1_13, gf1_14, hf_20, hf_21, hf_22, hf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * gf0_13[k]
                  - f_8 * gf1_13[k]
                  + pa_y[k] * hf_21[k];

        t_36[k] = f_5 * gf0_11[k]
                  - f_6 * gf1_11[k]
                  + pa_z[k] * hf_20[k];

        t_37[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_y[k] * hf_22[k];

        t_38[k] = pa_y[k] * hf_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, pb_z, hd_29, ip0_30, ip0_31, ip0_32, \
                         ip1_24, ip1_25, ip1_26, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * ip0_30[k]
                  - f_2 * ip1_24[k]
                  + pb_x[k] * id_39[k];

        t_40[k] = f_1 * ip0_31[k]
                  - f_2 * ip1_25[k]
                  + pb_y[k] * id_40[k];

        t_41[k] = f_0 * hd_29[k]
                  + f_1 * ip0_32[k]
                  - f_2 * ip1_26[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_74 = buffer.data(id + 74);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_0, hf_0, hf_1, \
                         ip0_0, ip1_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pa_y[k] * hf_0[k];

        t_2[k] = pa_z[k] * hf_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, pb_x, gf0_0, gf0_5, gf1_0, gf1_5, hd_4, \
                         hf_2, hf_5, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hd_4[k]
                 + pb_x[k] * id_10[k];

        t_5[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_5[k]
                 + pa_x[k] * hf_5[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_x, gf0_1, gf0_8, gf1_1, gf1_8, hd_6, \
                         hf_3, hf_8, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_16[k];

        t_8[k] = f_6 * gf0_8[k]
                 - f_7 * gf1_8[k]
                 + pa_x[k] * hf_8[k];

        t_9[k] = f_8 * gf0_1[k]
                 - f_9 * gf1_1[k]
                 + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_x, gf0_2, gf0_9, gf1_2, \
                         gf1_10, hd_8, hf_6, hf_11, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_10 * hd_8[k]
                  + pb_x[k] * id_18[k];

        t_11[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = pa_y[k] * hf_6[k];

        t_13[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_x, gf0_3, gf0_10, gf1_3, gf1_12, \
                         hd_11, hf_9, hf_15, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_10 * hd_11[k]
                  + pb_x[k] * id_28[k];

        t_15[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_15[k];

        t_16[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_x, gf0_6, gf0_11, gf1_6, gf1_13, \
                         hd_12, hf_12, hf_17, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_11 * hd_12[k]
                  + pb_x[k] * id_30[k];

        t_18[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pa_x[k] * hf_17[k];

        t_19[k] = f_3 * gf0_6[k]
                  - f_4 * gf1_6[k]
                  + pa_y[k] * hf_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pa_z, gf0_6, gf0_13, gf0_15, \
                         gf1_6, gf1_15, gf1_17, hf_13, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_15[k]
                  + pa_x[k] * hf_18[k];

        t_21[k] = f_3 * gf0_15[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_19[k];

        t_22[k] = pa_y[k] * hf_13[k];

        t_23[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, gf0_17, gf1_20, hd_15, \
                         hf_21, hf_22, hf_24, hf_26, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_11 * hd_15[k]
                  + pb_x[k] * id_43[k];

        t_25[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_20[k]
                  + pa_x[k] * hf_21[k];

        t_26[k] = pa_x[k] * hf_22[k];

        t_27[k] = pa_x[k] * hf_24[k];

        t_28[k] = pa_x[k] * hf_26[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pa_z, pb_y, hd_16, hf_22, hf_27, \
                         hf_29, hf_32, ip0_1, ip1_10, id_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_x[k] * hf_27[k];

        t_30[k] = pa_x[k] * hf_29[k];

        t_31[k] = pa_x[k] * hf_32[k];

        t_32[k] = f_0 * hd_16[k]
                  + f_1 * ip0_1[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_56[k];

        t_33[k] = pa_z[k] * hf_22[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_y, gf0_11, gf0_15, gf1_13, gf1_17, \
                         hd_19, hf_23, hf_26, id_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_13[k]
                  + pa_z[k] * hf_23[k];

        t_35[k] = f_5 * hd_19[k]
                  + pb_y[k] * id_63[k];

        t_36[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_17[k]
                  + pa_y[k] * hf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_y, gf0_12, gf0_16, gf1_14, gf1_19, \
                         hd_21, hf_24, hf_29, id_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_14[k]
                  + pa_z[k] * hf_24[k];

        t_38[k] = f_10 * hd_21[k]
                  + pb_y[k] * id_66[k];

        t_39[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_19[k]
                  + pa_y[k] * hf_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pb_y, gf0_13, gf0_17, gf1_15, \
                         gf1_20, hd_22, hf_27, hf_31, hf_32, id_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_15[k]
                  + pa_z[k] * hf_27[k];

        t_41[k] = f_11 * hd_22[k]
                  + pb_y[k] * id_69[k];

        t_42[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_20[k]
                  + pa_y[k] * hf_31[k];

        t_43[k] = pa_y[k] * hf_32[k];
    }

#pragma omp simd aligned(t_44, pb_z, hd_23, ip0_2, ip1_17, id_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * hd_23[k]
                  + f_1 * ip0_2[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_74[k];
    }
}

auto
compute_prim_if_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);

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
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_z[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hd_0, hd_1, hf_0, hf_1, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_1[k];

        t_8[k] = pa_z[k] * hf_0[k];

        t_9[k] = f_4 * hd_0[k]
                 + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_2, hd_8, hf_2, \
                         hf_3, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_2[k];

        t_11[k] = f_5 * gf0_0[k]
                  - f_6 * gf1_0[k]
                  + pa_y[k] * hf_3[k];

        t_12[k] = f_7 * hd_8[k]
                  + pb_x[k] * id_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_z, pb_z, gf0_0, gf0_5, gf1_0, gf1_5, hf_4, \
                         hf_7, ip0_3, ip1_3, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * gf0_5[k]
                  - f_9 * gf1_5[k]
                  + pa_x[k] * hf_7[k];

        t_14[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_9[k];

        t_15[k] = f_5 * gf0_0[k]
                  - f_6 * gf1_0[k]
                  + pa_z[k] * hf_4[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_y, pb_z, hd_5, hd_12, ip0_4, ip1_4, id_10, \
                         id_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_10 * hd_5[k]
                  + pb_z[k] * id_10[k];

        t_17[k] = f_7 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_18[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pb_x, gf0_1, gf0_8, gf1_1, gf1_8, \
                         hd_14, hf_5, hf_10, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_8 * gf0_8[k]
                  - f_9 * gf1_8[k]
                  + pa_x[k] * hf_10[k];

        t_20[k] = f_11 * gf0_1[k]
                  - f_12 * gf1_1[k]
                  + pa_y[k] * hf_5[k];

        t_21[k] = f_3 * hd_14[k]
                  + pb_x[k] * id_14[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, pb_z, gf0_10, gf1_10, hf_8, hf_13, \
                         ip0_5, ip1_5, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_11 * gf0_10[k]
                  - f_12 * gf1_10[k]
                  + pa_x[k] * hf_13[k];

        t_23[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_15[k];

        t_24[k] = pa_y[k] * hf_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_x, pb_z, gf0_2, gf1_2, hd_10, hd_19, hf_8, \
                         id_17, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * gf0_2[k]
                  - f_12 * gf1_2[k]
                  + pa_z[k] * hf_8[k];

        t_26[k] = f_3 * hd_10[k]
                  + pb_z[k] * id_17[k];

        t_27[k] = f_3 * hd_19[k]
                  + pb_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_y, gf0_3, gf0_12, gf1_3, gf1_12, \
                         hf_11, hf_17, ip0_6, ip1_6, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_18[k];

        t_29[k] = f_11 * gf0_12[k]
                  - f_12 * gf1_12[k]
                  + pa_x[k] * hf_17[k];

        t_30[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_3[k]
                  + pa_y[k] * hf_11[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pb_x, pb_z, gf0_13, gf1_13, hd_21, hf_19, \
                         ip0_7, ip1_7, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * hd_21[k]
                  + pb_x[k] * id_21[k];

        t_32[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_13[k]
                  + pa_x[k] * hf_19[k];

        t_33[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_22[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pa_y, gf0_6, gf0_15, gf0_17, gf1_6, \
                         gf1_15, gf1_17, hf_14, hf_15, hf_20, hf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * gf0_6[k]
                  - f_6 * gf1_6[k]
                  + pa_y[k] * hf_14[k];

        t_35[k] = f_5 * gf0_15[k]
                  - f_6 * gf1_15[k]
                  + pa_x[k] * hf_20[k];

        t_36[k] = f_5 * gf0_17[k]
                  - f_6 * gf1_17[k]
                  + pa_x[k] * hf_21[k];

        t_37[k] = pa_y[k] * hf_15[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_z, pb_x, pb_z, gf0_6, gf1_6, hd_17, hd_25, \
                         hf_15, id_27, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * gf0_6[k]
                  - f_9 * gf1_6[k]
                  + pa_z[k] * hf_15[k];

        t_39[k] = f_7 * hd_17[k]
                  + pb_z[k] * id_27[k];

        t_40[k] = f_10 * hd_25[k]
                  + pb_x[k] * id_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_x, pb_y, gf0_20, gf1_20, hd_26, hf_23, hf_24, \
                         ip0_8, ip1_8, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_28[k];

        t_42[k] = f_5 * gf0_20[k]
                  - f_6 * gf1_20[k]
                  + pa_x[k] * hf_23[k];

        t_43[k] = f_3 * hd_26[k]
                  + pa_x[k] * hf_24[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, pa_x, pb_x, hd_27, hf_25, hf_28, \
                         hf_30, hf_31, hf_33, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * hd_27[k]
                  + pb_x[k] * id_31[k];

        t_45[k] = pa_x[k] * hf_25[k];

        t_46[k] = pa_x[k] * hf_28[k];

        t_47[k] = pa_x[k] * hf_30[k];

        t_48[k] = pa_x[k] * hf_31[k];

        t_49[k] = pa_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pb_x, pb_z, hd_24, hd_39, hd_41, hf_36, \
                         hf_38, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * hd_39[k]
                  + pa_x[k] * hf_36[k];

        t_51[k] = f_13 * hd_24[k]
                  + pb_z[k] * id_36[k];

        t_52[k] = f_4 * hd_41[k]
                  + pb_x[k] * id_37[k];

        t_53[k] = pa_x[k] * hf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pb_x, pb_y, pb_z, hd_27, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_38[k];

        t_55[k] = pb_z[k] * id_38[k];

        t_56[k] = pb_x[k] * id_39[k];

        t_57[k] = f_0 * hd_27[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_39[k];

        t_58[k] = pb_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, hd_27, hd_28, hf_25, \
                         ip0_11, ip1_11, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * hd_28[k]
                  + pb_y[k] * id_40[k];

        t_60[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_40[k];

        t_61[k] = pa_z[k] * hf_25[k];

        t_62[k] = f_4 * hd_27[k]
                  + pb_z[k] * id_41[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_z, pb_x, pb_y, hd_28, hd_30, hf_26, ip0_12, \
                         ip1_12, id_42, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * hd_30[k]
                  + pb_y[k] * id_42[k];

        t_64[k] = f_3 * hd_28[k]
                  + pa_z[k] * hf_26[k];

        t_65[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_43[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_z, pb_y, pb_z, gf0_13, gf1_13, hd_29, hd_33, \
                         hf_27, id_44, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * gf0_13[k]
                  - f_6 * gf1_13[k]
                  + pa_z[k] * hf_27[k];

        t_67[k] = f_10 * hd_29[k]
                  + pb_z[k] * id_44[k];

        t_68[k] = f_7 * hd_33[k]
                  + pb_y[k] * id_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pa_z, pb_x, gf0_14, gf0_17, gf1_14, gf1_17, \
                         hf_28, hf_30, ip0_13, ip1_13, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_8 * gf0_17[k]
                  - f_9 * gf1_17[k]
                  + pa_y[k] * hf_30[k];

        t_70[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_46[k];

        t_71[k] = f_11 * gf0_14[k]
                  - f_12 * gf1_14[k]
                  + pa_z[k] * hf_28[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pb_y, pb_z, gf0_19, gf1_19, hd_32, hd_36, \
                         hf_33, id_47, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * hd_32[k]
                  + pb_z[k] * id_47[k];

        t_73[k] = f_3 * hd_36[k]
                  + pb_y[k] * id_48[k];

        t_74[k] = f_11 * gf0_19[k]
                  - f_12 * gf1_19[k]
                  + pa_y[k] * hf_33[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_x, pb_z, gf0_15, gf1_15, hd_35, hf_31, \
                         ip0_14, ip1_14, id_49, id_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_49[k];

        t_76[k] = f_8 * gf0_15[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_31[k];

        t_77[k] = f_7 * hd_35[k]
                  + pb_z[k] * id_50[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pb_y, pb_z, gf0_20, gf1_20, hd_37, \
                         hd_38, hd_40, hf_35, hf_37, id_51, id_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * hd_38[k]
                  + pb_y[k] * id_51[k];

        t_79[k] = f_5 * gf0_20[k]
                  - f_6 * gf1_20[k]
                  + pa_y[k] * hf_35[k];

        t_80[k] = f_3 * hd_40[k]
                  + pa_y[k] * hf_37[k];

        t_81[k] = f_13 * hd_37[k]
                  + pb_z[k] * id_52[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pb_x, pb_y, pb_z, hd_39, hd_41, \
                         hf_38, ip0_15, ip1_15, id_53, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_4 * hd_41[k]
                  + pb_y[k] * id_53[k];

        t_83[k] = pa_y[k] * hf_38[k];

        t_84[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_54[k];

        t_85[k] = pb_y[k] * id_54[k];

        t_86[k] = f_0 * hd_39[k]
                  + pb_z[k] * id_54[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pb_x, pb_y, pb_z, hd_40, hd_41, ip0_16, \
                         ip0_17, ip1_16, ip1_17, id_55, id_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * id_56[k];

        t_88[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_55[k];

        t_89[k] = f_0 * hd_40[k]
                  + pb_z[k] * id_55[k];

        t_90[k] = pb_y[k] * id_56[k];

        t_91[k] = f_0 * hd_41[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_56[k];
    }
}

auto
compute_prim_if_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_29 = buffer.data(hd + 29);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_1, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_x, gf0_0, gf0_5, gf1_0, gf1_5, hd_6, \
                         hf_2, hf_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_6[k];

        t_10[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_5[k];

        t_11[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pb_x, gf0_1, gf0_8, gf1_1, gf1_8, hd_8, \
                         hf_3, hf_8, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_13[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_8[k];

        t_14[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_y, pa_z, pb_x, gf0_2, gf0_10, gf1_2, \
                         gf1_10, hd_10, hf_6, hf_11, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_16[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_17[k] = pa_y[k] * hf_6[k];

        t_18[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pb_x, gf0_3, gf0_12, gf1_3, gf1_12, \
                         hd_13, hf_9, hf_15, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_10 * hd_13[k]
                  + pb_x[k] * id_13[k];

        t_20[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_15[k];

        t_21[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, pb_x, gf0_6, gf0_13, gf1_6, gf1_13, \
                         hd_14, hf_12, hf_16, id_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_11 * hd_14[k]
                  + pb_x[k] * id_15[k];

        t_23[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_13[k]
                  + pa_x[k] * hf_16[k];

        t_24[k] = f_3 * gf0_6[k]
                  - f_4 * gf1_6[k]
                  + pa_y[k] * hf_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_y, pa_z, gf0_6, gf0_15, gf0_17, \
                         gf1_6, gf1_15, gf1_17, hf_13, hf_17, hf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gf0_15[k]
                  - f_4 * gf1_15[k]
                  + pa_x[k] * hf_17[k];

        t_26[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_18[k];

        t_27[k] = pa_y[k] * hf_13[k];

        t_28[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pb_x, gf0_20, gf1_20, hd_17, \
                         hf_19, hf_20, hf_22, hf_24, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_21[k];

        t_30[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_20[k]
                  + pa_x[k] * hf_19[k];

        t_31[k] = pa_x[k] * hf_20[k];

        t_32[k] = pa_x[k] * hf_22[k];

        t_33[k] = pa_x[k] * hf_24[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, hf_25, hf_27, hf_29, \
                         ip0_9, ip1_9, id_28, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * hf_25[k];

        t_35[k] = pa_x[k] * hf_27[k];

        t_36[k] = pa_x[k] * hf_29[k];

        t_37[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_28[k];

        t_38[k] = pb_x[k] * id_29[k];

        t_39[k] = pb_x[k] * id_30[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_y, pb_z, hd_19, hf_20, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_29, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hd_19[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_29[k];

        t_41[k] = pb_z[k] * id_29[k];

        t_42[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_30[k];

        t_43[k] = pa_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_y, gf0_13, gf0_17, gf1_13, gf1_17, \
                         hd_23, hf_21, hf_24, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * gf0_13[k]
                  - f_4 * gf1_13[k]
                  + pa_z[k] * hf_21[k];

        t_45[k] = f_5 * hd_23[k]
                  + pb_y[k] * id_33[k];

        t_46[k] = f_6 * gf0_17[k]
                  - f_7 * gf1_17[k]
                  + pa_y[k] * hf_24[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pa_z, pb_y, gf0_14, gf0_19, gf1_14, gf1_19, \
                         hd_25, hf_22, hf_27, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_14[k]
                  + pa_z[k] * hf_22[k];

        t_48[k] = f_10 * hd_25[k]
                  + pb_y[k] * id_35[k];

        t_49[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_19[k]
                  + pa_y[k] * hf_27[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pa_z, pb_y, gf0_15, gf0_20, gf1_15, \
                         gf1_20, hd_26, hf_25, hf_28, hf_29, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_15[k]
                  + pa_z[k] * hf_25[k];

        t_51[k] = f_11 * hd_26[k]
                  + pb_y[k] * id_37[k];

        t_52[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_20[k]
                  + pa_y[k] * hf_28[k];

        t_53[k] = pa_y[k] * hf_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_55[k] = pb_x[k] * id_40[k];

        t_56[k] = pb_x[k] * id_41[k];

        t_57[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_58[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_59, pb_z, hd_29, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * hd_29[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_9 = buffer.data(gf1 + 9);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_11 = buffer.data(gf1 + 11);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_13 = buffer.data(gf1 + 13);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_17 = buffer.data(gf1 + 17);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_0, hf_0, hf_1, \
                         ip0_0, ip1_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pa_y[k] * hf_0[k];

        t_2[k] = pa_z[k] * hf_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, pb_x, gf0_0, gf0_5, gf1_0, gf1_5, hd_4, \
                         hf_2, hf_5, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hd_4[k]
                 + pb_x[k] * id_4[k];

        t_5[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_5[k]
                 + pa_x[k] * hf_5[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_x, gf0_1, gf0_8, gf1_1, gf1_8, hd_6, \
                         hf_3, hf_8, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_6[k];

        t_8[k] = f_6 * gf0_8[k]
                 - f_7 * gf1_8[k]
                 + pa_x[k] * hf_8[k];

        t_9[k] = f_8 * gf0_1[k]
                 - f_9 * gf1_1[k]
                 + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_x, gf0_2, gf0_9, gf1_2, gf1_9, hd_8, \
                         hf_6, hf_11, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_10 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_11[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_9[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pb_x, gf0_3, gf0_10, gf1_3, gf1_10, \
                         hd_10, hf_9, hf_14, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_14[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_14[k];

        t_15[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_z, pb_x, gf0_6, gf0_11, gf1_6, gf1_11, \
                         hd_11, hf_12, hf_15, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_11 * hd_11[k]
                  + pb_x[k] * id_12[k];

        t_17[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_x[k] * hf_15[k];

        t_18[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_12[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, gf0_17, gf1_17, hd_12, hf_16, \
                         hf_17, hf_26, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_11 * hd_12[k]
                  + pb_x[k] * id_14[k];

        t_20[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_16[k];

        t_21[k] = pa_x[k] * hf_17[k];

        t_22[k] = pa_x[k] * hf_26[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_y, gf0_11, gf1_11, hd_13, hd_16, \
                         hf_17, hf_18, ip0_1, ip1_1, id_17, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * hd_13[k]
                  + f_1 * ip0_1[k]
                  - f_2 * ip1_1[k]
                  + pb_y[k] * id_17[k];

        t_24[k] = pa_z[k] * hf_17[k];

        t_25[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_11[k]
                  + pa_z[k] * hf_18[k];

        t_26[k] = f_5 * hd_16[k]
                  + pb_y[k] * id_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pa_z, pb_y, gf0_12, gf0_15, gf1_12, gf1_15, \
                         hd_18, hf_19, hf_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_15[k]
                  + pa_y[k] * hf_21[k];

        t_28[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_z[k] * hf_19[k];

        t_29[k] = f_10 * hd_18[k]
                  + pb_y[k] * id_22[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pa_z, pb_y, gf0_13, gf0_16, gf1_13, gf1_16, \
                         hd_19, hf_22, hf_24, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_16[k]
                  + pa_y[k] * hf_24[k];

        t_31[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_13[k]
                  + pa_z[k] * hf_22[k];

        t_32[k] = f_11 * hd_19[k]
                  + pb_y[k] * id_24[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_z, gf0_17, gf1_17, hd_20, hf_25, hf_26, \
                         ip0_2, ip1_2, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_y[k] * hf_25[k];

        t_34[k] = pa_y[k] * hf_26[k];

        t_35[k] = f_0 * hd_20[k]
                  + f_1 * ip0_2[k]
                  - f_2 * ip1_2[k]
                  + pb_z[k] * id_26[k];
    }
}

auto
compute_prim_if_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_0, hf_0, hf_1, \
                         ip0_0, ip1_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pa_y[k] * hf_0[k];

        t_2[k] = pa_z[k] * hf_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, pb_x, gf0_0, gf0_5, gf1_0, gf1_5, hd_4, \
                         hf_2, hf_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hd_4[k]
                 + pb_x[k] * id_6[k];

        t_5[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_5[k]
                 + pa_x[k] * hf_5[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_x, gf0_1, gf0_8, gf1_1, gf1_8, hd_6, \
                         hf_3, hf_8, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_10[k];

        t_8[k] = f_6 * gf0_8[k]
                 - f_7 * gf1_8[k]
                 + pa_x[k] * hf_8[k];

        t_9[k] = f_8 * gf0_1[k]
                 - f_9 * gf1_1[k]
                 + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_x, gf0_2, gf0_9, gf1_2, gf1_10, \
                         hd_8, hf_6, hf_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_10 * hd_8[k]
                  + pb_x[k] * id_12[k];

        t_11[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pb_x, gf0_3, gf0_10, gf1_3, gf1_12, \
                         hd_10, hf_9, hf_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_16[k];

        t_14[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_15[k];

        t_15[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_z, pb_x, gf0_6, gf0_11, gf1_6, gf1_14, \
                         hd_11, hf_13, hf_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_11 * hd_11[k]
                  + pb_x[k] * id_18[k];

        t_17[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_17[k];

        t_18[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, gf0_17, gf1_23, hd_12, hf_21, \
                         hf_22, hf_32, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_11 * hd_12[k]
                  + pb_x[k] * id_22[k];

        t_20[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_21[k];

        t_21[k] = pa_x[k] * hf_22[k];

        t_22[k] = pa_x[k] * hf_32[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_y, gf0_11, gf1_14, hd_13, hd_16, \
                         hf_22, hf_23, ip0_1, ip1_10, id_26, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * hd_13[k]
                  + f_1 * ip0_1[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_24[k] = pa_z[k] * hf_22[k];

        t_25[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_14[k]
                  + pa_z[k] * hf_23[k];

        t_26[k] = f_5 * hd_16[k]
                  + pb_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pa_z, pb_y, gf0_12, gf0_15, gf1_15, gf1_18, \
                         hd_18, hf_24, hf_26, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_18[k]
                  + pa_y[k] * hf_26[k];

        t_28[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_24[k];

        t_29[k] = f_10 * hd_18[k]
                  + pb_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pa_z, pb_y, gf0_13, gf0_16, gf1_16, gf1_20, \
                         hd_19, hf_27, hf_29, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_20[k]
                  + pa_y[k] * hf_29[k];

        t_31[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_27[k];

        t_32[k] = f_11 * hd_19[k]
                  + pb_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_z, gf0_17, gf1_23, hd_20, hf_31, hf_32, \
                         ip0_2, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_23[k]
                  + pa_y[k] * hf_31[k];

        t_34[k] = pa_y[k] * hf_32[k];

        t_35[k] = f_0 * hd_20[k]
                  + f_1 * ip0_2[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_31(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_59 = buffer.data(hf + 59);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_z[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, hd_1, hf_0, hf_3, hf_5, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * id_2[k];

        t_6[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_7[k] = pa_y[k] * hf_0[k];

        t_8[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_9[k] = pa_y[k] * hf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_y, gf0_0, gf1_0, hd_2, \
                         hf_0, hf_3, hf_5, hf_6, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * hf_0[k];

        t_11[k] = pa_z[k] * hf_3[k];

        t_12[k] = pb_y[k] * id_5[k];

        t_13[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_5[k];

        t_14[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_y[k] * hf_6[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, gf0_5, gf1_5, hd_6, hf_12, \
                         id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * id_6[k];

        t_16[k] = f_6 * hd_6[k]
                  + pb_x[k] * id_7[k];

        t_17[k] = f_7 * gf0_5[k]
                  - f_8 * gf1_5[k]
                  + pa_x[k] * hf_12[k];

        t_18[k] = pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_7, hf_8, \
                         hf_9, ip0_3, ip1_3, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_8[k];

        t_20[k] = pa_z[k] * hf_7[k];

        t_21[k] = pa_y[k] * hf_9[k];

        t_22[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, hd_10, ip0_4, ip1_4, id_9, id_10, \
                         id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_y[k] * id_9[k];

        t_24[k] = f_6 * hd_10[k]
                  + pb_x[k] * id_11[k];

        t_25[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_10[k];

        t_26[k] = pb_y[k] * id_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, gf0_1, gf0_8, gf1_1, gf1_8, \
                         hf_10, hf_17, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_8[k]
                  + pa_x[k] * hf_17[k];

        t_28[k] = f_9 * gf0_1[k]
                  - f_10 * gf1_1[k]
                  + pa_y[k] * hf_10[k];

        t_29[k] = pb_z[k] * id_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, gf0_10, gf1_10, hd_12, \
                         hf_20, ip0_5, ip1_5, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hd_12[k]
                  + pb_x[k] * id_13[k];

        t_31[k] = f_9 * gf0_10[k]
                  - f_10 * gf1_10[k]
                  + pa_x[k] * hf_20[k];

        t_32[k] = pb_z[k] * id_13[k];

        t_33[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, hd_7, hd_9, hf_10, \
                         hf_12, hf_13, hf_14, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * hf_10[k];

        t_35[k] = pa_z[k] * hf_12[k];

        t_36[k] = f_3 * hd_7[k]
                  + pa_z[k] * hf_13[k];

        t_37[k] = pa_y[k] * hf_14[k];

        t_38[k] = f_3 * hd_9[k]
                  + pa_y[k] * hf_16[k];

        t_39[k] = pa_y[k] * hf_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_16, hf_14, \
                         ip0_6, ip1_6, id_15, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * gf0_2[k]
                  - f_10 * gf1_2[k]
                  + pa_z[k] * hf_14[k];

        t_41[k] = pb_y[k] * id_15[k];

        t_42[k] = f_3 * hd_16[k]
                  + pb_x[k] * id_17[k];

        t_43[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_3, \
                         gf1_12, hf_18, hf_26, id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_y[k] * id_17[k];

        t_45[k] = f_9 * gf0_12[k]
                  - f_10 * gf1_12[k]
                  + pa_x[k] * hf_26[k];

        t_46[k] = f_7 * gf0_3[k]
                  - f_8 * gf1_3[k]
                  + pa_y[k] * hf_18[k];

        t_47[k] = pb_z[k] * id_18[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_x, pb_z, gf0_14, gf1_14, hd_17, \
                         hf_29, ip0_7, ip1_7, id_19, id_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_19[k];

        t_49[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_14[k]
                  + pa_x[k] * hf_29[k];

        t_50[k] = pb_z[k] * id_19[k];

        t_51[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_20[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pa_z, gf0_6, gf1_6, hd_13, hf_18, \
                         hf_20, hf_21, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_z[k] * hf_18[k];

        t_53[k] = pa_z[k] * hf_20[k];

        t_54[k] = f_3 * hd_13[k]
                  + pa_z[k] * hf_21[k];

        t_55[k] = f_4 * gf0_6[k]
                  - f_5 * gf1_6[k]
                  + pa_y[k] * hf_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pa_y, gf0_16, gf0_18, gf1_16, gf1_18, \
                         hd_15, hf_23, hf_25, hf_30, hf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_4 * gf0_16[k]
                  - f_5 * gf1_16[k]
                  + pa_x[k] * hf_30[k];

        t_57[k] = f_4 * gf0_18[k]
                  - f_5 * gf1_18[k]
                  + pa_x[k] * hf_31[k];

        t_58[k] = pa_y[k] * hf_23[k];

        t_59[k] = f_3 * hd_15[k]
                  + pa_y[k] * hf_25[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_18, \
                         hf_23, hf_26, id_21, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * hf_26[k];

        t_61[k] = f_7 * gf0_6[k]
                  - f_8 * gf1_6[k]
                  + pa_z[k] * hf_23[k];

        t_62[k] = pb_y[k] * id_21[k];

        t_63[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_y, gf0_23, gf1_23, hd_19, hf_34, \
                         hf_35, ip0_8, ip1_8, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_22[k];

        t_65[k] = pb_y[k] * id_23[k];

        t_66[k] = f_4 * gf0_23[k]
                  - f_5 * gf1_23[k]
                  + pa_x[k] * hf_34[k];

        t_67[k] = f_3 * hd_19[k]
                  + pa_x[k] * hf_35[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_x, pa_z, pb_x, hd_20, hd_23, hf_27, \
                         hf_38, hf_43, hf_44, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_24[k];

        t_69[k] = pa_x[k] * hf_38[k];

        t_70[k] = pa_z[k] * hf_27[k];

        t_71[k] = f_3 * hd_23[k]
                  + pa_x[k] * hf_43[k];

        t_72[k] = pa_x[k] * hf_44[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, hd_26, hd_30, hf_46, hf_47, \
                         hf_48, hf_50, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_x[k] * hf_46[k];

        t_74[k] = f_3 * hd_26[k]
                  + pa_x[k] * hf_47[k];

        t_75[k] = pa_x[k] * hf_48[k];

        t_76[k] = pa_x[k] * hf_50[k];

        t_77[k] = f_3 * hd_30[k]
                  + pa_x[k] * hf_54[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_x, pb_x, pb_z, hd_32, hf_59, ip0_9, \
                         ip1_9, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_25[k];

        t_79[k] = pa_x[k] * hf_59[k];

        t_80[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_26[k];

        t_81[k] = pb_z[k] * id_26[k];

        t_82[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_x[k] * id_28[k];

        t_84[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_27[k];

        t_85[k] = pb_z[k] * id_27[k];

        t_86[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_28[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_z, pb_x, hd_21, hf_35, hf_38, hf_40, \
                         ip0_12, ip1_12, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * hf_35[k];

        t_88[k] = pb_x[k] * id_30[k];

        t_89[k] = pa_z[k] * hf_38[k];

        t_90[k] = f_3 * hd_21[k]
                  + pa_z[k] * hf_40[k];

        t_91[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pb_x, pb_y, gf0_14, gf1_14, hd_25, \
                         hf_41, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * id_32[k];

        t_93[k] = pb_x[k] * id_33[k];

        t_94[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_14[k]
                  + pa_z[k] * hf_41[k];

        t_95[k] = f_6 * hd_25[k]
                  + pb_y[k] * id_33[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pb_x, gf0_18, gf1_18, hf_46, ip0_13, \
                         ip1_13, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_7 * gf0_18[k]
                  - f_8 * gf1_18[k]
                  + pa_y[k] * hf_46[k];

        t_97[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_34[k];

        t_98[k] = pb_x[k] * id_35[k];

        t_99[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_y, pa_z, pb_y, gf0_15, gf0_20, gf1_15, \
                         gf1_20, hd_28, hf_44, hf_50, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_9 * gf0_15[k]
                   - f_10 * gf1_15[k]
                   + pa_z[k] * hf_44[k];

        t_101[k] = f_3 * hd_28[k]
                   + pb_y[k] * id_36[k];

        t_102[k] = f_9 * gf0_20[k]
                   - f_10 * gf1_20[k]
                   + pa_y[k] * hf_50[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_z, pb_x, gf0_16, gf1_16, hf_48, \
                         ip0_14, ip1_14, id_37, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * ip0_14[k]
                   - f_2 * ip1_14[k]
                   + pb_x[k] * id_37[k];

        t_104[k] = pb_x[k] * id_38[k];

        t_105[k] = pb_x[k] * id_39[k];

        t_106[k] = f_7 * gf0_16[k]
                   - f_8 * gf1_16[k]
                   + pa_z[k] * hf_48[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pb_x, pb_y, gf0_23, gf1_23, hd_29, \
                         hd_31, hf_53, hf_57, id_39, id_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_11 * hd_29[k]
                   + pb_y[k] * id_39[k];

        t_108[k] = f_4 * gf0_23[k]
                   - f_5 * gf1_23[k]
                   + pa_y[k] * hf_53[k];

        t_109[k] = pb_x[k] * id_40[k];

        t_110[k] = f_3 * hd_31[k]
                   + pa_y[k] * hf_57[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pa_y, pb_x, pb_y, hd_32, hf_59, \
                         ip0_15, ip1_15, id_41, id_42, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * hd_32[k]
                   + pb_y[k] * id_41[k];

        t_112[k] = pa_y[k] * hf_59[k];

        t_113[k] = f_1 * ip0_15[k]
                   - f_2 * ip1_15[k]
                   + pb_x[k] * id_42[k];

        t_114[k] = pb_y[k] * id_42[k];

        t_115[k] = pb_x[k] * id_43[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_x, pb_y, pb_z, hd_32, ip0_16, ip0_17, \
                         ip1_16, ip1_17, id_43, id_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pb_x[k] * id_44[k];

        t_117[k] = f_1 * ip0_16[k]
                   - f_2 * ip1_16[k]
                   + pb_y[k] * id_43[k];

        t_118[k] = pb_y[k] * id_44[k];

        t_119[k] = f_0 * hd_32[k]
                   + f_1 * ip0_17[k]
                   - f_2 * ip1_17[k]
                   + pb_z[k] * id_44[k];
    }
}

auto
compute_prim_if_electron_repulsion_32(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_29 = buffer.data(gf1 + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_47 = buffer.data(hf + 47);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_z[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_y, pb_z, hd_1, hf_0, hf_3, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * id_2[k];

        t_6[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_7[k] = pa_y[k] * hf_0[k];

        t_8[k] = f_3 * hd_1[k]
                 + pa_y[k] * hf_3[k];

        t_9[k] = pa_z[k] * hf_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, pb_z, gf0_0, gf1_0, hd_2, \
                         hd_6, hf_5, hf_6, id_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * hd_2[k]
                  + pa_z[k] * hf_5[k];

        t_11[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_y[k] * hf_6[k];

        t_12[k] = pb_z[k] * id_5[k];

        t_13[k] = f_6 * hd_6[k]
                  + pb_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_z, gf0_5, gf1_7, hf_10, ip0_3, ip1_3, \
                         id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * gf0_5[k]
                  - f_8 * gf1_7[k]
                  + pa_x[k] * hf_10[k];

        t_15[k] = pb_z[k] * id_6[k];

        t_16[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_18[k] = pb_y[k] * id_8[k];

        t_19[k] = f_6 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_20[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_3, \
                         gf1_10, hf_8, hf_13, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_y[k] * id_10[k];

        t_22[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_10[k]
                  + pa_x[k] * hf_13[k];

        t_23[k] = f_9 * gf0_1[k]
                  - f_10 * gf1_3[k]
                  + pa_y[k] * hf_8[k];

        t_24[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, pb_z, gf0_10, gf1_12, hd_12, \
                         hf_16, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_26[k] = f_9 * gf0_10[k]
                  - f_10 * gf1_12[k]
                  + pa_x[k] * hf_16[k];

        t_27[k] = pb_z[k] * id_12[k];

        t_28[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_x, pb_y, gf0_2, gf1_4, hd_16, \
                         hf_11, id_14, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * hf_11[k];

        t_30[k] = f_9 * gf0_2[k]
                  - f_10 * gf1_4[k]
                  + pa_z[k] * hf_11[k];

        t_31[k] = pb_y[k] * id_14[k];

        t_32[k] = f_3 * hd_16[k]
                  + pb_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pb_y, gf0_12, gf1_14, hf_20, ip0_6, ip1_6, \
                         id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];

        t_34[k] = pb_y[k] * id_16[k];

        t_35[k] = f_9 * gf0_12[k]
                  - f_10 * gf1_14[k]
                  + pa_x[k] * hf_20[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pb_x, pb_z, gf0_3, gf1_5, hd_17, hf_14, \
                         id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_7 * gf0_3[k]
                  - f_8 * gf1_5[k]
                  + pa_y[k] * hf_14[k];

        t_37[k] = pb_z[k] * id_17[k];

        t_38[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pb_z, gf0_14, gf1_17, hf_22, ip0_7, ip1_7, \
                         id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_17[k]
                  + pa_x[k] * hf_22[k];

        t_40[k] = pb_z[k] * id_18[k];

        t_41[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, gf0_6, gf0_16, gf0_18, gf1_8, \
                         gf1_20, gf1_22, hf_17, hf_18, hf_23, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * gf0_6[k]
                  - f_5 * gf1_8[k]
                  + pa_y[k] * hf_17[k];

        t_43[k] = f_4 * gf0_16[k]
                  - f_5 * gf1_20[k]
                  + pa_x[k] * hf_23[k];

        t_44[k] = f_4 * gf0_18[k]
                  - f_5 * gf1_22[k]
                  + pa_x[k] * hf_24[k];

        t_45[k] = pa_y[k] * hf_18[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_x, pb_y, gf0_6, gf1_8, hd_18, hf_18, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_7 * gf0_6[k]
                  - f_8 * gf1_8[k]
                  + pa_z[k] * hf_18[k];

        t_47[k] = pb_y[k] * id_20[k];

        t_48[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_49[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pb_x, pb_y, gf0_23, gf1_29, hd_19, \
                         hd_20, hf_26, hf_27, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_y[k] * id_22[k];

        t_51[k] = f_4 * gf0_23[k]
                  - f_5 * gf1_29[k]
                  + pa_x[k] * hf_26[k];

        t_52[k] = f_3 * hd_19[k]
                  + pa_x[k] * hf_27[k];

        t_53[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pa_x, hd_30, hf_30, hf_34, hf_36, \
                         hf_37, hf_39, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * hf_30[k];

        t_55[k] = pa_x[k] * hf_34[k];

        t_56[k] = pa_x[k] * hf_36[k];

        t_57[k] = pa_x[k] * hf_37[k];

        t_58[k] = pa_x[k] * hf_39[k];

        t_59[k] = f_3 * hd_30[k]
                  + pa_x[k] * hf_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, hd_32, hf_47, ip0_9, \
                         ip1_9, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_61[k] = pa_x[k] * hf_47[k];

        t_62[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_63[k] = pb_z[k] * id_25[k];

        t_64[k] = pb_x[k] * id_26[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_x[k] * id_27[k];

        t_66[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_67[k] = pb_z[k] * id_26[k];

        t_68[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_z, pb_x, hd_21, hf_30, hf_32, \
                         ip0_12, ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_z[k] * hf_30[k];

        t_70[k] = f_3 * hd_21[k]
                  + pa_z[k] * hf_32[k];

        t_71[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_72[k] = pb_x[k] * id_30[k];

        t_73[k] = pb_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pa_z, pb_y, gf0_14, gf0_18, gf1_17, gf1_22, \
                         hd_25, hf_33, hf_36, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_17[k]
                  + pa_z[k] * hf_33[k];

        t_75[k] = f_6 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_76[k] = f_7 * gf0_18[k]
                  - f_8 * gf1_22[k]
                  + pa_y[k] * hf_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_z, pb_x, gf0_15, gf1_19, hf_34, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_78[k] = pb_x[k] * id_33[k];

        t_79[k] = pb_x[k] * id_34[k];

        t_80[k] = f_9 * gf0_15[k]
                  - f_10 * gf1_19[k]
                  + pa_z[k] * hf_34[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pb_x, pb_y, gf0_20, gf1_24, hd_28, \
                         hf_39, ip0_14, ip1_14, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_82[k] = f_9 * gf0_20[k]
                  - f_10 * gf1_24[k]
                  + pa_y[k] * hf_39[k];

        t_83[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_84[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pa_z, pb_x, pb_y, gf0_16, gf0_23, \
                         gf1_20, gf1_29, hd_29, hf_37, hf_41, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_x[k] * id_37[k];

        t_86[k] = f_7 * gf0_16[k]
                  - f_8 * gf1_20[k]
                  + pa_z[k] * hf_37[k];

        t_87[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_88[k] = f_4 * gf0_23[k]
                  - f_5 * gf1_29[k]
                  + pa_y[k] * hf_41[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_y, pb_x, pb_y, hd_31, hd_32, hf_45, \
                         hf_47, ip0_15, ip1_15, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_3 * hd_31[k]
                  + pa_y[k] * hf_45[k];

        t_90[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_91[k] = pa_y[k] * hf_47[k];

        t_92[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_93[k] = pb_y[k] * id_39[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, pb_z, hd_32, ip0_16, \
                         ip0_17, ip1_16, ip1_17, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pb_x[k] * id_40[k];

        t_95[k] = pb_x[k] * id_41[k];

        t_96[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_97[k] = pb_y[k] * id_41[k];

        t_98[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_33(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_29 = buffer.data(gf0 + 29);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_26 = buffer.data(hd + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_1, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_x, gf0_0, gf0_7, gf1_0, gf1_5, hd_6, \
                         hf_2, hf_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_6[k];

        t_10[k] = f_6 * gf0_7[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_5[k];

        t_11[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pb_x, gf0_3, gf0_10, gf1_1, gf1_8, \
                         hd_8, hf_3, hf_8, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_13[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_8[k];

        t_14[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_y, pa_z, pb_x, gf0_4, gf0_12, gf1_2, \
                         gf1_10, hd_10, hf_6, hf_11, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_16[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_17[k] = pa_y[k] * hf_6[k];

        t_18[k] = f_8 * gf0_4[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pb_x, gf0_5, gf0_14, gf1_3, gf1_12, \
                         hd_12, hf_9, hf_15, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_20[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_15[k];

        t_21[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, pb_x, gf0_8, gf0_17, gf1_6, gf1_14, \
                         hd_13, hf_12, hf_16, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_11 * hd_13[k]
                  + pb_x[k] * id_14[k];

        t_23[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_16[k];

        t_24[k] = f_3 * gf0_8[k]
                  - f_4 * gf1_6[k]
                  + pa_y[k] * hf_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_y, pa_z, gf0_8, gf0_20, gf0_22, \
                         gf1_6, gf1_16, gf1_18, hf_13, hf_17, hf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gf0_20[k]
                  - f_4 * gf1_16[k]
                  + pa_x[k] * hf_17[k];

        t_26[k] = f_3 * gf0_22[k]
                  - f_4 * gf1_18[k]
                  + pa_x[k] * hf_18[k];

        t_27[k] = pa_y[k] * hf_13[k];

        t_28[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pb_x, gf0_29, gf1_23, hd_14, \
                         hf_19, hf_20, hf_22, hf_24, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * hd_14[k]
                  + pb_x[k] * id_16[k];

        t_30[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_19[k];

        t_31[k] = pa_x[k] * hf_20[k];

        t_32[k] = pa_x[k] * hf_22[k];

        t_33[k] = pa_x[k] * hf_24[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, hf_25, hf_27, hf_29, \
                         ip0_9, ip1_9, id_19, id_20, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * hf_25[k];

        t_35[k] = pa_x[k] * hf_27[k];

        t_36[k] = pa_x[k] * hf_29[k];

        t_37[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_19[k];

        t_38[k] = pb_x[k] * id_20[k];

        t_39[k] = pb_x[k] * id_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_y, pb_z, hd_16, hf_20, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_20, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hd_16[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_20[k];

        t_41[k] = pb_z[k] * id_20[k];

        t_42[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_21[k];

        t_43[k] = pa_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_y, gf0_17, gf0_22, gf1_14, gf1_18, \
                         hd_20, hf_21, hf_24, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_14[k]
                  + pa_z[k] * hf_21[k];

        t_45[k] = f_5 * hd_20[k]
                  + pb_y[k] * id_24[k];

        t_46[k] = f_6 * gf0_22[k]
                  - f_7 * gf1_18[k]
                  + pa_y[k] * hf_24[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pa_z, pb_y, gf0_19, gf0_24, gf1_15, gf1_20, \
                         hd_22, hf_22, hf_27, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_22[k];

        t_48[k] = f_10 * hd_22[k]
                  + pb_y[k] * id_26[k];

        t_49[k] = f_8 * gf0_24[k]
                  - f_9 * gf1_20[k]
                  + pa_y[k] * hf_27[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pa_z, pb_y, gf0_20, gf0_29, gf1_16, \
                         gf1_23, hd_23, hf_25, hf_28, hf_29, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * gf0_20[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_25[k];

        t_51[k] = f_11 * hd_23[k]
                  + pb_y[k] * id_28[k];

        t_52[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_23[k]
                  + pa_y[k] * hf_28[k];

        t_53[k] = pa_y[k] * hf_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_30[k];

        t_55[k] = pb_x[k] * id_31[k];

        t_56[k] = pb_x[k] * id_32[k];

        t_57[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_31[k];

        t_58[k] = pb_y[k] * id_32[k];
    }

#pragma omp simd aligned(t_59, pb_z, hd_26, ip0_17, ip1_17, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * hd_26[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_32[k];
    }
}

auto
compute_prim_if_electron_repulsion_34(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_9 = buffer.data(gf0 + 9);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_11 = buffer.data(gf0 + 11);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_13 = buffer.data(gf0 + 13);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_17 = buffer.data(gf0 + 17);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gf0_0, gf1_0, hd_0, hf_0, hf_1, \
                         ip0_0, ip1_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pa_y[k] * hf_0[k];

        t_2[k] = pa_z[k] * hf_0[k];

        t_3[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, pb_x, gf0_0, gf0_5, gf1_0, gf1_5, hd_4, \
                         hf_2, hf_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hd_4[k]
                 + pb_x[k] * id_6[k];

        t_5[k] = f_6 * gf0_5[k]
                 - f_7 * gf1_5[k]
                 + pa_x[k] * hf_5[k];

        t_6[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_x, gf0_1, gf0_8, gf1_1, gf1_8, hd_6, \
                         hf_3, hf_8, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_10[k];

        t_8[k] = f_6 * gf0_8[k]
                 - f_7 * gf1_8[k]
                 + pa_x[k] * hf_8[k];

        t_9[k] = f_8 * gf0_1[k]
                 - f_9 * gf1_1[k]
                 + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_x, gf0_2, gf0_9, gf1_2, gf1_10, \
                         hd_8, hf_6, hf_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_10 * hd_8[k]
                  + pb_x[k] * id_12[k];

        t_11[k] = f_8 * gf0_9[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pb_x, gf0_3, gf0_10, gf1_3, gf1_12, \
                         hd_10, hf_9, hf_14, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_16[k];

        t_14[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_14[k];

        t_15[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_z, pb_x, gf0_6, gf0_11, gf1_6, gf1_14, \
                         hd_11, hf_12, hf_16, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_11 * hd_11[k]
                  + pb_x[k] * id_18[k];

        t_17[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_16[k];

        t_18[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_12[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, gf0_17, gf1_23, hd_12, hf_18, \
                         hf_19, hf_29, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_11 * hd_12[k]
                  + pb_x[k] * id_22[k];

        t_20[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_18[k];

        t_21[k] = pa_x[k] * hf_19[k];

        t_22[k] = pa_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_y, gf0_11, gf1_14, hd_13, hd_16, \
                         hf_19, hf_20, ip0_1, ip1_10, id_26, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * hd_13[k]
                  + f_1 * ip0_1[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_24[k] = pa_z[k] * hf_19[k];

        t_25[k] = f_3 * gf0_11[k]
                  - f_4 * gf1_14[k]
                  + pa_z[k] * hf_20[k];

        t_26[k] = f_5 * hd_16[k]
                  + pb_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pa_z, pb_y, gf0_12, gf0_15, gf1_15, gf1_18, \
                         hd_18, hf_21, hf_23, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gf0_15[k]
                  - f_7 * gf1_18[k]
                  + pa_y[k] * hf_23[k];

        t_28[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_21[k];

        t_29[k] = f_10 * hd_18[k]
                  + pb_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pa_z, pb_y, gf0_13, gf0_16, gf1_16, gf1_20, \
                         hd_19, hf_24, hf_26, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * gf0_16[k]
                  - f_9 * gf1_20[k]
                  + pa_y[k] * hf_26[k];

        t_31[k] = f_6 * gf0_13[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_24[k];

        t_32[k] = f_11 * hd_19[k]
                  + pb_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_z, gf0_17, gf1_23, hd_20, hf_28, hf_29, \
                         ip0_2, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_23[k]
                  + pa_y[k] * hf_28[k];

        t_34[k] = pa_y[k] * hf_29[k];

        t_35[k] = f_0 * hd_20[k]
                  + f_1 * ip0_2[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_35(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_53 = buffer.data(hf + 53);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_7, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_7[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_5, hd_6, hf_11, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_11[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_8, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_8[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_1, \
                         gf1_8, hf_9, hf_16, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_16[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_9[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_10, gf1_10, hd_12, \
                         hf_19, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_19[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_16, hf_13, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_13[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_3, \
                         gf1_12, hf_17, hf_24, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_24[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_17[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_14, gf1_14, hd_17, \
                         hf_26, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_26[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_18, hf_21, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_21[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_23, gf1_23, hd_20, \
                         hf_28, hf_33, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_28[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_53, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_53[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_33, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_33[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_14, gf1_14, hf_36, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_z[k] * hf_36[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_18, gf1_18, hd_25, \
                         hf_40, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_18[k]
                  - f_7 * gf1_18[k]
                  + pa_y[k] * hf_40[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_15, gf0_20, \
                         gf1_15, gf1_20, hd_28, hf_38, hf_44, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_15[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_38[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_20[k]
                  - f_9 * gf1_20[k]
                  + pa_y[k] * hf_44[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_16, gf1_16, hf_42, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_16[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_42[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_23, gf1_23, hd_29, hd_32, \
                         hf_46, hf_53, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_23[k]
                  + pa_y[k] * hf_46[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_53[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_36(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / alpha;
    const auto f_8 = 1.5 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_29 = buffer.data(gf1 + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_z[k] * id_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_y, pb_z, hd_2, hf_0, hf_5, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * id_2[k];

        t_6[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_7[k] = pa_y[k] * hf_0[k];

        t_8[k] = pa_z[k] * hf_0[k];

        t_9[k] = f_3 * hd_2[k]
                 + pa_z[k] * hf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_x, pb_z, gf0_0, gf1_0, hd_6, hf_6, id_5, \
                         id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_y[k] * hf_6[k];

        t_11[k] = pb_z[k] * id_5[k];

        t_12[k] = f_6 * hd_6[k]
                  + pb_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_z, gf0_5, gf1_7, hf_10, ip0_3, ip1_3, \
                         id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * gf0_5[k]
                  - f_8 * gf1_7[k]
                  + pa_x[k] * hf_10[k];

        t_14[k] = pb_z[k] * id_6[k];

        t_15[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * gf0_0[k]
                  - f_5 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_17[k] = pb_y[k] * id_8[k];

        t_18[k] = f_6 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_19[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_3, \
                         gf1_10, hf_8, hf_15, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_y[k] * id_10[k];

        t_21[k] = f_7 * gf0_8[k]
                  - f_8 * gf1_10[k]
                  + pa_x[k] * hf_15[k];

        t_22[k] = f_9 * gf0_1[k]
                  - f_10 * gf1_3[k]
                  + pa_y[k] * hf_8[k];

        t_23[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, gf0_10, gf1_12, hd_12, \
                         hf_18, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_25[k] = f_9 * gf0_10[k]
                  - f_10 * gf1_12[k]
                  + pa_x[k] * hf_18[k];

        t_26[k] = pb_z[k] * id_12[k];

        t_27[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_x, pb_y, gf0_2, gf1_4, hd_16, hf_12, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_9 * gf0_2[k]
                  - f_10 * gf1_4[k]
                  + pa_z[k] * hf_12[k];

        t_29[k] = pb_y[k] * id_14[k];

        t_30[k] = f_3 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_31[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_5, \
                         gf1_14, hf_16, hf_23, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_y[k] * id_16[k];

        t_33[k] = f_9 * gf0_12[k]
                  - f_10 * gf1_14[k]
                  + pa_x[k] * hf_23[k];

        t_34[k] = f_7 * gf0_3[k]
                  - f_8 * gf1_5[k]
                  + pa_y[k] * hf_16[k];

        t_35[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_x, pb_z, gf0_14, gf1_17, hd_17, \
                         hf_25, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_37[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_17[k]
                  + pa_x[k] * hf_25[k];

        t_38[k] = pb_z[k] * id_18[k];

        t_39[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, pb_y, gf0_6, gf1_8, hd_18, hf_20, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * gf0_6[k]
                  - f_8 * gf1_8[k]
                  + pa_z[k] * hf_20[k];

        t_41[k] = pb_y[k] * id_20[k];

        t_42[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_43[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_y, gf0_23, gf1_29, hd_20, \
                         hf_27, hf_31, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_y[k] * id_22[k];

        t_45[k] = f_4 * gf0_23[k]
                  - f_5 * gf1_29[k]
                  + pa_x[k] * hf_27[k];

        t_46[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_47[k] = pa_x[k] * hf_31[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pb_x, pb_z, hd_32, hf_50, ip0_9, \
                         ip1_9, id_24, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_49[k] = pa_x[k] * hf_50[k];

        t_50[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_51[k] = pb_z[k] * id_25[k];

        t_52[k] = pb_x[k] * id_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, pb_z, hd_20, ip0_10, ip0_11, \
                         ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_x[k] * id_27[k];

        t_54[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_55[k] = pb_z[k] * id_26[k];

        t_56[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_z, pb_x, hd_21, hf_31, hf_33, \
                         ip0_12, ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * hf_31[k];

        t_58[k] = f_3 * hd_21[k]
                  + pa_z[k] * hf_33[k];

        t_59[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_60[k] = pb_x[k] * id_30[k];

        t_61[k] = pb_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pa_z, pb_y, gf0_14, gf0_18, gf1_17, gf1_22, \
                         hd_25, hf_34, hf_38, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * gf0_14[k]
                  - f_5 * gf1_17[k]
                  + pa_z[k] * hf_34[k];

        t_63[k] = f_6 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_64[k] = f_7 * gf0_18[k]
                  - f_8 * gf1_22[k]
                  + pa_y[k] * hf_38[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_x, gf0_15, gf1_19, hf_36, ip0_13, \
                         ip1_13, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_66[k] = pb_x[k] * id_33[k];

        t_67[k] = pb_x[k] * id_34[k];

        t_68[k] = f_9 * gf0_15[k]
                  - f_10 * gf1_19[k]
                  + pa_z[k] * hf_36[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_x, pb_y, gf0_20, gf1_24, hd_28, \
                         hf_42, ip0_14, ip1_14, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_70[k] = f_9 * gf0_20[k]
                  - f_10 * gf1_24[k]
                  + pa_y[k] * hf_42[k];

        t_71[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_72[k] = pb_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pa_z, pb_x, pb_y, gf0_16, gf0_23, \
                         gf1_20, gf1_29, hd_29, hf_40, hf_44, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_x[k] * id_37[k];

        t_74[k] = f_7 * gf0_16[k]
                  - f_8 * gf1_20[k]
                  + pa_z[k] * hf_40[k];

        t_75[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_76[k] = f_4 * gf0_23[k]
                  - f_5 * gf1_29[k]
                  + pa_y[k] * hf_44[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_y, pb_x, pb_y, hd_31, hd_32, hf_48, \
                         hf_50, ip0_15, ip1_15, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_3 * hd_31[k]
                  + pa_y[k] * hf_48[k];

        t_78[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_79[k] = pa_y[k] * hf_50[k];

        t_80[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_81[k] = pb_y[k] * id_39[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pb_x, pb_y, pb_z, hd_32, ip0_16, \
                         ip0_17, ip1_16, ip1_17, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_x[k] * id_40[k];

        t_83[k] = pb_x[k] * id_41[k];

        t_84[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_85[k] = pb_y[k] * id_41[k];

        t_86[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_37(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_29 = buffer.data(gf0 + 29);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_29 = buffer.data(gf1 + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_44 = buffer.data(hf + 44);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_7, gf1_7, hd_6, hf_10, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_7[k]
                  - f_7 * gf1_7[k]
                  + pa_x[k] * hf_10[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, gf1_3, \
                         gf1_10, hf_8, hf_13, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_10[k]
                  + pa_x[k] * hf_13[k];

        t_20[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_3[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_12, gf1_12, hd_12, \
                         hf_16, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_16[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_4, gf1_4, hd_16, hf_11, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_4[k]
                  - f_9 * gf1_4[k]
                  + pa_z[k] * hf_11[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_5, gf0_14, gf1_5, \
                         gf1_14, hf_14, hf_19, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_14[k]
                  + pa_x[k] * hf_19[k];

        t_32[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_y[k] * hf_14[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_17, gf1_17, hd_17, \
                         hf_21, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_21[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_8, gf1_8, hd_18, hf_17, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_z[k] * hf_17[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_29, gf1_29, hd_20, \
                         hf_23, hf_27, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_29[k]
                  + pa_x[k] * hf_23[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_27[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_44, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_44[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_27, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_17, gf1_17, hf_30, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_z[k] * hf_30[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_22, gf1_22, hd_25, \
                         hf_33, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_22[k]
                  - f_7 * gf1_22[k]
                  + pa_y[k] * hf_33[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_19, gf0_24, \
                         gf1_19, gf1_24, hd_28, hf_31, hf_36, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_19[k]
                  + pa_z[k] * hf_31[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_24[k]
                  - f_9 * gf1_24[k]
                  + pa_y[k] * hf_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_20, gf1_20, hf_34, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_20[k]
                  - f_7 * gf1_20[k]
                  + pa_z[k] * hf_34[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_29, gf1_29, hd_29, hd_32, \
                         hf_38, hf_44, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_29[k]
                  + pa_y[k] * hf_38[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_38(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_29 = buffer.data(gf0 + 29);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_26 = buffer.data(hd + 26);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_1, \
                         ip0_2, ip1_2, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_x, gf0_0, gf0_7, gf1_0, gf1_5, hd_6, \
                         hf_2, hf_5, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * hd_6[k]
                 + pb_x[k] * id_6[k];

        t_10[k] = f_6 * gf0_7[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_5[k];

        t_11[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pb_x, gf0_3, gf0_10, gf1_1, gf1_8, \
                         hd_8, hf_3, hf_8, id_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * hd_8[k]
                  + pb_x[k] * id_8[k];

        t_13[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_8[k];

        t_14[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_3[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_z, pb_x, gf0_4, gf0_12, gf1_2, gf1_10, \
                         hd_10, hf_6, hf_11, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_10 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_16[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_11[k];

        t_17[k] = f_8 * gf0_4[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_y, pb_x, gf0_5, gf0_14, gf1_3, gf1_12, \
                         hd_12, hf_9, hf_14, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_19[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_14[k];

        t_20[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_z, pb_x, gf0_8, gf0_17, gf1_6, gf1_14, \
                         hd_13, hf_12, hf_15, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * hd_13[k]
                  + pb_x[k] * id_14[k];

        t_22[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_15[k];

        t_23[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, gf0_29, gf1_23, hd_14, hf_16, \
                         hf_17, hf_26, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_11 * hd_14[k]
                  + pb_x[k] * id_16[k];

        t_25[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_16[k];

        t_26[k] = pa_x[k] * hf_17[k];

        t_27[k] = pa_x[k] * hf_26[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pb_x, pb_y, pb_z, hd_16, ip0_9, ip0_10, \
                         ip1_9, ip1_10, id_19, id_20, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_19[k];

        t_29[k] = pb_x[k] * id_20[k];

        t_30[k] = pb_x[k] * id_21[k];

        t_31[k] = f_0 * hd_16[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_20[k];

        t_32[k] = pb_z[k] * id_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, gf0_17, gf1_14, hd_20, \
                         hf_17, hf_18, ip0_11, ip1_11, id_21, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_21[k];

        t_34[k] = pa_z[k] * hf_17[k];

        t_35[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_14[k]
                  + pa_z[k] * hf_18[k];

        t_36[k] = f_5 * hd_20[k]
                  + pb_y[k] * id_24[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_y, gf0_19, gf0_22, gf1_15, gf1_18, \
                         hd_22, hf_19, hf_21, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_6 * gf0_22[k]
                  - f_7 * gf1_18[k]
                  + pa_y[k] * hf_21[k];

        t_38[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_19[k];

        t_39[k] = f_10 * hd_22[k]
                  + pb_y[k] * id_26[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_y, gf0_20, gf0_24, gf1_16, gf1_20, \
                         hd_23, hf_22, hf_24, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_8 * gf0_24[k]
                  - f_9 * gf1_20[k]
                  + pa_y[k] * hf_24[k];

        t_41[k] = f_6 * gf0_20[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_22[k];

        t_42[k] = f_11 * hd_23[k]
                  + pb_y[k] * id_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pb_x, gf0_29, gf1_23, hf_25, \
                         hf_26, ip0_15, ip1_15, id_30, id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_23[k]
                  + pa_y[k] * hf_25[k];

        t_44[k] = pa_y[k] * hf_26[k];

        t_45[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_30[k];

        t_46[k] = pb_x[k] * id_31[k];

        t_47[k] = pb_x[k] * id_32[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pb_z, hd_26, ip0_16, ip0_17, ip1_16, ip1_17, \
                         id_31, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_31[k];

        t_49[k] = pb_y[k] * id_32[k];

        t_50[k] = f_0 * hd_26[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_32[k];
    }
}

auto
compute_prim_if_electron_repulsion_39(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_1 = buffer.data(gf1 + 1);
    const auto *gf1_2 = buffer.data(gf1 + 2);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_6 = buffer.data(gf1 + 6);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_15 = buffer.data(gf1 + 15);
    const auto *gf1_16 = buffer.data(gf1 + 16);
    const auto *gf1_18 = buffer.data(gf1 + 18);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_23 = buffer.data(gf1 + 23);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_50 = buffer.data(hf + 50);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_5, hd_6, hf_10, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_x[k] * hf_10[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_1, \
                         gf1_8, hf_8, hf_15, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_x[k] * hf_15[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_1[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_10, gf1_10, hd_12, \
                         hf_18, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_10[k]
                  + pa_x[k] * hf_18[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_2, hd_16, hf_12, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_2[k]
                  + pa_z[k] * hf_12[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_3, \
                         gf1_12, hf_16, hf_23, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_23[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_3[k]
                  + pa_y[k] * hf_16[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_14, gf1_14, hd_17, \
                         hf_25, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_x[k] * hf_25[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_6, hd_18, hf_20, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_6[k]
                  + pa_z[k] * hf_20[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_23, gf1_23, hd_20, \
                         hf_27, hf_31, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_23[k]
                  + pa_x[k] * hf_27[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_31[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_50, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_50[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_31, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_14, gf1_14, hf_34, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_14[k]
                  + pa_z[k] * hf_34[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_18, gf1_18, hd_25, \
                         hf_38, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_18[k]
                  - f_7 * gf1_18[k]
                  + pa_y[k] * hf_38[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_15, gf0_20, \
                         gf1_15, gf1_20, hd_28, hf_36, hf_42, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_15[k]
                  - f_9 * gf1_15[k]
                  + pa_z[k] * hf_36[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_20[k]
                  - f_9 * gf1_20[k]
                  + pa_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_16, gf1_16, hf_40, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_16[k]
                  - f_7 * gf1_16[k]
                  + pa_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_23, gf1_23, hd_29, hd_32, \
                         hf_44, hf_50, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_23[k]
                  + pa_y[k] * hf_44[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_40(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_1 = buffer.data(gf0 + 1);
    const auto *gf0_2 = buffer.data(gf0 + 2);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_6 = buffer.data(gf0 + 6);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_15 = buffer.data(gf0 + 15);
    const auto *gf0_16 = buffer.data(gf0 + 16);
    const auto *gf0_18 = buffer.data(gf0 + 18);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_23 = buffer.data(gf0 + 23);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_29 = buffer.data(gf1 + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_50 = buffer.data(hf + 50);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_5, gf1_7, hd_6, hf_10, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_7[k]
                  + pa_x[k] * hf_10[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_1, gf0_8, gf1_3, \
                         gf1_10, hf_8, hf_15, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_10[k]
                  + pa_x[k] * hf_15[k];

        t_20[k] = f_8 * gf0_1[k]
                  - f_9 * gf1_3[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_10, gf1_12, hd_12, \
                         hf_18, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_10[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_18[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_2, gf1_4, hd_16, hf_12, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_2[k]
                  - f_9 * gf1_4[k]
                  + pa_z[k] * hf_12[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_12, gf1_5, \
                         gf1_14, hf_16, hf_23, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_14[k]
                  + pa_x[k] * hf_23[k];

        t_32[k] = f_6 * gf0_3[k]
                  - f_7 * gf1_5[k]
                  + pa_y[k] * hf_16[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_14, gf1_17, hd_17, \
                         hf_25, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_25[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_6, gf1_8, hd_18, hf_20, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_6[k]
                  - f_7 * gf1_8[k]
                  + pa_z[k] * hf_20[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_23, gf1_29, hd_20, \
                         hf_27, hf_31, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_29[k]
                  + pa_x[k] * hf_27[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_31[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_50, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_50[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_31, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_14, gf1_17, hf_34, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_14[k]
                  - f_4 * gf1_17[k]
                  + pa_z[k] * hf_34[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_18, gf1_22, hd_25, \
                         hf_38, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_18[k]
                  - f_7 * gf1_22[k]
                  + pa_y[k] * hf_38[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_15, gf0_20, \
                         gf1_19, gf1_24, hd_28, hf_36, hf_42, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_15[k]
                  - f_9 * gf1_19[k]
                  + pa_z[k] * hf_36[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_20[k]
                  - f_9 * gf1_24[k]
                  + pa_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_16, gf1_20, hf_40, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_16[k]
                  - f_7 * gf1_20[k]
                  + pa_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_23, gf1_29, hd_29, hd_32, \
                         hf_44, hf_50, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_23[k]
                  - f_4 * gf1_29[k]
                  + pa_y[k] * hf_44[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

auto
compute_prim_if_electron_repulsion_41(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gf0, const size_t gf1,
                                      const size_t hd, const size_t hf, const size_t ip0,
                                      const size_t ip1, const size_t id, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / p;
    const auto f_11 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf0_0 = buffer.data(gf0 + 0);
    const auto *gf0_3 = buffer.data(gf0 + 3);
    const auto *gf0_4 = buffer.data(gf0 + 4);
    const auto *gf0_5 = buffer.data(gf0 + 5);
    const auto *gf0_7 = buffer.data(gf0 + 7);
    const auto *gf0_8 = buffer.data(gf0 + 8);
    const auto *gf0_10 = buffer.data(gf0 + 10);
    const auto *gf0_12 = buffer.data(gf0 + 12);
    const auto *gf0_14 = buffer.data(gf0 + 14);
    const auto *gf0_17 = buffer.data(gf0 + 17);
    const auto *gf0_19 = buffer.data(gf0 + 19);
    const auto *gf0_20 = buffer.data(gf0 + 20);
    const auto *gf0_22 = buffer.data(gf0 + 22);
    const auto *gf0_24 = buffer.data(gf0 + 24);
    const auto *gf0_29 = buffer.data(gf0 + 29);

    const auto *gf1_0 = buffer.data(gf1 + 0);
    const auto *gf1_3 = buffer.data(gf1 + 3);
    const auto *gf1_4 = buffer.data(gf1 + 4);
    const auto *gf1_5 = buffer.data(gf1 + 5);
    const auto *gf1_7 = buffer.data(gf1 + 7);
    const auto *gf1_8 = buffer.data(gf1 + 8);
    const auto *gf1_10 = buffer.data(gf1 + 10);
    const auto *gf1_12 = buffer.data(gf1 + 12);
    const auto *gf1_14 = buffer.data(gf1 + 14);
    const auto *gf1_17 = buffer.data(gf1 + 17);
    const auto *gf1_19 = buffer.data(gf1 + 19);
    const auto *gf1_20 = buffer.data(gf1 + 20);
    const auto *gf1_22 = buffer.data(gf1 + 22);
    const auto *gf1_24 = buffer.data(gf1 + 24);
    const auto *gf1_29 = buffer.data(gf1 + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_50 = buffer.data(hf + 50);

    const auto *ip0_0 = buffer.data(ip0 + 0);
    const auto *ip0_1 = buffer.data(ip0 + 1);
    const auto *ip0_2 = buffer.data(ip0 + 2);
    const auto *ip0_3 = buffer.data(ip0 + 3);
    const auto *ip0_4 = buffer.data(ip0 + 4);
    const auto *ip0_5 = buffer.data(ip0 + 5);
    const auto *ip0_6 = buffer.data(ip0 + 6);
    const auto *ip0_7 = buffer.data(ip0 + 7);
    const auto *ip0_8 = buffer.data(ip0 + 8);
    const auto *ip0_9 = buffer.data(ip0 + 9);
    const auto *ip0_10 = buffer.data(ip0 + 10);
    const auto *ip0_11 = buffer.data(ip0 + 11);
    const auto *ip0_12 = buffer.data(ip0 + 12);
    const auto *ip0_13 = buffer.data(ip0 + 13);
    const auto *ip0_14 = buffer.data(ip0 + 14);
    const auto *ip0_15 = buffer.data(ip0 + 15);
    const auto *ip0_16 = buffer.data(ip0 + 16);
    const auto *ip0_17 = buffer.data(ip0 + 17);

    const auto *ip1_0 = buffer.data(ip1 + 0);
    const auto *ip1_1 = buffer.data(ip1 + 1);
    const auto *ip1_2 = buffer.data(ip1 + 2);
    const auto *ip1_3 = buffer.data(ip1 + 3);
    const auto *ip1_4 = buffer.data(ip1 + 4);
    const auto *ip1_5 = buffer.data(ip1 + 5);
    const auto *ip1_6 = buffer.data(ip1 + 6);
    const auto *ip1_7 = buffer.data(ip1 + 7);
    const auto *ip1_8 = buffer.data(ip1 + 8);
    const auto *ip1_9 = buffer.data(ip1 + 9);
    const auto *ip1_10 = buffer.data(ip1 + 10);
    const auto *ip1_11 = buffer.data(ip1 + 11);
    const auto *ip1_12 = buffer.data(ip1 + 12);
    const auto *ip1_13 = buffer.data(ip1 + 13);
    const auto *ip1_14 = buffer.data(ip1 + 14);
    const auto *ip1_15 = buffer.data(ip1 + 15);
    const auto *ip1_16 = buffer.data(ip1 + 16);
    const auto *ip1_17 = buffer.data(ip1 + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, ip0_0, ip0_1, ip1_0, \
                         ip1_1, id_0, id_1, id_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip0_0[k]
                 - f_2 * ip1_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_1 * ip0_1[k]
                 - f_2 * ip1_1[k]
                 + pb_y[k] * id_1[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gf0_0, gf1_0, hf_0, hf_6, \
                         ip0_2, ip1_2, id_2, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * ip0_2[k]
                 - f_2 * ip1_2[k]
                 + pb_z[k] * id_2[k];

        t_6[k] = pa_y[k] * hf_0[k];

        t_7[k] = pa_z[k] * hf_0[k];

        t_8[k] = f_3 * gf0_0[k]
                 - f_4 * gf1_0[k]
                 + pa_y[k] * hf_6[k];

        t_9[k] = pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, gf0_7, gf1_7, hd_6, hf_10, \
                         ip0_3, ip1_3, id_6, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd_6[k]
                  + pb_x[k] * id_6[k];

        t_11[k] = f_6 * gf0_7[k]
                  - f_7 * gf1_7[k]
                  + pa_x[k] * hf_10[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_1 * ip0_3[k]
                  - f_2 * ip1_3[k]
                  + pb_z[k] * id_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gf0_0, gf1_0, hd_10, hf_7, \
                         ip0_4, ip1_4, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gf0_0[k]
                  - f_4 * gf1_0[k]
                  + pa_z[k] * hf_7[k];

        t_15[k] = pb_y[k] * id_8[k];

        t_16[k] = f_5 * hd_10[k]
                  + pb_x[k] * id_10[k];

        t_17[k] = f_1 * ip0_4[k]
                  - f_2 * ip1_4[k]
                  + pb_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, gf0_3, gf0_10, gf1_3, \
                         gf1_10, hf_8, hf_15, id_10, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * id_10[k];

        t_19[k] = f_6 * gf0_10[k]
                  - f_7 * gf1_10[k]
                  + pa_x[k] * hf_15[k];

        t_20[k] = f_8 * gf0_3[k]
                  - f_9 * gf1_3[k]
                  + pa_y[k] * hf_8[k];

        t_21[k] = pb_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, gf0_12, gf1_12, hd_12, \
                         hf_18, ip0_5, ip1_5, id_12, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hd_12[k]
                  + pb_x[k] * id_12[k];

        t_23[k] = f_8 * gf0_12[k]
                  - f_9 * gf1_12[k]
                  + pa_x[k] * hf_18[k];

        t_24[k] = pb_z[k] * id_12[k];

        t_25[k] = f_1 * ip0_5[k]
                  - f_2 * ip1_5[k]
                  + pb_z[k] * id_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, gf0_4, gf1_4, hd_16, hf_12, \
                         ip0_6, ip1_6, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gf0_4[k]
                  - f_9 * gf1_4[k]
                  + pa_z[k] * hf_12[k];

        t_27[k] = pb_y[k] * id_14[k];

        t_28[k] = f_10 * hd_16[k]
                  + pb_x[k] * id_16[k];

        t_29[k] = f_1 * ip0_6[k]
                  - f_2 * ip1_6[k]
                  + pb_y[k] * id_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, gf0_5, gf0_14, gf1_5, \
                         gf1_14, hf_16, hf_23, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * id_16[k];

        t_31[k] = f_8 * gf0_14[k]
                  - f_9 * gf1_14[k]
                  + pa_x[k] * hf_23[k];

        t_32[k] = f_6 * gf0_5[k]
                  - f_7 * gf1_5[k]
                  + pa_y[k] * hf_16[k];

        t_33[k] = pb_z[k] * id_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, gf0_17, gf1_17, hd_17, \
                         hf_25, ip0_7, ip1_7, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_11 * hd_17[k]
                  + pb_x[k] * id_18[k];

        t_35[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_x[k] * hf_25[k];

        t_36[k] = pb_z[k] * id_18[k];

        t_37[k] = f_1 * ip0_7[k]
                  - f_2 * ip1_7[k]
                  + pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, gf0_8, gf1_8, hd_18, hf_20, \
                         ip0_8, ip1_8, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * gf0_8[k]
                  - f_7 * gf1_8[k]
                  + pa_z[k] * hf_20[k];

        t_39[k] = pb_y[k] * id_20[k];

        t_40[k] = f_11 * hd_18[k]
                  + pb_x[k] * id_22[k];

        t_41[k] = f_1 * ip0_8[k]
                  - f_2 * ip1_8[k]
                  + pb_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_x, pb_y, gf0_29, gf1_29, hd_20, \
                         hf_27, hf_31, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * id_22[k];

        t_43[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_29[k]
                  + pa_x[k] * hf_27[k];

        t_44[k] = f_12 * hd_20[k]
                  + pb_x[k] * id_23[k];

        t_45[k] = pa_x[k] * hf_31[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, hd_32, hf_50, ip0_9, ip1_9, \
                         id_24, id_25, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * hd_32[k]
                  + pb_x[k] * id_24[k];

        t_47[k] = pa_x[k] * hf_50[k];

        t_48[k] = f_1 * ip0_9[k]
                  - f_2 * ip1_9[k]
                  + pb_x[k] * id_25[k];

        t_49[k] = pb_x[k] * id_26[k];

        t_50[k] = pb_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, hd_20, hf_31, ip0_10, \
                         ip0_11, ip1_10, ip1_11, id_26, id_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * hd_20[k]
                  + f_1 * ip0_10[k]
                  - f_2 * ip1_10[k]
                  + pb_y[k] * id_26[k];

        t_52[k] = pb_z[k] * id_26[k];

        t_53[k] = f_1 * ip0_11[k]
                  - f_2 * ip1_11[k]
                  + pb_z[k] * id_27[k];

        t_54[k] = pa_z[k] * hf_31[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_x, gf0_17, gf1_17, hf_34, ip0_12, \
                         ip1_12, id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ip0_12[k]
                  - f_2 * ip1_12[k]
                  + pb_x[k] * id_29[k];

        t_56[k] = pb_x[k] * id_30[k];

        t_57[k] = pb_x[k] * id_31[k];

        t_58[k] = f_3 * gf0_17[k]
                  - f_4 * gf1_17[k]
                  + pa_z[k] * hf_34[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, gf0_22, gf1_22, hd_25, \
                         hf_38, ip0_13, ip1_13, id_31, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_5 * hd_25[k]
                  + pb_y[k] * id_31[k];

        t_60[k] = f_6 * gf0_22[k]
                  - f_7 * gf1_22[k]
                  + pa_y[k] * hf_38[k];

        t_61[k] = f_1 * ip0_13[k]
                  - f_2 * ip1_13[k]
                  + pb_x[k] * id_32[k];

        t_62[k] = pb_x[k] * id_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pa_z, pb_x, pb_y, gf0_19, gf0_24, \
                         gf1_19, gf1_24, hd_28, hf_36, hf_42, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * id_34[k];

        t_64[k] = f_8 * gf0_19[k]
                  - f_9 * gf1_19[k]
                  + pa_z[k] * hf_36[k];

        t_65[k] = f_10 * hd_28[k]
                  + pb_y[k] * id_34[k];

        t_66[k] = f_8 * gf0_24[k]
                  - f_9 * gf1_24[k]
                  + pa_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, gf0_20, gf1_20, hf_40, ip0_14, \
                         ip1_14, id_35, id_36, id_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ip0_14[k]
                  - f_2 * ip1_14[k]
                  + pb_x[k] * id_35[k];

        t_68[k] = pb_x[k] * id_36[k];

        t_69[k] = pb_x[k] * id_37[k];

        t_70[k] = f_6 * gf0_20[k]
                  - f_7 * gf1_20[k]
                  + pa_z[k] * hf_40[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, gf0_29, gf1_29, hd_29, hd_32, \
                         hf_44, hf_50, id_37, id_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_11 * hd_29[k]
                  + pb_y[k] * id_37[k];

        t_72[k] = f_3 * gf0_29[k]
                  - f_4 * gf1_29[k]
                  + pa_y[k] * hf_44[k];

        t_73[k] = f_12 * hd_32[k]
                  + pb_y[k] * id_38[k];

        t_74[k] = pa_y[k] * hf_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pb_x, pb_y, ip0_15, ip0_16, ip1_15, \
                         ip1_16, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ip0_15[k]
                  - f_2 * ip1_15[k]
                  + pb_x[k] * id_39[k];

        t_76[k] = pb_x[k] * id_40[k];

        t_77[k] = pb_x[k] * id_41[k];

        t_78[k] = f_1 * ip0_16[k]
                  - f_2 * ip1_16[k]
                  + pb_y[k] * id_40[k];

        t_79[k] = pb_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_80, pb_z, hd_32, ip0_17, ip1_17, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hd_32[k]
                  + f_1 * ip0_17[k]
                  - f_2 * ip1_17[k]
                  + pb_z[k] * id_41[k];
    }
}

}  // namespace simdt2ceri
