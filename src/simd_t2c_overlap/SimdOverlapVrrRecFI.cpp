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


#include "SimdOverlapVrrRecFI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_fi_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pi, const size_t dh, const size_t di,
                          const size_t fg, const size_t fh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_12 = buffer.data(pi + 12);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dh_0, fg_0, fh_0, \
                         fh_1, fh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_0[k]
                 + f_1 * fg_0[k]
                 + pb_x[k] * fh_0[k];

        t_1[k] = pb_y[k] * fh_0[k];

        t_2[k] = pb_z[k] * fh_0[k];

        t_3[k] = f_2 * fg_0[k]
                 + pb_y[k] * fh_1[k];

        t_4[k] = pb_y[k] * fh_2[k];

        t_5[k] = f_2 * fg_0[k]
                 + pb_z[k] * fh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, fg_1, fg_2, fg_3, fh_3, \
                         fh_4, fh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fg_1[k]
                 + pb_y[k] * fh_3[k];

        t_7[k] = pb_z[k] * fh_3[k];

        t_8[k] = pb_y[k] * fh_4[k];

        t_9[k] = f_3 * fg_2[k]
                 + pb_z[k] * fh_4[k];

        t_10[k] = f_0 * fg_3[k]
                  + pb_y[k] * fh_5[k];

        t_11[k] = pb_z[k] * fh_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, dh_9, fg_4, fh_6, \
                         fh_7, fh_8, fh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * fg_4[k]
                  + pb_y[k] * fh_6[k];

        t_13[k] = pb_y[k] * fh_7[k];

        t_14[k] = f_0 * fg_4[k]
                  + pb_z[k] * fh_7[k];

        t_15[k] = f_0 * dh_9[k]
                  + pb_x[k] * fh_10[k];

        t_16[k] = pb_z[k] * fh_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, dh_10, dh_11, dh_12, fg_5, \
                         fh_9, fh_10, fh_11, fh_12, fh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dh_10[k]
                  + pb_x[k] * fh_11[k];

        t_18[k] = f_0 * dh_11[k]
                  + pb_x[k] * fh_12[k];

        t_19[k] = pb_y[k] * fh_9[k];

        t_20[k] = f_0 * dh_12[k]
                  + pb_x[k] * fh_14[k];

        t_21[k] = f_1 * fg_5[k]
                  + pb_y[k] * fh_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, fg_6, fg_7, fg_8, \
                         fh_10, fh_11, fh_12, fh_13, fh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * fh_10[k];

        t_23[k] = f_0 * fg_6[k]
                  + pb_y[k] * fh_11[k];

        t_24[k] = f_3 * fg_7[k]
                  + pb_y[k] * fh_12[k];

        t_25[k] = f_2 * fg_8[k]
                  + pb_y[k] * fh_13[k];

        t_26[k] = pb_y[k] * fh_14[k];

        t_27[k] = f_1 * fg_8[k]
                  + pb_z[k] * fh_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pa_y, pb_y, pb_z, dh_0, dh_1, \
                         di_0, di_1, di_2, fh_15, fh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * di_0[k];

        t_29[k] = f_2 * dh_0[k]
                  + pb_y[k] * fh_15[k];

        t_30[k] = pb_z[k] * fh_15[k];

        t_31[k] = f_3 * dh_1[k]
                  + pa_y[k] * di_1[k];

        t_32[k] = pb_z[k] * fh_16[k];

        t_33[k] = pa_y[k] * di_2[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, dh_3, dh_4, dh_5, \
                         di_3, di_5, di_6, fh_17, fh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dh_3[k]
                  + pa_y[k] * di_3[k];

        t_35[k] = pb_z[k] * fh_17[k];

        t_36[k] = f_2 * dh_4[k]
                  + pb_y[k] * fh_18[k];

        t_37[k] = pa_y[k] * di_5[k];

        t_38[k] = f_4 * dh_5[k]
                  + pa_y[k] * di_6[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_y, pb_z, dh_8, di_9, fg_10, fh_19, \
                         fh_20, fh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_z[k] * fh_19[k];

        t_40[k] = f_2 * fg_10[k]
                  + pb_z[k] * fh_20[k];

        t_41[k] = f_2 * dh_8[k]
                  + pb_y[k] * fh_21[k];

        t_42[k] = pa_y[k] * di_9[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_z, dh_18, dh_19, dh_20, dh_21, \
                         fh_22, fh_23, fh_25, fh_26, fh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * dh_18[k]
                  + pb_x[k] * fh_23[k];

        t_44[k] = pb_z[k] * fh_22[k];

        t_45[k] = f_3 * dh_19[k]
                  + pb_x[k] * fh_25[k];

        t_46[k] = f_3 * dh_20[k]
                  + pb_x[k] * fh_26[k];

        t_47[k] = f_3 * dh_21[k]
                  + pb_x[k] * fh_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pa_y, pb_z, pi_4, di_11, di_19, \
                         fg_11, fg_12, fh_23, fh_24, fh_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * di_11[k];

        t_49[k] = f_2 * pi_4[k]
                  + pa_x[k] * di_19[k];

        t_50[k] = pb_z[k] * fh_23[k];

        t_51[k] = f_2 * fg_11[k]
                  + pb_z[k] * fh_24[k];

        t_52[k] = f_3 * fg_12[k]
                  + pb_z[k] * fh_25[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pa_z, pb_y, pb_z, dh_12, di_0, \
                         di_13, fg_13, fh_26, fh_28, fh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * fg_13[k]
                  + pb_z[k] * fh_26[k];

        t_54[k] = f_2 * dh_12[k]
                  + pb_y[k] * fh_28[k];

        t_55[k] = pa_y[k] * di_13[k];

        t_56[k] = pa_z[k] * di_0[k];

        t_57[k] = pb_y[k] * fh_29[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, dh_0, dh_2, di_1, \
                         di_2, di_3, fh_29, fh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * dh_0[k]
                  + pb_z[k] * fh_29[k];

        t_59[k] = pa_z[k] * di_1[k];

        t_60[k] = pb_y[k] * fh_30[k];

        t_61[k] = f_3 * dh_2[k]
                  + pa_z[k] * di_2[k];

        t_62[k] = pa_z[k] * di_3[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_z, pb_y, dh_4, di_5, di_6, fg_16, \
                         fg_17, fh_31, fh_32, fh_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * fg_16[k]
                  + pb_y[k] * fh_31[k];

        t_64[k] = pb_y[k] * fh_32[k];

        t_65[k] = f_0 * dh_4[k]
                  + pa_z[k] * di_5[k];

        t_66[k] = pa_z[k] * di_6[k];

        t_67[k] = f_3 * fg_17[k]
                  + pb_y[k] * fh_33[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_x, pb_y, dh_8, dh_26, di_9, \
                         di_10, fg_18, fh_34, fh_35, fh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * fg_18[k]
                  + pb_y[k] * fh_34[k];

        t_69[k] = pb_y[k] * fh_35[k];

        t_70[k] = f_4 * dh_8[k]
                  + pa_z[k] * di_9[k];

        t_71[k] = pa_z[k] * di_10[k];

        t_72[k] = f_3 * dh_26[k]
                  + pb_x[k] * fh_37[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, dh_27, dh_28, dh_29, \
                         di_12, fh_36, fh_38, fh_39, fh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * dh_27[k]
                  + pb_x[k] * fh_38[k];

        t_74[k] = f_3 * dh_28[k]
                  + pb_x[k] * fh_39[k];

        t_75[k] = pb_y[k] * fh_36[k];

        t_76[k] = f_3 * dh_29[k]
                  + pb_x[k] * fh_41[k];

        t_77[k] = pa_z[k] * di_12[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_y, fg_19, fg_20, fg_21, fg_22, \
                         fh_37, fh_38, fh_39, fh_40, fh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * fg_19[k]
                  + pb_y[k] * fh_37[k];

        t_79[k] = f_0 * fg_20[k]
                  + pb_y[k] * fh_38[k];

        t_80[k] = f_3 * fg_21[k]
                  + pb_y[k] * fh_39[k];

        t_81[k] = f_2 * fg_22[k]
                  + pb_y[k] * fh_40[k];

        t_82[k] = pb_y[k] * fh_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_x, pb_y, pb_z, pi_12, dh_13, dh_30, \
                         dh_32, di_26, di_27, di_29, fh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * pi_12[k]
                  + pa_x[k] * di_26[k];

        t_84[k] = f_5 * dh_30[k]
                  + pa_x[k] * di_27[k];

        t_85[k] = f_3 * dh_13[k]
                  + pb_y[k] * fh_42[k];

        t_86[k] = pb_z[k] * fh_42[k];

        t_87[k] = f_4 * dh_32[k]
                  + pa_x[k] * di_29[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_x, pb_y, pb_z, dh_15, dh_35, di_32, \
                         fg_23, fh_43, fh_44, fh_45, fh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_z[k] * fh_43[k];

        t_89[k] = f_2 * fg_23[k]
                  + pb_z[k] * fh_44[k];

        t_90[k] = f_0 * dh_35[k]
                  + pa_x[k] * di_32[k];

        t_91[k] = pb_z[k] * fh_45[k];

        t_92[k] = f_3 * dh_15[k]
                  + pb_y[k] * fh_46[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_x, pb_z, dh_39, di_36, fg_24, fg_25, \
                         fh_46, fh_47, fh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_3 * fg_24[k]
                  + pb_z[k] * fh_46[k];

        t_94[k] = f_3 * dh_39[k]
                  + pa_x[k] * di_36[k];

        t_95[k] = pb_z[k] * fh_47[k];

        t_96[k] = f_2 * fg_25[k]
                  + pb_z[k] * fh_48[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pb_x, pb_y, pb_z, dh_17, dh_43, \
                         dh_45, fg_26, fh_49, fh_50, fh_51, fh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * dh_17[k]
                  + pb_y[k] * fh_49[k];

        t_98[k] = f_0 * fg_26[k]
                  + pb_z[k] * fh_49[k];

        t_99[k] = f_2 * dh_43[k]
                  + pb_x[k] * fh_51[k];

        t_100[k] = pb_z[k] * fh_50[k];

        t_101[k] = f_2 * dh_45[k]
                   + pb_x[k] * fh_52[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_x, pb_x, pb_z, dh_46, dh_47, \
                         dh_48, di_41, fh_51, fh_53, fh_54, fh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_2 * dh_46[k]
                   + pb_x[k] * fh_53[k];

        t_103[k] = f_2 * dh_47[k]
                   + pb_x[k] * fh_54[k];

        t_104[k] = f_2 * dh_48[k]
                   + pb_x[k] * fh_55[k];

        t_105[k] = pa_x[k] * di_41[k];

        t_106[k] = pb_z[k] * fh_51[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, pa_x, pa_y, di_20, di_42, \
                         di_43, di_44, di_45, di_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * di_42[k];

        t_108[k] = pa_x[k] * di_43[k];

        t_109[k] = pa_x[k] * di_44[k];

        t_110[k] = pa_x[k] * di_45[k];

        t_111[k] = pa_x[k] * di_46[k];

        t_112[k] = pa_y[k] * di_20[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_y, pa_z, pb_y, dh_23, \
                         di_14, di_15, di_16, di_21, di_22, fh_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * di_14[k];

        t_114[k] = pa_y[k] * di_21[k];

        t_115[k] = pa_z[k] * di_15[k];

        t_116[k] = f_2 * dh_23[k]
                   + pb_y[k] * fh_56[k];

        t_117[k] = pa_y[k] * di_22[k];

        t_118[k] = pa_z[k] * di_16[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_y, pa_z, pb_y, pb_z, dh_14, dh_24, \
                         di_17, di_23, fh_57, fh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * dh_14[k]
                   + pb_z[k] * fh_57[k];

        t_120[k] = f_2 * dh_24[k]
                   + pb_y[k] * fh_58[k];

        t_121[k] = pa_y[k] * di_23[k];

        t_122[k] = pa_z[k] * di_17[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pa_y, pb_y, pb_z, dh_16, dh_25, \
                         dh_49, di_24, di_47, fh_59, fh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_2 * dh_16[k]
                   + pb_z[k] * fh_59[k];

        t_124[k] = f_3 * dh_49[k]
                   + pa_x[k] * di_47[k];

        t_125[k] = f_2 * dh_25[k]
                   + pb_y[k] * fh_60[k];

        t_126[k] = pa_y[k] * di_24[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_z, pb_x, dh_51, dh_52, dh_53, \
                         dh_54, di_18, fh_61, fh_62, fh_63, fh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_z[k] * di_18[k];

        t_128[k] = f_2 * dh_51[k]
                   + pb_x[k] * fh_61[k];

        t_129[k] = f_2 * dh_52[k]
                   + pb_x[k] * fh_62[k];

        t_130[k] = f_2 * dh_53[k]
                   + pb_x[k] * fh_63[k];

        t_131[k] = f_2 * dh_54[k]
                   + pb_x[k] * fh_64[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, t_138, pa_x, pa_y, di_25, \
                         di_48, di_49, di_50, di_51, di_52, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * di_25[k];

        t_133[k] = pa_x[k] * di_48[k];

        t_134[k] = pa_x[k] * di_49[k];

        t_135[k] = pa_x[k] * di_50[k];

        t_136[k] = pa_x[k] * di_51[k];

        t_137[k] = pa_x[k] * di_52[k];

        t_138[k] = pa_x[k] * di_53[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_y, pb_z, dh_22, dh_56, \
                         di_54, di_55, fg_28, fh_65, fh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_x[k] * di_54[k];

        t_140[k] = f_5 * dh_56[k]
                   + pa_x[k] * di_55[k];

        t_141[k] = pb_y[k] * fh_65[k];

        t_142[k] = f_3 * dh_22[k]
                   + pb_z[k] * fh_65[k];

        t_143[k] = f_2 * fg_28[k]
                   + pb_y[k] * fh_66[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pa_x, pb_y, dh_61, di_60, fg_29, \
                         fg_30, fh_67, fh_68, fh_69, fh_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_y[k] * fh_67[k];

        t_145[k] = f_4 * dh_61[k]
                   + pa_x[k] * di_60[k];

        t_146[k] = f_3 * fg_29[k]
                   + pb_y[k] * fh_68[k];

        t_147[k] = f_2 * fg_30[k]
                   + pb_y[k] * fh_69[k];

        t_148[k] = pb_y[k] * fh_70[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pa_x, pb_y, dh_65, di_64, fg_31, \
                         fg_32, fg_33, fh_71, fh_72, fh_73, fh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * dh_65[k]
                   + pa_x[k] * di_64[k];

        t_150[k] = f_0 * fg_31[k]
                   + pb_y[k] * fh_71[k];

        t_151[k] = f_3 * fg_32[k]
                   + pb_y[k] * fh_72[k];

        t_152[k] = f_2 * fg_33[k]
                   + pb_y[k] * fh_73[k];

        t_153[k] = pb_y[k] * fh_74[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, dh_69, dh_70, dh_71, dh_72, \
                         di_69, fh_76, fh_77, fh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_3 * dh_69[k]
                   + pa_x[k] * di_69[k];

        t_155[k] = f_2 * dh_70[k]
                   + pb_x[k] * fh_76[k];

        t_156[k] = f_2 * dh_71[k]
                   + pb_x[k] * fh_77[k];

        t_157[k] = f_2 * dh_72[k]
                   + pb_x[k] * fh_78[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_x, pb_x, pb_y, dh_73, dh_75, \
                         di_70, di_71, fh_75, fh_79, fh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_2 * dh_73[k]
                   + pb_x[k] * fh_79[k];

        t_159[k] = pb_y[k] * fh_75[k];

        t_160[k] = f_2 * dh_75[k]
                   + pb_x[k] * fh_80[k];

        t_161[k] = pa_x[k] * di_70[k];

        t_162[k] = pa_x[k] * di_71[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pb_x, pb_y, di_72, \
                         di_73, di_74, di_75, fg_35, fh_80, fh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * di_72[k];

        t_164[k] = pa_x[k] * di_73[k];

        t_165[k] = pa_x[k] * di_74[k];

        t_166[k] = pb_y[k] * fh_80[k];

        t_167[k] = pa_x[k] * di_75[k];

        t_168[k] = f_1 * fg_35[k]
                   + pb_x[k] * fh_81[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pb_x, pb_z, fg_36, fg_37, fg_38, \
                         fh_81, fh_82, fh_83, fh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_4 * fg_36[k]
                   + pb_x[k] * fh_82[k];

        t_170[k] = pb_z[k] * fh_81[k];

        t_171[k] = f_0 * fg_37[k]
                   + pb_x[k] * fh_83[k];

        t_172[k] = pb_z[k] * fh_82[k];

        t_173[k] = f_0 * fg_38[k]
                   + pb_x[k] * fh_84[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pb_x, pb_z, fg_39, fg_40, fg_41, \
                         fg_42, fh_83, fh_85, fh_86, fh_87, fh_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_3 * fg_39[k]
                   + pb_x[k] * fh_85[k];

        t_175[k] = pb_z[k] * fh_83[k];

        t_176[k] = f_3 * fg_40[k]
                   + pb_x[k] * fh_86[k];

        t_177[k] = f_3 * fg_41[k]
                   + pb_x[k] * fh_87[k];

        t_178[k] = f_2 * fg_42[k]
                   + pb_x[k] * fh_88[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pb_x, pb_z, fg_44, fg_45, fg_46, \
                         fh_85, fh_89, fh_90, fh_91, fh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pb_z[k] * fh_85[k];

        t_180[k] = f_2 * fg_44[k]
                   + pb_x[k] * fh_89[k];

        t_181[k] = f_2 * fg_45[k]
                   + pb_x[k] * fh_90[k];

        t_182[k] = f_2 * fg_46[k]
                   + pb_x[k] * fh_91[k];

        t_183[k] = pb_x[k] * fh_92[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, t_189, pb_x, pb_y, dh_43, fg_42, \
                         fh_92, fh_93, fh_94, fh_95, fh_96, fh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_x[k] * fh_93[k];

        t_185[k] = pb_x[k] * fh_94[k];

        t_186[k] = pb_x[k] * fh_95[k];

        t_187[k] = pb_x[k] * fh_96[k];

        t_188[k] = pb_x[k] * fh_97[k];

        t_189[k] = f_0 * dh_43[k]
                   + f_1 * fg_42[k]
                   + pb_y[k] * fh_92[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pb_y, pb_z, dh_48, fg_42, fg_43, \
                         fg_44, fh_92, fh_93, fh_94, fh_95, fh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pb_z[k] * fh_92[k];

        t_191[k] = f_2 * fg_42[k]
                   + pb_z[k] * fh_93[k];

        t_192[k] = f_3 * fg_43[k]
                   + pb_z[k] * fh_94[k];

        t_193[k] = f_0 * fg_44[k]
                   + pb_z[k] * fh_95[k];

        t_194[k] = f_0 * dh_48[k]
                   + pb_y[k] * fh_97[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pa_z, pb_x, pb_z, di_27, di_28, \
                         di_29, fg_46, fg_47, fh_97, fh_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * fg_46[k]
                   + pb_z[k] * fh_97[k];

        t_196[k] = pa_z[k] * di_27[k];

        t_197[k] = pa_z[k] * di_28[k];

        t_198[k] = f_4 * fg_47[k]
                   + pb_x[k] * fh_98[k];

        t_199[k] = pa_z[k] * di_29[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pa_z, pb_x, di_32, fg_48, fg_49, \
                         fg_50, fg_51, fh_99, fh_100, fh_101, fh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_0 * fg_48[k]
                   + pb_x[k] * fh_99[k];

        t_201[k] = f_0 * fg_49[k]
                   + pb_x[k] * fh_100[k];

        t_202[k] = pa_z[k] * di_32[k];

        t_203[k] = f_3 * fg_50[k]
                   + pb_x[k] * fh_101[k];

        t_204[k] = f_3 * fg_51[k]
                   + pb_x[k] * fh_102[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, pa_z, pb_x, di_36, fg_52, fg_54, \
                         fg_55, fg_56, fh_103, fh_104, fh_105, fh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_3 * fg_52[k]
                   + pb_x[k] * fh_103[k];

        t_206[k] = pa_z[k] * di_36[k];

        t_207[k] = f_2 * fg_54[k]
                   + pb_x[k] * fh_104[k];

        t_208[k] = f_2 * fg_55[k]
                   + pb_x[k] * fh_105[k];

        t_209[k] = f_2 * fg_56[k]
                   + pb_x[k] * fh_106[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, pb_x, fg_57, fh_107, \
                         fh_108, fh_109, fh_110, fh_111, fh_112, \
                         fh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_2 * fg_57[k]
                   + pb_x[k] * fh_107[k];

        t_211[k] = pb_x[k] * fh_108[k];

        t_212[k] = pb_x[k] * fh_109[k];

        t_213[k] = pb_x[k] * fh_110[k];

        t_214[k] = pb_x[k] * fh_111[k];

        t_215[k] = pb_x[k] * fh_112[k];

        t_216[k] = pb_x[k] * fh_113[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_z, dh_43, dh_44, dh_45, \
                         dh_46, di_41, di_42, di_43, di_44, fh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_z[k] * di_41[k];

        t_218[k] = f_2 * dh_43[k]
                   + pb_z[k] * fh_108[k];

        t_219[k] = f_3 * dh_44[k]
                   + pa_z[k] * di_42[k];

        t_220[k] = f_0 * dh_45[k]
                   + pa_z[k] * di_43[k];

        t_221[k] = f_4 * dh_46[k]
                   + pa_z[k] * di_44[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pa_y, pb_x, pb_y, pi_12, dh_55, \
                         di_54, di_55, di_57, fg_58, fh_113, fh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_3 * dh_55[k]
                   + pb_y[k] * fh_113[k];

        t_223[k] = f_2 * pi_12[k]
                   + pa_y[k] * di_54[k];

        t_224[k] = pa_y[k] * di_55[k];

        t_225[k] = f_4 * fg_58[k]
                   + pb_x[k] * fh_114[k];

        t_226[k] = pa_y[k] * di_57[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, pa_y, pb_x, di_60, fg_59, fg_60, \
                         fg_61, fg_62, fh_115, fh_116, fh_117, fh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_0 * fg_59[k]
                   + pb_x[k] * fh_115[k];

        t_228[k] = f_0 * fg_60[k]
                   + pb_x[k] * fh_116[k];

        t_229[k] = pa_y[k] * di_60[k];

        t_230[k] = f_3 * fg_61[k]
                   + pb_x[k] * fh_117[k];

        t_231[k] = f_3 * fg_62[k]
                   + pb_x[k] * fh_118[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pa_y, pb_x, di_64, fg_63, fg_64, \
                         fg_65, fg_66, fh_119, fh_120, fh_121, fh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * fg_63[k]
                   + pb_x[k] * fh_119[k];

        t_233[k] = pa_y[k] * di_64[k];

        t_234[k] = f_2 * fg_64[k]
                   + pb_x[k] * fh_120[k];

        t_235[k] = f_2 * fg_65[k]
                   + pb_x[k] * fh_121[k];

        t_236[k] = f_2 * fg_66[k]
                   + pb_x[k] * fh_122[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, t_242, pa_y, pb_x, di_69, fg_67, \
                         fh_123, fh_124, fh_125, fh_126, fh_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_2 * fg_67[k]
                   + pb_x[k] * fh_123[k];

        t_238[k] = pa_y[k] * di_69[k];

        t_239[k] = pb_x[k] * fh_124[k];

        t_240[k] = pb_x[k] * fh_125[k];

        t_241[k] = pb_x[k] * fh_126[k];

        t_242[k] = pb_x[k] * fh_127[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, pa_y, pb_x, pb_z, dh_50, dh_70, \
                         dh_72, di_70, di_72, fh_124, fh_128, fh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pb_x[k] * fh_128[k];

        t_244[k] = pb_x[k] * fh_129[k];

        t_245[k] = f_5 * dh_70[k]
                   + pa_y[k] * di_70[k];

        t_246[k] = f_3 * dh_50[k]
                   + pb_z[k] * fh_124[k];

        t_247[k] = f_4 * dh_72[k]
                   + pa_y[k] * di_72[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_y, pb_y, dh_73, dh_74, dh_75, di_73, \
                         di_74, di_75, fh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_0 * dh_73[k]
                   + pa_y[k] * di_73[k];

        t_249[k] = f_3 * dh_74[k]
                   + pa_y[k] * di_74[k];

        t_250[k] = f_2 * dh_75[k]
                   + pb_y[k] * fh_129[k];

        t_251[k] = pa_y[k] * di_75[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, t_257, pb_x, pb_y, fg_69, fg_70, \
                         fg_71, fg_72, fh_130, fh_131, fh_132, fh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * fg_69[k]
                   + pb_x[k] * fh_130[k];

        t_253[k] = pb_y[k] * fh_130[k];

        t_254[k] = f_4 * fg_70[k]
                   + pb_x[k] * fh_131[k];

        t_255[k] = f_0 * fg_71[k]
                   + pb_x[k] * fh_132[k];

        t_256[k] = pb_y[k] * fh_131[k];

        t_257[k] = f_0 * fg_72[k]
                   + pb_x[k] * fh_133[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pb_x, pb_y, fg_73, fg_74, fg_75, \
                         fg_76, fh_133, fh_134, fh_135, fh_136, \
                         fh_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_3 * fg_73[k]
                   + pb_x[k] * fh_134[k];

        t_259[k] = f_3 * fg_74[k]
                   + pb_x[k] * fh_135[k];

        t_260[k] = pb_y[k] * fh_133[k];

        t_261[k] = f_3 * fg_75[k]
                   + pb_x[k] * fh_136[k];

        t_262[k] = f_2 * fg_76[k]
                   + pb_x[k] * fh_137[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, pb_x, pb_y, fg_77, fg_78, fg_80, \
                         fh_136, fh_138, fh_139, fh_140, fh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_2 * fg_77[k]
                   + pb_x[k] * fh_138[k];

        t_264[k] = f_2 * fg_78[k]
                   + pb_x[k] * fh_139[k];

        t_265[k] = pb_y[k] * fh_136[k];

        t_266[k] = f_2 * fg_80[k]
                   + pb_x[k] * fh_140[k];

        t_267[k] = pb_x[k] * fh_141[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, t_273, pb_x, pb_y, fg_76, fh_141, \
                         fh_142, fh_143, fh_144, fh_145, fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = pb_x[k] * fh_142[k];

        t_269[k] = pb_x[k] * fh_143[k];

        t_270[k] = pb_x[k] * fh_144[k];

        t_271[k] = pb_x[k] * fh_145[k];

        t_272[k] = pb_x[k] * fh_146[k];

        t_273[k] = f_1 * fg_76[k]
                   + pb_y[k] * fh_141[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_y, fg_77, fg_78, fg_79, fg_80, \
                         fh_142, fh_143, fh_144, fh_145, fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_4 * fg_77[k]
                   + pb_y[k] * fh_142[k];

        t_275[k] = f_0 * fg_78[k]
                   + pb_y[k] * fh_143[k];

        t_276[k] = f_3 * fg_79[k]
                   + pb_y[k] * fh_144[k];

        t_277[k] = f_2 * fg_80[k]
                   + pb_y[k] * fh_145[k];

        t_278[k] = pb_y[k] * fh_146[k];
    }

#pragma omp simd aligned(t_279, pb_z, dh_75, fg_80, fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_0 * dh_75[k]
                   + f_1 * fg_80[k]
                   + pb_z[k] * fh_146[k];
    }
}

auto
compute_prim_fi_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pi, const size_t dh, const size_t di,
                          const size_t fg, const size_t fh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_4 = buffer.data(pi + 4);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_15 = buffer.data(fg + 15);
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
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
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
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dh_0, fg_0, fg_1, fh_0, \
                         fh_1, fh_2, fh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_0[k]
                 + f_1 * fg_0[k]
                 + pb_x[k] * fh_0[k];

        t_1[k] = f_2 * fg_0[k]
                 + pb_y[k] * fh_1[k];

        t_2[k] = f_2 * fg_0[k]
                 + pb_z[k] * fh_2[k];

        t_3[k] = f_3 * fg_1[k]
                 + pb_y[k] * fh_3[k];

        t_4[k] = pb_z[k] * fh_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fg_2, fg_3, fg_4, fh_4, fh_5, \
                         fh_7, fh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * fg_2[k]
                 + pb_z[k] * fh_4[k];

        t_6[k] = f_0 * fg_3[k]
                 + pb_y[k] * fh_5[k];

        t_7[k] = pb_z[k] * fh_5[k];

        t_8[k] = f_2 * fg_4[k]
                 + pb_y[k] * fh_7[k];

        t_9[k] = f_0 * fg_4[k]
                 + pb_z[k] * fh_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, dh_9, dh_10, \
                         di_0, fg_5, fg_8, fh_9, fh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * dh_9[k]
                  + pb_x[k] * fh_9[k];

        t_11[k] = f_0 * dh_10[k]
                  + pb_x[k] * fh_12[k];

        t_12[k] = f_1 * fg_5[k]
                  + pb_y[k] * fh_9[k];

        t_13[k] = f_1 * fg_8[k]
                  + pb_z[k] * fh_12[k];

        t_14[k] = pa_y[k] * di_0[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pb_y, dh_0, dh_1, dh_3, dh_5, di_1, \
                         di_3, di_5, fh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * dh_0[k]
                  + pb_y[k] * fh_13[k];

        t_16[k] = f_3 * dh_1[k]
                  + pa_y[k] * di_1[k];

        t_17[k] = f_0 * dh_3[k]
                  + pa_y[k] * di_3[k];

        t_18[k] = f_4 * dh_5[k]
                  + pa_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pb_x, pb_z, pi_0, dh_0, dh_14, \
                         di_0, di_8, fh_18, fh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * dh_14[k]
                  + pb_x[k] * fh_18[k];

        t_20[k] = f_2 * pi_0[k]
                  + pa_x[k] * di_8[k];

        t_21[k] = pa_z[k] * di_0[k];

        t_22[k] = f_2 * dh_0[k]
                  + pb_z[k] * fh_23[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_y, dh_2, dh_4, dh_8, di_2, di_4, \
                         di_7, fg_15, fh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * dh_2[k]
                  + pa_z[k] * di_2[k];

        t_24[k] = f_0 * dh_4[k]
                  + pa_z[k] * di_4[k];

        t_25[k] = f_2 * fg_15[k]
                  + pb_y[k] * fh_26[k];

        t_26[k] = f_4 * dh_8[k]
                  + pa_z[k] * di_7[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_x, pb_y, pi_4, dh_11, dh_19, dh_20, \
                         di_9, di_10, fh_31, fh_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * dh_19[k]
                  + pb_x[k] * fh_31[k];

        t_28[k] = f_2 * pi_4[k]
                  + pa_x[k] * di_9[k];

        t_29[k] = f_5 * dh_20[k]
                  + pa_x[k] * di_10[k];

        t_30[k] = f_3 * dh_11[k]
                  + pb_y[k] * fh_32[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, dh_21, dh_23, dh_26, dh_28, \
                         di_11, di_12, di_13, di_15, fh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_4 * dh_21[k]
                  + pa_x[k] * di_11[k];

        t_32[k] = f_0 * dh_23[k]
                  + pa_x[k] * di_12[k];

        t_33[k] = f_3 * dh_26[k]
                  + pa_x[k] * di_13[k];

        t_34[k] = f_2 * dh_28[k]
                  + pb_x[k] * fh_38[k];

        t_35[k] = pa_x[k] * di_15[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_x, pb_z, dh_15, dh_40, di_19, di_20, \
                         di_21, di_23, fh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * di_19[k];

        t_37[k] = pa_x[k] * di_20[k];

        t_38[k] = pa_x[k] * di_21[k];

        t_39[k] = f_5 * dh_40[k]
                  + pa_x[k] * di_23[k];

        t_40[k] = f_3 * dh_15[k]
                  + pb_z[k] * fh_52[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_x, pb_x, dh_43, dh_46, dh_49, dh_55, \
                         di_25, di_27, di_30, di_35, fh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * dh_43[k]
                  + pa_x[k] * di_25[k];

        t_42[k] = f_0 * dh_46[k]
                  + pa_x[k] * di_27[k];

        t_43[k] = f_3 * dh_49[k]
                  + pa_x[k] * di_30[k];

        t_44[k] = f_2 * dh_55[k]
                  + pb_x[k] * fh_61[k];

        t_45[k] = pa_x[k] * di_35[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_z, fg_30, fg_31, fg_32, fg_33, \
                         fh_62, fh_63, fh_64, fh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * fg_30[k]
                  + pb_x[k] * fh_62[k];

        t_47[k] = f_4 * fg_31[k]
                  + pb_x[k] * fh_63[k];

        t_48[k] = f_0 * fg_32[k]
                  + pb_x[k] * fh_64[k];

        t_49[k] = pb_z[k] * fh_63[k];

        t_50[k] = f_0 * fg_33[k]
                  + pb_x[k] * fh_66[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pb_x, pb_z, fg_34, fg_35, fg_36, fg_37, \
                         fh_64, fh_67, fh_69, fh_70, fh_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fg_34[k]
                  + pb_x[k] * fh_67[k];

        t_52[k] = pb_z[k] * fh_64[k];

        t_53[k] = f_3 * fg_35[k]
                  + pb_x[k] * fh_69[k];

        t_54[k] = f_3 * fg_36[k]
                  + pb_x[k] * fh_70[k];

        t_55[k] = f_2 * fg_37[k]
                  + pb_x[k] * fh_71[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_z, fg_39, fg_40, fg_41, fh_67, \
                         fh_72, fh_73, fh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * fh_67[k];

        t_57[k] = f_2 * fg_39[k]
                  + pb_x[k] * fh_72[k];

        t_58[k] = f_2 * fg_40[k]
                  + pb_x[k] * fh_73[k];

        t_59[k] = f_2 * fg_41[k]
                  + pb_x[k] * fh_74[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_y, pb_z, dh_28, fg_37, fg_38, fg_39, \
                         fh_75, fh_76, fh_77, fh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * dh_28[k]
                  + f_1 * fg_37[k]
                  + pb_y[k] * fh_75[k];

        t_61[k] = f_2 * fg_37[k]
                  + pb_z[k] * fh_76[k];

        t_62[k] = f_3 * fg_38[k]
                  + pb_z[k] * fh_77[k];

        t_63[k] = f_0 * fg_39[k]
                  + pb_z[k] * fh_78[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_x, pb_y, pb_z, dh_33, fg_41, fg_42, fg_43, \
                         fh_80, fh_81, fh_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * dh_33[k]
                  + pb_y[k] * fh_80[k];

        t_65[k] = f_1 * fg_41[k]
                  + pb_z[k] * fh_80[k];

        t_66[k] = f_0 * fg_42[k]
                  + pb_x[k] * fh_81[k];

        t_67[k] = f_3 * fg_43[k]
                  + pb_x[k] * fh_82[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_x, pb_z, dh_28, di_15, fg_45, fg_46, \
                         fh_83, fh_84, fh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * fg_45[k]
                  + pb_x[k] * fh_83[k];

        t_69[k] = f_2 * fg_46[k]
                  + pb_x[k] * fh_84[k];

        t_70[k] = pa_z[k] * di_15[k];

        t_71[k] = f_2 * dh_28[k]
                  + pb_z[k] * fh_85[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_y, dh_29, dh_30, dh_31, dh_39, \
                         di_16, di_17, di_18, fh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * dh_29[k]
                  + pa_z[k] * di_16[k];

        t_73[k] = f_0 * dh_30[k]
                  + pa_z[k] * di_17[k];

        t_74[k] = f_4 * dh_31[k]
                  + pa_z[k] * di_18[k];

        t_75[k] = f_3 * dh_39[k]
                  + pb_y[k] * fh_90[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pi_4, di_22, fg_47, fg_48, fg_49, \
                         fh_91, fh_92, fh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * pi_4[k]
                  + pa_y[k] * di_22[k];

        t_77[k] = f_0 * fg_47[k]
                  + pb_x[k] * fh_91[k];

        t_78[k] = f_3 * fg_48[k]
                  + pb_x[k] * fh_92[k];

        t_79[k] = f_2 * fg_49[k]
                  + pb_x[k] * fh_93[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pb_x, pb_z, dh_34, dh_50, dh_52, di_31, \
                         di_32, fg_50, fh_94, fh_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_2 * fg_50[k]
                  + pb_x[k] * fh_94[k];

        t_81[k] = f_5 * dh_50[k]
                  + pa_y[k] * di_31[k];

        t_82[k] = f_3 * dh_34[k]
                  + pb_z[k] * fh_95[k];

        t_83[k] = f_4 * dh_52[k]
                  + pa_y[k] * di_32[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_y, pb_y, dh_53, dh_54, dh_55, di_33, \
                         di_34, di_35, fh_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_0 * dh_53[k]
                  + pa_y[k] * di_33[k];

        t_85[k] = f_3 * dh_54[k]
                  + pa_y[k] * di_34[k];

        t_86[k] = f_2 * dh_55[k]
                  + pb_y[k] * fh_100[k];

        t_87[k] = pa_y[k] * di_35[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, pb_x, pb_y, fg_53, fg_54, fg_55, \
                         fg_56, fh_101, fh_103, fh_104, fh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_1 * fg_53[k]
                  + pb_x[k] * fh_101[k];

        t_89[k] = pb_y[k] * fh_101[k];

        t_90[k] = f_4 * fg_54[k]
                  + pb_x[k] * fh_103[k];

        t_91[k] = f_0 * fg_55[k]
                  + pb_x[k] * fh_104[k];

        t_92[k] = pb_y[k] * fh_103[k];

        t_93[k] = f_0 * fg_56[k]
                  + pb_x[k] * fh_106[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, fg_57, fg_58, fg_59, fg_60, \
                         fh_106, fh_107, fh_108, fh_110, fh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_3 * fg_57[k]
                  + pb_x[k] * fh_107[k];

        t_95[k] = f_3 * fg_58[k]
                  + pb_x[k] * fh_108[k];

        t_96[k] = pb_y[k] * fh_106[k];

        t_97[k] = f_3 * fg_59[k]
                  + pb_x[k] * fh_110[k];

        t_98[k] = f_2 * fg_60[k]
                  + pb_x[k] * fh_111[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pb_x, pb_y, fg_60, fg_61, fg_62, \
                         fg_64, fh_110, fh_112, fh_113, fh_114, \
                         fh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_2 * fg_61[k]
                  + pb_x[k] * fh_112[k];

        t_100[k] = f_2 * fg_62[k]
                   + pb_x[k] * fh_113[k];

        t_101[k] = pb_y[k] * fh_110[k];

        t_102[k] = f_2 * fg_64[k]
                   + pb_x[k] * fh_114[k];

        t_103[k] = f_1 * fg_60[k]
                   + pb_y[k] * fh_115[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, fg_61, fg_62, fg_63, fg_64, fh_116, \
                         fh_117, fh_118, fh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_4 * fg_61[k]
                   + pb_y[k] * fh_116[k];

        t_105[k] = f_0 * fg_62[k]
                   + pb_y[k] * fh_117[k];

        t_106[k] = f_3 * fg_63[k]
                   + pb_y[k] * fh_118[k];

        t_107[k] = f_2 * fg_64[k]
                   + pb_y[k] * fh_119[k];
    }

#pragma omp simd aligned(t_108, pb_z, dh_55, fg_64, fh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * dh_55[k]
                   + f_1 * fg_64[k]
                   + pb_z[k] * fh_120[k];
    }
}

auto
compute_prim_fi_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pi, const size_t dh, const size_t di,
                          const size_t fg, const size_t fh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_4 = buffer.data(pi + 4);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dh_0, fg_0, fg_1, fh_0, fh_1, \
                         fh_2, fh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_0[k]
                 + f_1 * fg_0[k]
                 + pb_x[k] * fh_0[k];

        t_1[k] = f_2 * fg_0[k]
                 + pb_y[k] * fh_1[k];

        t_2[k] = f_2 * fg_0[k]
                 + pb_z[k] * fh_2[k];

        t_3[k] = f_3 * fg_1[k]
                 + pb_y[k] * fh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, di_0, fg_2, fg_3, fg_4, \
                         fh_4, fh_5, fh_7, fh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * fg_2[k]
                 + pb_z[k] * fh_4[k];

        t_5[k] = f_0 * fg_3[k]
                 + pb_y[k] * fh_5[k];

        t_6[k] = f_2 * fg_4[k]
                 + pb_y[k] * fh_7[k];

        t_7[k] = f_0 * fg_4[k]
                 + pb_z[k] * fh_8[k];

        t_8[k] = pa_y[k] * di_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_x, pa_z, pi_0, dh_1, dh_3, dh_5, \
                         di_0, di_1, di_2, di_3, di_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * pi_0[k]
                 + pa_x[k] * di_4[k];

        t_10[k] = pa_z[k] * di_0[k];

        t_11[k] = f_3 * dh_1[k]
                  + pa_z[k] * di_1[k];

        t_12[k] = f_0 * dh_3[k]
                  + pa_z[k] * di_2[k];

        t_13[k] = f_4 * dh_5[k]
                  + pa_z[k] * di_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pa_x, pi_4, di_5, di_6, di_10, \
                         di_11, di_12, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_2 * pi_4[k]
                  + pa_x[k] * di_5[k];

        t_15[k] = pa_x[k] * di_6[k];

        t_16[k] = pa_x[k] * di_10[k];

        t_17[k] = pa_x[k] * di_11[k];

        t_18[k] = pa_x[k] * di_12[k];

        t_19[k] = pa_x[k] * di_18[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_x, fg_20, fg_21, fg_23, fg_25, \
                         fg_27, fh_33, fh_34, fh_36, fh_39, fh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * fg_20[k]
                  + pb_x[k] * fh_33[k];

        t_21[k] = f_0 * fg_21[k]
                  + pb_x[k] * fh_34[k];

        t_22[k] = f_3 * fg_23[k]
                  + pb_x[k] * fh_36[k];

        t_23[k] = f_2 * fg_25[k]
                  + pb_x[k] * fh_39[k];

        t_24[k] = f_2 * fg_27[k]
                  + pb_x[k] * fh_40[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_y, pb_z, dh_11, fg_25, fg_26, fg_27, \
                         fh_41, fh_42, fh_43, fh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * dh_11[k]
                  + f_1 * fg_25[k]
                  + pb_y[k] * fh_41[k];

        t_26[k] = f_2 * fg_25[k]
                  + pb_z[k] * fh_42[k];

        t_27[k] = f_3 * fg_26[k]
                  + pb_z[k] * fh_43[k];

        t_28[k] = f_0 * fg_27[k]
                  + pb_z[k] * fh_44[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_y, pa_z, pi_4, dh_12, dh_13, dh_14, \
                         di_6, di_7, di_8, di_9, di_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * di_6[k];

        t_30[k] = f_3 * dh_12[k]
                  + pa_z[k] * di_7[k];

        t_31[k] = f_0 * dh_13[k]
                  + pa_z[k] * di_8[k];

        t_32[k] = f_4 * dh_14[k]
                  + pa_z[k] * di_9[k];

        t_33[k] = f_2 * pi_4[k]
                  + pa_y[k] * di_13[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, dh_23, dh_24, dh_25, dh_26, \
                         di_14, di_15, di_16, di_17, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * dh_23[k]
                  + pa_y[k] * di_14[k];

        t_35[k] = f_4 * dh_24[k]
                  + pa_y[k] * di_15[k];

        t_36[k] = f_0 * dh_25[k]
                  + pa_y[k] * di_16[k];

        t_37[k] = f_3 * dh_26[k]
                  + pa_y[k] * di_17[k];

        t_38[k] = pa_y[k] * di_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pb_x, fg_39, fg_41, fg_42, fg_43, \
                         fg_44, fh_63, fh_65, fh_66, fh_67, fh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * fg_39[k]
                  + pb_x[k] * fh_63[k];

        t_40[k] = f_0 * fg_41[k]
                  + pb_x[k] * fh_65[k];

        t_41[k] = f_0 * fg_42[k]
                  + pb_x[k] * fh_66[k];

        t_42[k] = f_3 * fg_43[k]
                  + pb_x[k] * fh_67[k];

        t_43[k] = f_3 * fg_44[k]
                  + pb_x[k] * fh_69[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pb_x, pb_y, fg_45, fg_47, fg_49, fh_70, \
                         fh_71, fh_72, fh_73, fh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * fg_45[k]
                  + pb_x[k] * fh_70[k];

        t_45[k] = f_2 * fg_47[k]
                  + pb_x[k] * fh_71[k];

        t_46[k] = f_2 * fg_49[k]
                  + pb_x[k] * fh_72[k];

        t_47[k] = f_1 * fg_45[k]
                  + pb_y[k] * fh_73[k];

        t_48[k] = f_0 * fg_47[k]
                  + pb_y[k] * fh_75[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_y, pb_z, dh_27, fg_48, fg_49, fh_76, fh_77, \
                         fh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * fg_48[k]
                  + pb_y[k] * fh_76[k];

        t_50[k] = f_2 * fg_49[k]
                  + pb_y[k] * fh_77[k];

        t_51[k] = f_0 * dh_27[k]
                  + f_1 * fg_49[k]
                  + pb_z[k] * fh_78[k];
    }
}

auto
compute_prim_fi_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pi, const size_t dh, const size_t di,
                          const size_t fg, const size_t fh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_4 = buffer.data(pi + 4);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, dh_0, fg_0, fg_2, fg_4, fh_0, fh_1, \
                         fh_3, fh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_0[k]
                 + f_1 * fg_0[k]
                 + pb_x[k] * fh_0[k];

        t_1[k] = f_2 * fg_0[k]
                 + pb_z[k] * fh_1[k];

        t_2[k] = f_3 * fg_2[k]
                 + pb_z[k] * fh_3[k];

        t_3[k] = f_0 * fg_4[k]
                 + pb_z[k] * fh_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, pi_0, dh_1, dh_2, di_0, \
                         di_1, di_2, di_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_y[k] * di_0[k];

        t_5[k] = f_2 * pi_0[k]
                 + pa_x[k] * di_4[k];

        t_6[k] = pa_z[k] * di_0[k];

        t_7[k] = f_3 * dh_1[k]
                 + pa_z[k] * di_1[k];

        t_8[k] = f_0 * dh_2[k]
                 + pa_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, t_14, pa_x, pa_z, pi_4, dh_3, di_3, \
                         di_5, di_6, di_10, di_11, di_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dh_3[k]
                 + pa_z[k] * di_3[k];

        t_10[k] = f_2 * pi_4[k]
                  + pa_x[k] * di_5[k];

        t_11[k] = pa_x[k] * di_6[k];

        t_12[k] = pa_x[k] * di_10[k];

        t_13[k] = pa_x[k] * di_11[k];

        t_14[k] = pa_x[k] * di_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_y, pb_z, dh_8, di_18, fg_21, fg_22, \
                         fh_25, fh_26, fh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * di_18[k];

        t_16[k] = f_0 * dh_8[k]
                  + f_1 * fg_21[k]
                  + pb_y[k] * fh_25[k];

        t_17[k] = f_2 * fg_21[k]
                  + pb_z[k] * fh_26[k];

        t_18[k] = f_3 * fg_22[k]
                  + pb_z[k] * fh_27[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_z, pb_z, dh_9, dh_10, dh_11, di_6, \
                         di_7, di_8, di_9, fg_23, fh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * fg_23[k]
                  + pb_z[k] * fh_28[k];

        t_20[k] = pa_z[k] * di_6[k];

        t_21[k] = f_3 * dh_9[k]
                  + pa_z[k] * di_7[k];

        t_22[k] = f_0 * dh_10[k]
                  + pa_z[k] * di_8[k];

        t_23[k] = f_4 * dh_11[k]
                  + pa_z[k] * di_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pi_4, dh_17, dh_18, dh_19, dh_20, \
                         di_13, di_14, di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_2 * pi_4[k]
                  + pa_y[k] * di_13[k];

        t_25[k] = f_5 * dh_17[k]
                  + pa_y[k] * di_14[k];

        t_26[k] = f_4 * dh_18[k]
                  + pa_y[k] * di_15[k];

        t_27[k] = f_0 * dh_19[k]
                  + pa_y[k] * di_16[k];

        t_28[k] = f_3 * dh_20[k]
                  + pa_y[k] * di_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_y, pb_y, di_18, fg_37, fg_38, fg_39, \
                         fg_40, fh_43, fh_44, fh_45, fh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * di_18[k];

        t_30[k] = f_1 * fg_37[k]
                  + pb_y[k] * fh_43[k];

        t_31[k] = f_0 * fg_38[k]
                  + pb_y[k] * fh_44[k];

        t_32[k] = f_3 * fg_39[k]
                  + pb_y[k] * fh_45[k];

        t_33[k] = f_2 * fg_40[k]
                  + pb_y[k] * fh_46[k];
    }

#pragma omp simd aligned(t_34, pb_z, dh_21, fg_40, fh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dh_21[k]
                  + f_1 * fg_40[k]
                  + pb_z[k] * fh_47[k];
    }
}

}  // namespace simdovl
