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


#include "SimdElectronRepulsionVrrRecDL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_dl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pk, const size_t pl,
                                     const size_t di0, const size_t di1, const size_t dk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);
    const auto *pl_24 = buffer.data(pl + 24);
    const auto *pl_25 = buffer.data(pl + 25);
    const auto *pl_26 = buffer.data(pl + 26);
    const auto *pl_27 = buffer.data(pl + 27);
    const auto *pl_28 = buffer.data(pl + 28);
    const auto *pl_29 = buffer.data(pl + 29);
    const auto *pl_30 = buffer.data(pl + 30);
    const auto *pl_31 = buffer.data(pl + 31);
    const auto *pl_32 = buffer.data(pl + 32);
    const auto *pl_33 = buffer.data(pl + 33);
    const auto *pl_34 = buffer.data(pl + 34);
    const auto *pl_35 = buffer.data(pl + 35);
    const auto *pl_36 = buffer.data(pl + 36);
    const auto *pl_37 = buffer.data(pl + 37);
    const auto *pl_38 = buffer.data(pl + 38);
    const auto *pl_39 = buffer.data(pl + 39);
    const auto *pl_40 = buffer.data(pl + 40);
    const auto *pl_41 = buffer.data(pl + 41);
    const auto *pl_42 = buffer.data(pl + 42);
    const auto *pl_43 = buffer.data(pl + 43);
    const auto *pl_44 = buffer.data(pl + 44);
    const auto *pl_45 = buffer.data(pl + 45);
    const auto *pl_46 = buffer.data(pl + 46);
    const auto *pl_47 = buffer.data(pl + 47);
    const auto *pl_48 = buffer.data(pl + 48);
    const auto *pl_49 = buffer.data(pl + 49);
    const auto *pl_50 = buffer.data(pl + 50);
    const auto *pl_51 = buffer.data(pl + 51);
    const auto *pl_52 = buffer.data(pl + 52);
    const auto *pl_53 = buffer.data(pl + 53);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);
    const auto *di0_3 = buffer.data(di0 + 3);
    const auto *di0_4 = buffer.data(di0 + 4);
    const auto *di0_5 = buffer.data(di0 + 5);
    const auto *di0_6 = buffer.data(di0 + 6);
    const auto *di0_7 = buffer.data(di0 + 7);
    const auto *di0_8 = buffer.data(di0 + 8);
    const auto *di0_9 = buffer.data(di0 + 9);
    const auto *di0_10 = buffer.data(di0 + 10);
    const auto *di0_11 = buffer.data(di0 + 11);
    const auto *di0_12 = buffer.data(di0 + 12);
    const auto *di0_13 = buffer.data(di0 + 13);
    const auto *di0_14 = buffer.data(di0 + 14);
    const auto *di0_15 = buffer.data(di0 + 15);
    const auto *di0_16 = buffer.data(di0 + 16);
    const auto *di0_17 = buffer.data(di0 + 17);
    const auto *di0_18 = buffer.data(di0 + 18);
    const auto *di0_19 = buffer.data(di0 + 19);
    const auto *di0_20 = buffer.data(di0 + 20);
    const auto *di0_21 = buffer.data(di0 + 21);
    const auto *di0_22 = buffer.data(di0 + 22);
    const auto *di0_23 = buffer.data(di0 + 23);
    const auto *di0_24 = buffer.data(di0 + 24);
    const auto *di0_25 = buffer.data(di0 + 25);
    const auto *di0_26 = buffer.data(di0 + 26);
    const auto *di0_27 = buffer.data(di0 + 27);
    const auto *di0_28 = buffer.data(di0 + 28);
    const auto *di0_29 = buffer.data(di0 + 29);
    const auto *di0_30 = buffer.data(di0 + 30);
    const auto *di0_31 = buffer.data(di0 + 31);
    const auto *di0_32 = buffer.data(di0 + 32);
    const auto *di0_33 = buffer.data(di0 + 33);
    const auto *di0_34 = buffer.data(di0 + 34);
    const auto *di0_35 = buffer.data(di0 + 35);
    const auto *di0_36 = buffer.data(di0 + 36);
    const auto *di0_37 = buffer.data(di0 + 37);
    const auto *di0_38 = buffer.data(di0 + 38);
    const auto *di0_39 = buffer.data(di0 + 39);
    const auto *di0_40 = buffer.data(di0 + 40);
    const auto *di0_41 = buffer.data(di0 + 41);
    const auto *di0_42 = buffer.data(di0 + 42);
    const auto *di0_43 = buffer.data(di0 + 43);
    const auto *di0_44 = buffer.data(di0 + 44);
    const auto *di0_45 = buffer.data(di0 + 45);
    const auto *di0_46 = buffer.data(di0 + 46);
    const auto *di0_47 = buffer.data(di0 + 47);
    const auto *di0_48 = buffer.data(di0 + 48);
    const auto *di0_49 = buffer.data(di0 + 49);
    const auto *di0_50 = buffer.data(di0 + 50);
    const auto *di0_51 = buffer.data(di0 + 51);
    const auto *di0_52 = buffer.data(di0 + 52);
    const auto *di0_53 = buffer.data(di0 + 53);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_1 = buffer.data(di1 + 1);
    const auto *di1_2 = buffer.data(di1 + 2);
    const auto *di1_3 = buffer.data(di1 + 3);
    const auto *di1_4 = buffer.data(di1 + 4);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_6 = buffer.data(di1 + 6);
    const auto *di1_7 = buffer.data(di1 + 7);
    const auto *di1_8 = buffer.data(di1 + 8);
    const auto *di1_9 = buffer.data(di1 + 9);
    const auto *di1_10 = buffer.data(di1 + 10);
    const auto *di1_11 = buffer.data(di1 + 11);
    const auto *di1_12 = buffer.data(di1 + 12);
    const auto *di1_13 = buffer.data(di1 + 13);
    const auto *di1_14 = buffer.data(di1 + 14);
    const auto *di1_15 = buffer.data(di1 + 15);
    const auto *di1_16 = buffer.data(di1 + 16);
    const auto *di1_17 = buffer.data(di1 + 17);
    const auto *di1_18 = buffer.data(di1 + 18);
    const auto *di1_19 = buffer.data(di1 + 19);
    const auto *di1_20 = buffer.data(di1 + 20);
    const auto *di1_21 = buffer.data(di1 + 21);
    const auto *di1_22 = buffer.data(di1 + 22);
    const auto *di1_23 = buffer.data(di1 + 23);
    const auto *di1_24 = buffer.data(di1 + 24);
    const auto *di1_25 = buffer.data(di1 + 25);
    const auto *di1_26 = buffer.data(di1 + 26);
    const auto *di1_27 = buffer.data(di1 + 27);
    const auto *di1_28 = buffer.data(di1 + 28);
    const auto *di1_29 = buffer.data(di1 + 29);
    const auto *di1_30 = buffer.data(di1 + 30);
    const auto *di1_31 = buffer.data(di1 + 31);
    const auto *di1_32 = buffer.data(di1 + 32);
    const auto *di1_33 = buffer.data(di1 + 33);
    const auto *di1_34 = buffer.data(di1 + 34);
    const auto *di1_35 = buffer.data(di1 + 35);
    const auto *di1_36 = buffer.data(di1 + 36);
    const auto *di1_37 = buffer.data(di1 + 37);
    const auto *di1_38 = buffer.data(di1 + 38);
    const auto *di1_39 = buffer.data(di1 + 39);
    const auto *di1_40 = buffer.data(di1 + 40);
    const auto *di1_41 = buffer.data(di1 + 41);
    const auto *di1_42 = buffer.data(di1 + 42);
    const auto *di1_43 = buffer.data(di1 + 43);
    const auto *di1_44 = buffer.data(di1 + 44);
    const auto *di1_45 = buffer.data(di1 + 45);
    const auto *di1_46 = buffer.data(di1 + 46);
    const auto *di1_47 = buffer.data(di1 + 47);
    const auto *di1_48 = buffer.data(di1 + 48);
    const auto *di1_49 = buffer.data(di1 + 49);
    const auto *di1_50 = buffer.data(di1 + 50);
    const auto *di1_51 = buffer.data(di1 + 51);
    const auto *di1_52 = buffer.data(di1 + 52);
    const auto *di1_53 = buffer.data(di1 + 53);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, pk_0, di0_0, di1_0, \
                         dk_0, dk_1, dk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = pb_y[k] * dk_0[k];

        t_2[k] = pb_z[k] * dk_0[k];

        t_3[k] = f_3 * di0_0[k]
                 - f_4 * di1_0[k]
                 + pb_y[k] * dk_1[k];

        t_4[k] = pb_y[k] * dk_2[k];

        t_5[k] = f_3 * di0_0[k]
                 - f_4 * di1_0[k]
                 + pb_z[k] * dk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, di0_1, di0_2, di0_3, di1_1, \
                         di1_2, di1_3, dk_3, dk_4, dk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * di0_1[k]
                 - f_6 * di1_1[k]
                 + pb_y[k] * dk_3[k];

        t_7[k] = pb_z[k] * dk_3[k];

        t_8[k] = pb_y[k] * dk_4[k];

        t_9[k] = f_5 * di0_2[k]
                 - f_6 * di1_2[k]
                 + pb_z[k] * dk_4[k];

        t_10[k] = f_7 * di0_3[k]
                  - f_8 * di1_3[k]
                  + pb_y[k] * dk_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, di0_4, di0_5, di1_4, \
                         di1_5, dk_5, dk_6, dk_7, dk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * dk_5[k];

        t_12[k] = f_3 * di0_4[k]
                  - f_4 * di1_4[k]
                  + pb_y[k] * dk_6[k];

        t_13[k] = pb_y[k] * dk_7[k];

        t_14[k] = f_7 * di0_4[k]
                  - f_8 * di1_4[k]
                  + pb_z[k] * dk_7[k];

        t_15[k] = f_9 * di0_5[k]
                  - f_10 * di1_5[k]
                  + pb_y[k] * dk_8[k];

        t_16[k] = pb_z[k] * dk_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, di0_6, di0_7, di1_6, di1_7, dk_9, \
                         dk_10, dk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * di0_6[k]
                  - f_6 * di1_6[k]
                  + pb_y[k] * dk_9[k];

        t_18[k] = f_3 * di0_7[k]
                  - f_4 * di1_7[k]
                  + pb_y[k] * dk_10[k];

        t_19[k] = pb_y[k] * dk_11[k];

        t_20[k] = f_9 * di0_7[k]
                  - f_10 * di1_7[k]
                  + pb_z[k] * dk_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, di0_8, di0_9, di0_10, di1_8, \
                         di1_9, di1_10, dk_12, dk_13, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * di0_8[k]
                  - f_12 * di1_8[k]
                  + pb_y[k] * dk_12[k];

        t_22[k] = pb_z[k] * dk_12[k];

        t_23[k] = f_7 * di0_9[k]
                  - f_8 * di1_9[k]
                  + pb_y[k] * dk_13[k];

        t_24[k] = f_5 * di0_10[k]
                  - f_6 * di1_10[k]
                  + pb_y[k] * dk_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, pk_9, di0_11, di1_11, \
                         dk_15, dk_16, dk_17, dk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * di0_11[k]
                  - f_4 * di1_11[k]
                  + pb_y[k] * dk_15[k];

        t_26[k] = pb_y[k] * dk_16[k];

        t_27[k] = f_11 * di0_11[k]
                  - f_12 * di1_11[k]
                  + pb_z[k] * dk_16[k];

        t_28[k] = f_0 * pk_9[k]
                  + pb_x[k] * dk_19[k];

        t_29[k] = pb_z[k] * dk_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, pk_10, pk_11, pk_12, pk_13, \
                         dk_18, dk_20, dk_21, dk_22, dk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * pk_10[k]
                  + pb_x[k] * dk_20[k];

        t_31[k] = f_0 * pk_11[k]
                  + pb_x[k] * dk_21[k];

        t_32[k] = f_0 * pk_12[k]
                  + pb_x[k] * dk_22[k];

        t_33[k] = f_0 * pk_13[k]
                  + pb_x[k] * dk_23[k];

        t_34[k] = pb_y[k] * dk_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, pk_14, di0_12, di0_13, \
                         di1_12, di1_13, dk_19, dk_20, dk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * pk_14[k]
                  + pb_x[k] * dk_25[k];

        t_36[k] = f_1 * di0_12[k]
                  - f_2 * di1_12[k]
                  + pb_y[k] * dk_19[k];

        t_37[k] = pb_z[k] * dk_19[k];

        t_38[k] = f_11 * di0_13[k]
                  - f_12 * di1_13[k]
                  + pb_y[k] * dk_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, di0_14, di0_15, di0_16, di1_14, di1_15, \
                         di1_16, dk_21, dk_22, dk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * di0_14[k]
                  - f_10 * di1_14[k]
                  + pb_y[k] * dk_21[k];

        t_40[k] = f_7 * di0_15[k]
                  - f_8 * di1_15[k]
                  + pb_y[k] * dk_22[k];

        t_41[k] = f_5 * di0_16[k]
                  - f_6 * di1_16[k]
                  + pb_y[k] * dk_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, pk_0, pl_0, \
                         di0_17, di1_17, dk_24, dk_25, dk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * di0_17[k]
                  - f_4 * di1_17[k]
                  + pb_y[k] * dk_24[k];

        t_43[k] = pb_y[k] * dk_25[k];

        t_44[k] = f_1 * di0_17[k]
                  - f_2 * di1_17[k]
                  + pb_z[k] * dk_25[k];

        t_45[k] = pa_y[k] * pl_0[k];

        t_46[k] = f_13 * pk_0[k]
                  + pb_y[k] * dk_26[k];

        t_47[k] = pb_z[k] * dk_26[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pa_y, pb_z, pk_16, pk_18, pl_2, \
                         pl_14, pl_15, dk_27, dk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * pk_16[k]
                  + pa_x[k] * pl_14[k];

        t_49[k] = pb_z[k] * dk_27[k];

        t_50[k] = pa_y[k] * pl_2[k];

        t_51[k] = f_15 * pk_18[k]
                  + pa_x[k] * pl_15[k];

        t_52[k] = pb_z[k] * dk_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pa_y, pb_y, pb_z, pk_2, pk_20, pl_4, \
                         pl_16, dk_29, dk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * pk_2[k]
                  + pb_y[k] * dk_29[k];

        t_54[k] = pa_y[k] * pl_4[k];

        t_55[k] = f_16 * pk_20[k]
                  + pa_x[k] * pl_16[k];

        t_56[k] = pb_z[k] * dk_30[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pa_y, pb_y, pk_4, pk_21, pk_23, pl_6, \
                         pl_17, pl_18, dk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_16 * pk_21[k]
                  + pa_x[k] * pl_17[k];

        t_58[k] = f_13 * pk_4[k]
                  + pb_y[k] * dk_31[k];

        t_59[k] = pa_y[k] * pl_6[k];

        t_60[k] = f_17 * pk_23[k]
                  + pa_x[k] * pl_18[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_y, pb_z, pk_6, pk_24, pk_25, pl_19, \
                         pl_20, dk_32, dk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * dk_32[k];

        t_62[k] = f_17 * pk_24[k]
                  + pa_x[k] * pl_19[k];

        t_63[k] = f_17 * pk_25[k]
                  + pa_x[k] * pl_20[k];

        t_64[k] = f_13 * pk_6[k]
                  + pb_y[k] * dk_33[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pa_y, pb_z, pk_27, pk_28, pk_29, \
                         pl_8, pl_21, pl_22, pl_23, dk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_y[k] * pl_8[k];

        t_66[k] = f_0 * pk_27[k]
                  + pa_x[k] * pl_21[k];

        t_67[k] = pb_z[k] * dk_34[k];

        t_68[k] = f_0 * pk_28[k]
                  + pa_x[k] * pl_22[k];

        t_69[k] = f_0 * pk_29[k]
                  + pa_x[k] * pl_23[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pa_y, pb_x, pb_y, pk_8, pk_30, pk_31, \
                         pl_10, pl_24, dk_35, dk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * pk_30[k]
                  + pa_x[k] * pl_24[k];

        t_71[k] = f_13 * pk_8[k]
                  + pb_y[k] * dk_35[k];

        t_72[k] = pa_y[k] * pl_10[k];

        t_73[k] = f_13 * pk_31[k]
                  + pb_x[k] * dk_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pb_x, pb_z, pk_32, pk_33, pk_34, pk_35, \
                         dk_36, dk_38, dk_39, dk_40, dk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_z[k] * dk_36[k];

        t_75[k] = f_13 * pk_32[k]
                  + pb_x[k] * dk_38[k];

        t_76[k] = f_13 * pk_33[k]
                  + pb_x[k] * dk_39[k];

        t_77[k] = f_13 * pk_34[k]
                  + pb_x[k] * dk_40[k];

        t_78[k] = f_13 * pk_35[k]
                  + pb_x[k] * dk_41[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_x, pa_y, pb_x, pb_z, pk_36, pl_12, \
                         pl_25, pl_26, dk_37, dk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_13 * pk_36[k]
                  + pb_x[k] * dk_42[k];

        t_80[k] = pa_y[k] * pl_12[k];

        t_81[k] = pa_x[k] * pl_25[k];

        t_82[k] = pb_z[k] * dk_37[k];

        t_83[k] = pa_x[k] * pl_26[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, t_90, pa_x, pa_z, pl_0, pl_27, \
                         pl_28, pl_29, pl_30, pl_31, pl_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_x[k] * pl_27[k];

        t_85[k] = pa_x[k] * pl_28[k];

        t_86[k] = pa_x[k] * pl_29[k];

        t_87[k] = pa_x[k] * pl_30[k];

        t_88[k] = pa_x[k] * pl_31[k];

        t_89[k] = pa_x[k] * pl_32[k];

        t_90[k] = pa_z[k] * pl_0[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pa_x, pa_z, pb_y, pb_z, pk_0, pk_41, \
                         pl_1, pl_35, dk_43, dk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_y[k] * dk_43[k];

        t_92[k] = f_13 * pk_0[k]
                  + pb_z[k] * dk_43[k];

        t_93[k] = pa_z[k] * pl_1[k];

        t_94[k] = pb_y[k] * dk_44[k];

        t_95[k] = f_14 * pk_41[k]
                  + pa_x[k] * pl_35[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_x, pa_z, pb_y, pb_z, pk_1, pk_44, \
                         pl_3, pl_5, pl_36, dk_45, dk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * pl_3[k];

        t_97[k] = f_13 * pk_1[k]
                  + pb_z[k] * dk_45[k];

        t_98[k] = pb_y[k] * dk_46[k];

        t_99[k] = f_15 * pk_44[k]
                  + pa_x[k] * pl_36[k];

        t_100[k] = pa_z[k] * pl_5[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pb_y, pb_z, pk_3, pk_46, pk_48, \
                         pl_37, pl_38, dk_47, dk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_13 * pk_3[k]
                   + pb_z[k] * dk_47[k];

        t_102[k] = f_16 * pk_46[k]
                   + pa_x[k] * pl_37[k];

        t_103[k] = pb_y[k] * dk_48[k];

        t_104[k] = f_16 * pk_48[k]
                   + pa_x[k] * pl_38[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pa_z, pb_z, pk_5, pk_50, pk_51, \
                         pl_7, pl_39, pl_40, dk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * pl_7[k];

        t_106[k] = f_13 * pk_5[k]
                   + pb_z[k] * dk_49[k];

        t_107[k] = f_17 * pk_50[k]
                   + pa_x[k] * pl_39[k];

        t_108[k] = f_17 * pk_51[k]
                   + pa_x[k] * pl_40[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_z, pb_y, pb_z, pk_7, pk_53, \
                         pl_9, pl_41, dk_50, dk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_y[k] * dk_50[k];

        t_110[k] = f_17 * pk_53[k]
                   + pa_x[k] * pl_41[k];

        t_111[k] = pa_z[k] * pl_9[k];

        t_112[k] = f_13 * pk_7[k]
                   + pb_z[k] * dk_51[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_x, pb_y, pk_54, pk_55, pk_56, \
                         pk_57, pl_42, pl_43, pl_44, pl_45, dk_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_0 * pk_54[k]
                   + pa_x[k] * pl_42[k];

        t_114[k] = f_0 * pk_55[k]
                   + pa_x[k] * pl_43[k];

        t_115[k] = f_0 * pk_56[k]
                   + pa_x[k] * pl_44[k];

        t_116[k] = pb_y[k] * dk_52[k];

        t_117[k] = f_0 * pk_57[k]
                   + pa_x[k] * pl_45[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, pk_59, pk_60, pk_61, \
                         pk_62, pl_11, dk_54, dk_55, dk_56, dk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * pl_11[k];

        t_119[k] = f_13 * pk_59[k]
                   + pb_x[k] * dk_54[k];

        t_120[k] = f_13 * pk_60[k]
                   + pb_x[k] * dk_55[k];

        t_121[k] = f_13 * pk_61[k]
                   + pb_x[k] * dk_56[k];

        t_122[k] = f_13 * pk_62[k]
                   + pb_x[k] * dk_57[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, pa_x, pb_x, pb_y, pk_63, pk_65, \
                         pl_46, pl_47, dk_53, dk_58, dk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_13 * pk_63[k]
                   + pb_x[k] * dk_58[k];

        t_124[k] = pb_y[k] * dk_53[k];

        t_125[k] = f_13 * pk_65[k]
                   + pb_x[k] * dk_59[k];

        t_126[k] = pa_x[k] * pl_46[k];

        t_127[k] = pa_x[k] * pl_47[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, t_134, pa_x, pb_y, pl_48, \
                         pl_49, pl_50, pl_51, pl_52, pl_53, dk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pa_x[k] * pl_48[k];

        t_129[k] = pa_x[k] * pl_49[k];

        t_130[k] = pa_x[k] * pl_50[k];

        t_131[k] = pa_x[k] * pl_51[k];

        t_132[k] = pa_x[k] * pl_52[k];

        t_133[k] = pb_y[k] * dk_59[k];

        t_134[k] = pa_x[k] * pl_53[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pb_x, pb_y, pb_z, pk_15, di0_18, \
                         di0_19, di1_18, di1_19, dk_60, dk_61, dk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * di0_18[k]
                   - f_2 * di1_18[k]
                   + pb_x[k] * dk_60[k];

        t_136[k] = f_0 * pk_15[k]
                   + pb_y[k] * dk_60[k];

        t_137[k] = pb_z[k] * dk_60[k];

        t_138[k] = f_11 * di0_19[k]
                   - f_12 * di1_19[k]
                   + pb_x[k] * dk_62[k];

        t_139[k] = pb_z[k] * dk_61[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, pb_y, pb_z, pk_17, di0_20, di0_21, \
                         di1_20, di1_21, dk_62, dk_63, dk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_11 * di0_20[k]
                   - f_12 * di1_20[k]
                   + pb_x[k] * dk_63[k];

        t_141[k] = f_9 * di0_21[k]
                   - f_10 * di1_21[k]
                   + pb_x[k] * dk_64[k];

        t_142[k] = pb_z[k] * dk_62[k];

        t_143[k] = f_0 * pk_17[k]
                   + pb_y[k] * dk_63[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_z, di0_22, di0_23, di0_24, \
                         di1_22, di1_23, di1_24, dk_64, dk_65, dk_66, \
                         dk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_9 * di0_22[k]
                   - f_10 * di1_22[k]
                   + pb_x[k] * dk_65[k];

        t_145[k] = f_7 * di0_23[k]
                   - f_8 * di1_23[k]
                   + pb_x[k] * dk_66[k];

        t_146[k] = pb_z[k] * dk_64[k];

        t_147[k] = f_7 * di0_24[k]
                   - f_8 * di1_24[k]
                   + pb_x[k] * dk_67[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pb_y, pb_z, pk_19, di0_25, di0_26, \
                         di1_25, di1_26, dk_65, dk_66, dk_68, dk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * pk_19[k]
                   + pb_y[k] * dk_65[k];

        t_149[k] = f_7 * di0_25[k]
                   - f_8 * di1_25[k]
                   + pb_x[k] * dk_68[k];

        t_150[k] = f_5 * di0_26[k]
                   - f_6 * di1_26[k]
                   + pb_x[k] * dk_69[k];

        t_151[k] = pb_z[k] * dk_66[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_x, pb_y, pk_22, di0_27, di0_28, di1_27, \
                         di1_28, dk_68, dk_70, dk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * di0_27[k]
                   - f_6 * di1_27[k]
                   + pb_x[k] * dk_70[k];

        t_153[k] = f_5 * di0_28[k]
                   - f_6 * di1_28[k]
                   + pb_x[k] * dk_71[k];

        t_154[k] = f_0 * pk_22[k]
                   + pb_y[k] * dk_68[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pb_x, pb_z, di0_29, di0_30, di0_32, \
                         di1_29, di1_30, di1_32, dk_69, dk_72, dk_73, \
                         dk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_5 * di0_29[k]
                   - f_6 * di1_29[k]
                   + pb_x[k] * dk_72[k];

        t_156[k] = f_3 * di0_30[k]
                   - f_4 * di1_30[k]
                   + pb_x[k] * dk_73[k];

        t_157[k] = pb_z[k] * dk_69[k];

        t_158[k] = f_3 * di0_32[k]
                   - f_4 * di1_32[k]
                   + pb_x[k] * dk_74[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_x, pb_y, pk_26, di0_33, di0_34, di1_33, \
                         di1_34, dk_72, dk_75, dk_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_3 * di0_33[k]
                   - f_4 * di1_33[k]
                   + pb_x[k] * dk_75[k];

        t_160[k] = f_3 * di0_34[k]
                   - f_4 * di1_34[k]
                   + pb_x[k] * dk_76[k];

        t_161[k] = f_0 * pk_26[k]
                   + pb_y[k] * dk_72[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, pb_x, di0_35, di1_35, \
                         dk_77, dk_78, dk_79, dk_80, dk_81, dk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_3 * di0_35[k]
                   - f_4 * di1_35[k]
                   + pb_x[k] * dk_77[k];

        t_163[k] = pb_x[k] * dk_78[k];

        t_164[k] = pb_x[k] * dk_79[k];

        t_165[k] = pb_x[k] * dk_80[k];

        t_166[k] = pb_x[k] * dk_81[k];

        t_167[k] = pb_x[k] * dk_82[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pb_x, pb_y, pb_z, pk_31, di0_30, \
                         di1_30, dk_78, dk_83, dk_84, dk_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pb_x[k] * dk_83[k];

        t_169[k] = pb_x[k] * dk_84[k];

        t_170[k] = pb_x[k] * dk_85[k];

        t_171[k] = f_0 * pk_31[k]
                   + f_1 * di0_30[k]
                   - f_2 * di1_30[k]
                   + pb_y[k] * dk_78[k];

        t_172[k] = pb_z[k] * dk_78[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pb_z, di0_30, di0_31, di0_32, di1_30, di1_31, \
                         di1_32, dk_79, dk_80, dk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_3 * di0_30[k]
                   - f_4 * di1_30[k]
                   + pb_z[k] * dk_79[k];

        t_174[k] = f_5 * di0_31[k]
                   - f_6 * di1_31[k]
                   + pb_z[k] * dk_80[k];

        t_175[k] = f_7 * di0_32[k]
                   - f_8 * di1_32[k]
                   + pb_z[k] * dk_81[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pb_z, pk_37, di0_33, di0_34, \
                         di0_35, di1_33, di1_34, di1_35, dk_82, dk_83, \
                         dk_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * di0_33[k]
                   - f_10 * di1_33[k]
                   + pb_z[k] * dk_82[k];

        t_177[k] = f_11 * di0_34[k]
                   - f_12 * di1_34[k]
                   + pb_z[k] * dk_83[k];

        t_178[k] = f_0 * pk_37[k]
                   + pb_y[k] * dk_85[k];

        t_179[k] = f_1 * di0_35[k]
                   - f_2 * di1_35[k]
                   + pb_z[k] * dk_85[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pa_y, pa_z, pb_y, pk_39, \
                         pl_13, pl_14, pl_33, pl_34, pl_35, dk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * pl_33[k];

        t_181[k] = pa_z[k] * pl_13[k];

        t_182[k] = pa_y[k] * pl_34[k];

        t_183[k] = pa_z[k] * pl_14[k];

        t_184[k] = f_13 * pk_39[k]
                   + pb_y[k] * dk_86[k];

        t_185[k] = pa_y[k] * pl_35[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, pk_16, \
                         pk_41, pl_15, pl_16, pl_36, dk_87, dk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * pl_15[k];

        t_187[k] = f_13 * pk_16[k]
                   + pb_z[k] * dk_87[k];

        t_188[k] = f_13 * pk_41[k]
                   + pb_y[k] * dk_88[k];

        t_189[k] = pa_y[k] * pl_36[k];

        t_190[k] = pa_z[k] * pl_16[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, pk_18, pk_43, pk_44, \
                         pl_37, pl_38, dk_89, dk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * pk_18[k]
                   + pb_z[k] * dk_89[k];

        t_192[k] = f_0 * pk_43[k]
                   + pa_y[k] * pl_37[k];

        t_193[k] = f_13 * pk_44[k]
                   + pb_y[k] * dk_90[k];

        t_194[k] = pa_y[k] * pl_38[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pa_z, pb_z, pk_20, pk_46, pk_47, \
                         pl_18, pl_39, pl_40, dk_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * pl_18[k];

        t_196[k] = f_13 * pk_20[k]
                   + pb_z[k] * dk_91[k];

        t_197[k] = f_17 * pk_46[k]
                   + pa_y[k] * pl_39[k];

        t_198[k] = f_0 * pk_47[k]
                   + pa_y[k] * pl_40[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, pk_23, pk_48, \
                         pl_21, pl_41, dk_92, dk_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * pk_48[k]
                   + pb_y[k] * dk_92[k];

        t_200[k] = pa_y[k] * pl_41[k];

        t_201[k] = pa_z[k] * pl_21[k];

        t_202[k] = f_13 * pk_23[k]
                   + pb_z[k] * dk_93[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pa_y, pb_y, pk_50, pk_51, pk_52, \
                         pk_53, pl_42, pl_43, pl_44, pl_45, dk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_16 * pk_50[k]
                   + pa_y[k] * pl_42[k];

        t_204[k] = f_17 * pk_51[k]
                   + pa_y[k] * pl_43[k];

        t_205[k] = f_0 * pk_52[k]
                   + pa_y[k] * pl_44[k];

        t_206[k] = f_13 * pk_53[k]
                   + pb_y[k] * dk_94[k];

        t_207[k] = pa_y[k] * pl_45[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, t_213, t_214, pb_x, dk_95, dk_96, \
                         dk_97, dk_98, dk_99, dk_100, dk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pb_x[k] * dk_95[k];

        t_209[k] = pb_x[k] * dk_96[k];

        t_210[k] = pb_x[k] * dk_97[k];

        t_211[k] = pb_x[k] * dk_98[k];

        t_212[k] = pb_x[k] * dk_99[k];

        t_213[k] = pb_x[k] * dk_100[k];

        t_214[k] = pb_x[k] * dk_101[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pa_y, pa_z, pb_x, pb_z, pk_31, pk_60, \
                         pl_25, pl_48, dk_95, dk_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pb_x[k] * dk_102[k];

        t_216[k] = pa_z[k] * pl_25[k];

        t_217[k] = f_13 * pk_31[k]
                   + pb_z[k] * dk_95[k];

        t_218[k] = f_14 * pk_60[k]
                   + pa_y[k] * pl_48[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pa_y, pk_61, pk_62, pk_63, pk_64, pl_49, \
                         pl_50, pl_51, pl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_15 * pk_61[k]
                   + pa_y[k] * pl_49[k];

        t_220[k] = f_16 * pk_62[k]
                   + pa_y[k] * pl_50[k];

        t_221[k] = f_17 * pk_63[k]
                   + pa_y[k] * pl_51[k];

        t_222[k] = f_0 * pk_64[k]
                   + pa_y[k] * pl_52[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_y, pb_x, pb_y, pb_z, pk_38, \
                         pk_65, pl_53, di0_36, di1_36, dk_102, dk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_13 * pk_65[k]
                   + pb_y[k] * dk_102[k];

        t_224[k] = pa_y[k] * pl_53[k];

        t_225[k] = f_1 * di0_36[k]
                   - f_2 * di1_36[k]
                   + pb_x[k] * dk_103[k];

        t_226[k] = pb_y[k] * dk_103[k];

        t_227[k] = f_0 * pk_38[k]
                   + pb_z[k] * dk_103[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, pb_y, di0_37, di0_38, di0_39, \
                         di1_37, di1_38, di1_39, dk_104, dk_105, dk_106, \
                         dk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_11 * di0_37[k]
                   - f_12 * di1_37[k]
                   + pb_x[k] * dk_105[k];

        t_229[k] = pb_y[k] * dk_104[k];

        t_230[k] = f_11 * di0_38[k]
                   - f_12 * di1_38[k]
                   + pb_x[k] * dk_106[k];

        t_231[k] = f_9 * di0_39[k]
                   - f_10 * di1_39[k]
                   + pb_x[k] * dk_107[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pb_x, pb_y, pb_z, pk_40, di0_40, di0_41, \
                         di1_40, di1_41, dk_105, dk_106, dk_108, \
                         dk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_0 * pk_40[k]
                   + pb_z[k] * dk_105[k];

        t_233[k] = pb_y[k] * dk_106[k];

        t_234[k] = f_9 * di0_40[k]
                   - f_10 * di1_40[k]
                   + pb_x[k] * dk_108[k];

        t_235[k] = f_7 * di0_41[k]
                   - f_8 * di1_41[k]
                   + pb_x[k] * dk_109[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_x, pb_y, pb_z, pk_42, di0_42, di0_43, \
                         di1_42, di1_43, dk_107, dk_108, dk_110, \
                         dk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * pk_42[k]
                   + pb_z[k] * dk_107[k];

        t_237[k] = f_7 * di0_42[k]
                   - f_8 * di1_42[k]
                   + pb_x[k] * dk_110[k];

        t_238[k] = pb_y[k] * dk_108[k];

        t_239[k] = f_7 * di0_43[k]
                   - f_8 * di1_43[k]
                   + pb_x[k] * dk_111[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pb_x, pb_z, pk_45, di0_44, di0_45, di1_44, \
                         di1_45, dk_109, dk_112, dk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_5 * di0_44[k]
                   - f_6 * di1_44[k]
                   + pb_x[k] * dk_112[k];

        t_241[k] = f_0 * pk_45[k]
                   + pb_z[k] * dk_109[k];

        t_242[k] = f_5 * di0_45[k]
                   - f_6 * di1_45[k]
                   + pb_x[k] * dk_113[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_x, pb_y, di0_46, di0_47, di0_48, \
                         di1_46, di1_47, di1_48, dk_111, dk_114, dk_115, \
                         dk_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_5 * di0_46[k]
                   - f_6 * di1_46[k]
                   + pb_x[k] * dk_114[k];

        t_244[k] = pb_y[k] * dk_111[k];

        t_245[k] = f_5 * di0_47[k]
                   - f_6 * di1_47[k]
                   + pb_x[k] * dk_115[k];

        t_246[k] = f_3 * di0_48[k]
                   - f_4 * di1_48[k]
                   + pb_x[k] * dk_116[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pb_x, pb_z, pk_49, di0_49, di0_50, di1_49, \
                         di1_50, dk_112, dk_117, dk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * pk_49[k]
                   + pb_z[k] * dk_112[k];

        t_248[k] = f_3 * di0_49[k]
                   - f_4 * di1_49[k]
                   + pb_x[k] * dk_117[k];

        t_249[k] = f_3 * di0_50[k]
                   - f_4 * di1_50[k]
                   + pb_x[k] * dk_118[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, pb_x, pb_y, di0_51, di0_53, \
                         di1_51, di1_53, dk_115, dk_119, dk_120, dk_121, \
                         dk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_3 * di0_51[k]
                   - f_4 * di1_51[k]
                   + pb_x[k] * dk_119[k];

        t_251[k] = pb_y[k] * dk_115[k];

        t_252[k] = f_3 * di0_53[k]
                   - f_4 * di1_53[k]
                   + pb_x[k] * dk_120[k];

        t_253[k] = pb_x[k] * dk_121[k];

        t_254[k] = pb_x[k] * dk_122[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, pb_x, dk_123, dk_124, \
                         dk_125, dk_126, dk_127, dk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pb_x[k] * dk_123[k];

        t_256[k] = pb_x[k] * dk_124[k];

        t_257[k] = pb_x[k] * dk_125[k];

        t_258[k] = pb_x[k] * dk_126[k];

        t_259[k] = pb_x[k] * dk_127[k];

        t_260[k] = pb_x[k] * dk_128[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_y, pb_z, pk_58, di0_48, di0_49, \
                         di0_50, di1_48, di1_49, di1_50, dk_121, dk_123, \
                         dk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * di0_48[k]
                   - f_2 * di1_48[k]
                   + pb_y[k] * dk_121[k];

        t_262[k] = f_0 * pk_58[k]
                   + pb_z[k] * dk_121[k];

        t_263[k] = f_11 * di0_49[k]
                   - f_12 * di1_49[k]
                   + pb_y[k] * dk_123[k];

        t_264[k] = f_9 * di0_50[k]
                   - f_10 * di1_50[k]
                   + pb_y[k] * dk_124[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_y, di0_51, di0_52, di0_53, di1_51, \
                         di1_52, di1_53, dk_125, dk_126, dk_127, \
                         dk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_7 * di0_51[k]
                   - f_8 * di1_51[k]
                   + pb_y[k] * dk_125[k];

        t_266[k] = f_5 * di0_52[k]
                   - f_6 * di1_52[k]
                   + pb_y[k] * dk_126[k];

        t_267[k] = f_3 * di0_53[k]
                   - f_4 * di1_53[k]
                   + pb_y[k] * dk_127[k];

        t_268[k] = pb_y[k] * dk_128[k];
    }

#pragma omp simd aligned(t_269, pb_z, pk_65, di0_53, di1_53, dk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_0 * pk_65[k]
                   + f_1 * di0_53[k]
                   - f_2 * di1_53[k]
                   + pb_z[k] * dk_128[k];
    }
}

auto
compute_prim_dl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pk, const size_t pl,
                                     const size_t di0, const size_t di1, const size_t dk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 1.5 / p;

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

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);
    const auto *di0_3 = buffer.data(di0 + 3);
    const auto *di0_4 = buffer.data(di0 + 4);
    const auto *di0_5 = buffer.data(di0 + 5);
    const auto *di0_6 = buffer.data(di0 + 6);
    const auto *di0_7 = buffer.data(di0 + 7);
    const auto *di0_8 = buffer.data(di0 + 8);
    const auto *di0_9 = buffer.data(di0 + 9);
    const auto *di0_10 = buffer.data(di0 + 10);
    const auto *di0_11 = buffer.data(di0 + 11);
    const auto *di0_12 = buffer.data(di0 + 12);
    const auto *di0_13 = buffer.data(di0 + 13);
    const auto *di0_14 = buffer.data(di0 + 14);
    const auto *di0_15 = buffer.data(di0 + 15);
    const auto *di0_16 = buffer.data(di0 + 16);
    const auto *di0_17 = buffer.data(di0 + 17);
    const auto *di0_18 = buffer.data(di0 + 18);
    const auto *di0_19 = buffer.data(di0 + 19);
    const auto *di0_20 = buffer.data(di0 + 20);
    const auto *di0_21 = buffer.data(di0 + 21);
    const auto *di0_22 = buffer.data(di0 + 22);
    const auto *di0_23 = buffer.data(di0 + 23);
    const auto *di0_24 = buffer.data(di0 + 24);
    const auto *di0_25 = buffer.data(di0 + 25);
    const auto *di0_26 = buffer.data(di0 + 26);
    const auto *di0_27 = buffer.data(di0 + 27);
    const auto *di0_28 = buffer.data(di0 + 28);
    const auto *di0_29 = buffer.data(di0 + 29);
    const auto *di0_30 = buffer.data(di0 + 30);
    const auto *di0_31 = buffer.data(di0 + 31);
    const auto *di0_32 = buffer.data(di0 + 32);
    const auto *di0_33 = buffer.data(di0 + 33);
    const auto *di0_34 = buffer.data(di0 + 34);
    const auto *di0_35 = buffer.data(di0 + 35);
    const auto *di0_36 = buffer.data(di0 + 36);
    const auto *di0_37 = buffer.data(di0 + 37);
    const auto *di0_38 = buffer.data(di0 + 38);
    const auto *di0_39 = buffer.data(di0 + 39);
    const auto *di0_40 = buffer.data(di0 + 40);
    const auto *di0_41 = buffer.data(di0 + 41);
    const auto *di0_42 = buffer.data(di0 + 42);
    const auto *di0_43 = buffer.data(di0 + 43);
    const auto *di0_44 = buffer.data(di0 + 44);
    const auto *di0_45 = buffer.data(di0 + 45);
    const auto *di0_46 = buffer.data(di0 + 46);
    const auto *di0_47 = buffer.data(di0 + 47);
    const auto *di0_48 = buffer.data(di0 + 48);
    const auto *di0_49 = buffer.data(di0 + 49);
    const auto *di0_50 = buffer.data(di0 + 50);
    const auto *di0_51 = buffer.data(di0 + 51);
    const auto *di0_52 = buffer.data(di0 + 52);
    const auto *di0_53 = buffer.data(di0 + 53);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_1 = buffer.data(di1 + 1);
    const auto *di1_2 = buffer.data(di1 + 2);
    const auto *di1_3 = buffer.data(di1 + 3);
    const auto *di1_4 = buffer.data(di1 + 4);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_6 = buffer.data(di1 + 6);
    const auto *di1_7 = buffer.data(di1 + 7);
    const auto *di1_8 = buffer.data(di1 + 8);
    const auto *di1_9 = buffer.data(di1 + 9);
    const auto *di1_10 = buffer.data(di1 + 10);
    const auto *di1_11 = buffer.data(di1 + 11);
    const auto *di1_12 = buffer.data(di1 + 12);
    const auto *di1_14 = buffer.data(di1 + 14);
    const auto *di1_15 = buffer.data(di1 + 15);
    const auto *di1_16 = buffer.data(di1 + 16);
    const auto *di1_17 = buffer.data(di1 + 17);
    const auto *di1_18 = buffer.data(di1 + 18);
    const auto *di1_23 = buffer.data(di1 + 23);
    const auto *di1_25 = buffer.data(di1 + 25);
    const auto *di1_26 = buffer.data(di1 + 26);
    const auto *di1_27 = buffer.data(di1 + 27);
    const auto *di1_28 = buffer.data(di1 + 28);
    const auto *di1_29 = buffer.data(di1 + 29);
    const auto *di1_30 = buffer.data(di1 + 30);
    const auto *di1_31 = buffer.data(di1 + 31);
    const auto *di1_32 = buffer.data(di1 + 32);
    const auto *di1_33 = buffer.data(di1 + 33);
    const auto *di1_34 = buffer.data(di1 + 34);
    const auto *di1_35 = buffer.data(di1 + 35);
    const auto *di1_36 = buffer.data(di1 + 36);
    const auto *di1_37 = buffer.data(di1 + 37);
    const auto *di1_38 = buffer.data(di1 + 38);
    const auto *di1_39 = buffer.data(di1 + 39);
    const auto *di1_40 = buffer.data(di1 + 40);
    const auto *di1_41 = buffer.data(di1 + 41);
    const auto *di1_47 = buffer.data(di1 + 47);
    const auto *di1_49 = buffer.data(di1 + 49);
    const auto *di1_50 = buffer.data(di1 + 50);
    const auto *di1_51 = buffer.data(di1 + 51);
    const auto *di1_52 = buffer.data(di1 + 52);
    const auto *di1_53 = buffer.data(di1 + 53);
    const auto *di1_54 = buffer.data(di1 + 54);
    const auto *di1_55 = buffer.data(di1 + 55);
    const auto *di1_56 = buffer.data(di1 + 56);
    const auto *di1_57 = buffer.data(di1 + 57);
    const auto *di1_58 = buffer.data(di1 + 58);
    const auto *di1_59 = buffer.data(di1 + 59);
    const auto *di1_60 = buffer.data(di1 + 60);
    const auto *di1_61 = buffer.data(di1 + 61);
    const auto *di1_62 = buffer.data(di1 + 62);
    const auto *di1_63 = buffer.data(di1 + 63);
    const auto *di1_64 = buffer.data(di1 + 64);
    const auto *di1_65 = buffer.data(di1 + 65);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pk_0, di0_0, di0_1, di1_0, \
                         di1_1, dk_0, dk_1, dk_2, dk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_3 * di0_0[k]
                 - f_4 * di1_0[k]
                 + pb_y[k] * dk_1[k];

        t_2[k] = f_3 * di0_0[k]
                 - f_4 * di1_0[k]
                 + pb_z[k] * dk_2[k];

        t_3[k] = f_5 * di0_1[k]
                 - f_6 * di1_1[k]
                 + pb_y[k] * dk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, di0_2, di0_3, di0_4, di1_2, di1_3, \
                         di1_4, dk_4, dk_5, dk_6, dk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * di0_2[k]
                 - f_6 * di1_2[k]
                 + pb_z[k] * dk_4[k];

        t_5[k] = f_7 * di0_3[k]
                 - f_8 * di1_3[k]
                 + pb_y[k] * dk_5[k];

        t_6[k] = f_3 * di0_4[k]
                 - f_4 * di1_4[k]
                 + pb_y[k] * dk_6[k];

        t_7[k] = f_7 * di0_4[k]
                 - f_8 * di1_4[k]
                 + pb_z[k] * dk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, di0_5, di0_6, di0_7, di1_5, di1_6, \
                         di1_7, dk_8, dk_9, dk_10, dk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * di0_5[k]
                 - f_10 * di1_5[k]
                 + pb_y[k] * dk_8[k];

        t_9[k] = f_5 * di0_6[k]
                 - f_6 * di1_6[k]
                 + pb_y[k] * dk_9[k];

        t_10[k] = f_3 * di0_7[k]
                  - f_4 * di1_7[k]
                  + pb_y[k] * dk_10[k];

        t_11[k] = f_9 * di0_7[k]
                  - f_10 * di1_7[k]
                  + pb_z[k] * dk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, di0_8, di0_9, di0_10, di1_8, di1_9, di1_10, \
                         dk_12, dk_13, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * di0_8[k]
                  - f_12 * di1_8[k]
                  + pb_y[k] * dk_12[k];

        t_13[k] = f_7 * di0_9[k]
                  - f_8 * di1_9[k]
                  + pb_y[k] * dk_13[k];

        t_14[k] = f_5 * di0_10[k]
                  - f_6 * di1_10[k]
                  + pb_y[k] * dk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pk_9, pk_10, di0_11, \
                         di1_11, dk_15, dk_16, dk_17, dk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * di0_11[k]
                  - f_4 * di1_11[k]
                  + pb_y[k] * dk_15[k];

        t_16[k] = f_11 * di0_11[k]
                  - f_12 * di1_11[k]
                  + pb_z[k] * dk_16[k];

        t_17[k] = f_0 * pk_9[k]
                  + pb_x[k] * dk_17[k];

        t_18[k] = f_0 * pk_10[k]
                  + pb_x[k] * dk_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, di0_12, di0_13, di0_14, di1_12, di1_14, \
                         di1_15, dk_17, dk_18, dk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * di0_12[k]
                  - f_2 * di1_12[k]
                  + pb_y[k] * dk_17[k];

        t_20[k] = f_11 * di0_13[k]
                  - f_12 * di1_14[k]
                  + pb_y[k] * dk_18[k];

        t_21[k] = f_9 * di0_14[k]
                  - f_10 * di1_15[k]
                  + pb_y[k] * dk_19[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_y, pb_z, di0_15, di0_16, di0_17, di1_16, \
                         di1_17, di1_18, dk_20, dk_21, dk_22, dk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * di0_15[k]
                  - f_8 * di1_16[k]
                  + pb_y[k] * dk_20[k];

        t_23[k] = f_5 * di0_16[k]
                  - f_6 * di1_17[k]
                  + pb_y[k] * dk_21[k];

        t_24[k] = f_3 * di0_17[k]
                  - f_4 * di1_18[k]
                  + pb_y[k] * dk_22[k];

        t_25[k] = f_1 * di0_17[k]
                  - f_2 * di1_18[k]
                  + pb_z[k] * dk_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pb_y, pk_0, pk_12, pk_13, pk_14, pl_1, \
                         pl_2, pl_3, dk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_13 * pk_0[k]
                  + pb_y[k] * dk_24[k];

        t_27[k] = f_14 * pk_12[k]
                  + pa_x[k] * pl_1[k];

        t_28[k] = f_15 * pk_13[k]
                  + pa_x[k] * pl_2[k];

        t_29[k] = f_16 * pk_14[k]
                  + pa_x[k] * pl_3[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_z, pb_x, pk_15, pk_16, pk_17, pl_0, \
                         pl_4, pl_5, dk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_17 * pk_15[k]
                  + pa_x[k] * pl_4[k];

        t_31[k] = f_0 * pk_16[k]
                  + pa_x[k] * pl_5[k];

        t_32[k] = f_13 * pk_17[k]
                  + pb_x[k] * dk_29[k];

        t_33[k] = pa_z[k] * pl_0[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_z, pk_0, pk_26, pk_28, pk_31, pl_7, \
                         pl_8, pl_10, dk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * pk_0[k]
                  + pb_z[k] * dk_30[k];

        t_35[k] = f_14 * pk_26[k]
                  + pa_x[k] * pl_7[k];

        t_36[k] = f_15 * pk_28[k]
                  + pa_x[k] * pl_8[k];

        t_37[k] = f_16 * pk_31[k]
                  + pa_x[k] * pl_10[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_x, pk_35, pk_36, pk_44, pl_13, \
                         pl_17, di0_18, di1_23, dk_36, dk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_17 * pk_35[k]
                  + pa_x[k] * pl_13[k];

        t_39[k] = f_0 * pk_36[k]
                  + pa_x[k] * pl_17[k];

        t_40[k] = f_13 * pk_44[k]
                  + pb_x[k] * dk_36[k];

        t_41[k] = f_1 * di0_18[k]
                  - f_2 * di1_23[k]
                  + pb_x[k] * dk_37[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, pk_11, di0_19, di0_20, di1_25, di1_26, \
                         dk_37, dk_38, dk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_0 * pk_11[k]
                  + pb_y[k] * dk_37[k];

        t_43[k] = f_11 * di0_19[k]
                  - f_12 * di1_25[k]
                  + pb_x[k] * dk_38[k];

        t_44[k] = f_11 * di0_20[k]
                  - f_12 * di1_26[k]
                  + pb_x[k] * dk_39[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, di0_21, di0_22, di0_23, di1_27, di1_28, \
                         di1_29, dk_40, dk_41, dk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * di0_21[k]
                  - f_10 * di1_27[k]
                  + pb_x[k] * dk_40[k];

        t_46[k] = f_9 * di0_22[k]
                  - f_10 * di1_28[k]
                  + pb_x[k] * dk_41[k];

        t_47[k] = f_7 * di0_23[k]
                  - f_8 * di1_29[k]
                  + pb_x[k] * dk_42[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, di0_24, di0_25, di0_26, di1_30, di1_31, \
                         di1_32, dk_43, dk_44, dk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_7 * di0_24[k]
                  - f_8 * di1_30[k]
                  + pb_x[k] * dk_43[k];

        t_49[k] = f_7 * di0_25[k]
                  - f_8 * di1_31[k]
                  + pb_x[k] * dk_44[k];

        t_50[k] = f_5 * di0_26[k]
                  - f_6 * di1_32[k]
                  + pb_x[k] * dk_45[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, di0_27, di0_28, di0_29, di1_33, di1_34, \
                         di1_35, dk_46, dk_47, dk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * di0_27[k]
                  - f_6 * di1_33[k]
                  + pb_x[k] * dk_46[k];

        t_52[k] = f_5 * di0_28[k]
                  - f_6 * di1_34[k]
                  + pb_x[k] * dk_47[k];

        t_53[k] = f_5 * di0_29[k]
                  - f_6 * di1_35[k]
                  + pb_x[k] * dk_48[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, di0_30, di0_32, di0_33, di1_36, di1_38, \
                         di1_39, dk_49, dk_50, dk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * di0_30[k]
                  - f_4 * di1_36[k]
                  + pb_x[k] * dk_49[k];

        t_55[k] = f_3 * di0_32[k]
                  - f_4 * di1_38[k]
                  + pb_x[k] * dk_50[k];

        t_56[k] = f_3 * di0_33[k]
                  - f_4 * di1_39[k]
                  + pb_x[k] * dk_51[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, pb_y, pk_17, di0_30, di0_34, di0_35, di1_36, \
                         di1_40, di1_41, dk_52, dk_53, dk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * di0_34[k]
                  - f_4 * di1_40[k]
                  + pb_x[k] * dk_52[k];

        t_58[k] = f_3 * di0_35[k]
                  - f_4 * di1_41[k]
                  + pb_x[k] * dk_53[k];

        t_59[k] = f_0 * pk_17[k]
                  + f_1 * di0_30[k]
                  - f_2 * di1_36[k]
                  + pb_y[k] * dk_54[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_z, di0_30, di0_31, di0_32, di1_36, di1_37, \
                         di1_38, dk_55, dk_56, dk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * di0_30[k]
                  - f_4 * di1_36[k]
                  + pb_z[k] * dk_55[k];

        t_61[k] = f_5 * di0_31[k]
                  - f_6 * di1_37[k]
                  + pb_z[k] * dk_56[k];

        t_62[k] = f_7 * di0_32[k]
                  - f_8 * di1_38[k]
                  + pb_z[k] * dk_57[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_y, pb_z, pk_23, di0_33, di0_34, di0_35, \
                         di1_39, di1_40, di1_41, dk_58, dk_59, dk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * di0_33[k]
                  - f_10 * di1_39[k]
                  + pb_z[k] * dk_58[k];

        t_64[k] = f_11 * di0_34[k]
                  - f_12 * di1_40[k]
                  + pb_z[k] * dk_59[k];

        t_65[k] = f_0 * pk_23[k]
                  + pb_y[k] * dk_61[k];

        t_66[k] = f_1 * di0_35[k]
                  - f_2 * di1_41[k]
                  + pb_z[k] * dk_61[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pa_y, pk_27, pk_29, pk_30, pk_32, \
                         pk_33, pl_9, pl_11, pl_12, pl_14, pl_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * pk_27[k]
                  + pa_y[k] * pl_9[k];

        t_68[k] = f_17 * pk_29[k]
                  + pa_y[k] * pl_11[k];

        t_69[k] = f_0 * pk_30[k]
                  + pa_y[k] * pl_12[k];

        t_70[k] = f_16 * pk_32[k]
                  + pa_y[k] * pl_14[k];

        t_71[k] = f_17 * pk_33[k]
                  + pa_y[k] * pl_15[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pa_z, pb_z, pk_17, pk_34, pk_39, pl_6, \
                         pl_16, pl_18, dk_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * pk_34[k]
                  + pa_y[k] * pl_16[k];

        t_73[k] = pa_z[k] * pl_6[k];

        t_74[k] = f_13 * pk_17[k]
                  + pb_z[k] * dk_62[k];

        t_75[k] = f_14 * pk_39[k]
                  + pa_y[k] * pl_18[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pk_40, pk_41, pk_42, pk_43, pl_19, \
                         pl_20, pl_21, pl_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_15 * pk_40[k]
                  + pa_y[k] * pl_19[k];

        t_77[k] = f_16 * pk_41[k]
                  + pa_y[k] * pl_20[k];

        t_78[k] = f_17 * pk_42[k]
                  + pa_y[k] * pl_21[k];

        t_79[k] = f_0 * pk_43[k]
                  + pa_y[k] * pl_22[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pb_x, pb_y, pb_z, pk_24, pk_44, pl_23, \
                         di0_36, di1_47, dk_69, dk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_13 * pk_44[k]
                  + pb_y[k] * dk_69[k];

        t_81[k] = pa_y[k] * pl_23[k];

        t_82[k] = f_1 * di0_36[k]
                  - f_2 * di1_47[k]
                  + pb_x[k] * dk_70[k];

        t_83[k] = f_0 * pk_24[k]
                  + pb_z[k] * dk_70[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, di0_37, di0_38, di0_39, di1_49, di1_50, \
                         di1_51, dk_72, dk_73, dk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * di0_37[k]
                  - f_12 * di1_49[k]
                  + pb_x[k] * dk_72[k];

        t_85[k] = f_11 * di0_38[k]
                  - f_12 * di1_50[k]
                  + pb_x[k] * dk_73[k];

        t_86[k] = f_9 * di0_39[k]
                  - f_10 * di1_51[k]
                  + pb_x[k] * dk_74[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, di0_40, di0_41, di0_42, di1_52, di1_53, \
                         di1_54, dk_75, dk_76, dk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_9 * di0_40[k]
                  - f_10 * di1_52[k]
                  + pb_x[k] * dk_75[k];

        t_88[k] = f_7 * di0_41[k]
                  - f_8 * di1_53[k]
                  + pb_x[k] * dk_76[k];

        t_89[k] = f_7 * di0_42[k]
                  - f_8 * di1_54[k]
                  + pb_x[k] * dk_77[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_x, di0_43, di0_44, di0_45, di1_55, di1_56, \
                         di1_57, dk_78, dk_79, dk_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_7 * di0_43[k]
                  - f_8 * di1_55[k]
                  + pb_x[k] * dk_78[k];

        t_91[k] = f_5 * di0_44[k]
                  - f_6 * di1_56[k]
                  + pb_x[k] * dk_79[k];

        t_92[k] = f_5 * di0_45[k]
                  - f_6 * di1_57[k]
                  + pb_x[k] * dk_80[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, di0_46, di0_47, di0_48, di1_58, di1_59, \
                         di1_60, dk_81, dk_82, dk_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_5 * di0_46[k]
                  - f_6 * di1_58[k]
                  + pb_x[k] * dk_81[k];

        t_94[k] = f_5 * di0_47[k]
                  - f_6 * di1_59[k]
                  + pb_x[k] * dk_82[k];

        t_95[k] = f_3 * di0_48[k]
                  - f_4 * di1_60[k]
                  + pb_x[k] * dk_83[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_x, di0_49, di0_50, di0_51, di1_61, di1_62, \
                         di1_63, dk_84, dk_85, dk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_3 * di0_49[k]
                  - f_4 * di1_61[k]
                  + pb_x[k] * dk_84[k];

        t_97[k] = f_3 * di0_50[k]
                  - f_4 * di1_62[k]
                  + pb_x[k] * dk_85[k];

        t_98[k] = f_3 * di0_51[k]
                  - f_4 * di1_63[k]
                  + pb_x[k] * dk_86[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_y, pb_z, pk_37, di0_48, di0_53, di1_60, \
                         di1_65, dk_87, dk_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_3 * di0_53[k]
                  - f_4 * di1_65[k]
                  + pb_x[k] * dk_87[k];

        t_100[k] = f_1 * di0_48[k]
                   - f_2 * di1_60[k]
                   + pb_y[k] * dk_88[k];

        t_101[k] = f_0 * pk_37[k]
                   + pb_z[k] * dk_88[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_y, di0_49, di0_50, di0_51, di1_61, di1_62, \
                         di1_63, dk_90, dk_91, dk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_11 * di0_49[k]
                   - f_12 * di1_61[k]
                   + pb_y[k] * dk_90[k];

        t_103[k] = f_9 * di0_50[k]
                   - f_10 * di1_62[k]
                   + pb_y[k] * dk_91[k];

        t_104[k] = f_7 * di0_51[k]
                   - f_8 * di1_63[k]
                   + pb_y[k] * dk_92[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pb_y, pb_z, pk_44, di0_52, di0_53, di1_64, \
                         di1_65, dk_93, dk_94, dk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_5 * di0_52[k]
                   - f_6 * di1_64[k]
                   + pb_y[k] * dk_93[k];

        t_106[k] = f_3 * di0_53[k]
                   - f_4 * di1_65[k]
                   + pb_y[k] * dk_94[k];

        t_107[k] = f_0 * pk_44[k]
                   + f_1 * di0_53[k]
                   - f_2 * di1_65[k]
                   + pb_z[k] * dk_95[k];
    }
}

auto
compute_prim_dl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_1 = buffer.data(di1 + 1);
    const auto *di1_2 = buffer.data(di1 + 2);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_1, di1_0, di1_1, dk_0, \
                         dk_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_1[k]
                 - f_2 * di1_1[k]
                 + pb_y[k] * dk_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_2, di1_2, dk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_2[k]
                 - f_2 * di1_2[k]
                 + pb_z[k] * dk_2[k];
    }
}

auto
compute_prim_dl_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_37 = buffer.data(di1 + 37);
    const auto *di1_65 = buffer.data(di1 + 65);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_95 = buffer.data(dk + 95);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_1, di1_0, di1_37, dk_0, \
                         dk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_1[k]
                 - f_2 * di1_37[k]
                 + pb_y[k] * dk_54[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_2, di1_65, dk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_2[k]
                 - f_2 * di1_65[k]
                 + pb_z[k] * dk_95[k];
    }
}

auto
compute_prim_dl_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_14 = buffer.data(pk + 14);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);
    const auto *di0_3 = buffer.data(di0 + 3);
    const auto *di0_4 = buffer.data(di0 + 4);
    const auto *di0_5 = buffer.data(di0 + 5);
    const auto *di0_7 = buffer.data(di0 + 7);
    const auto *di0_8 = buffer.data(di0 + 8);
    const auto *di0_9 = buffer.data(di0 + 9);
    const auto *di0_11 = buffer.data(di0 + 11);
    const auto *di0_12 = buffer.data(di0 + 12);
    const auto *di0_13 = buffer.data(di0 + 13);
    const auto *di0_14 = buffer.data(di0 + 14);
    const auto *di0_16 = buffer.data(di0 + 16);
    const auto *di0_17 = buffer.data(di0 + 17);
    const auto *di0_18 = buffer.data(di0 + 18);
    const auto *di0_19 = buffer.data(di0 + 19);
    const auto *di0_20 = buffer.data(di0 + 20);
    const auto *di0_22 = buffer.data(di0 + 22);
    const auto *di0_24 = buffer.data(di0 + 24);
    const auto *di0_25 = buffer.data(di0 + 25);
    const auto *di0_26 = buffer.data(di0 + 26);
    const auto *di0_28 = buffer.data(di0 + 28);
    const auto *di0_29 = buffer.data(di0 + 29);
    const auto *di0_31 = buffer.data(di0 + 31);
    const auto *di0_32 = buffer.data(di0 + 32);
    const auto *di0_33 = buffer.data(di0 + 33);
    const auto *di0_34 = buffer.data(di0 + 34);
    const auto *di0_35 = buffer.data(di0 + 35);
    const auto *di0_36 = buffer.data(di0 + 36);
    const auto *di0_37 = buffer.data(di0 + 37);
    const auto *di0_38 = buffer.data(di0 + 38);
    const auto *di0_39 = buffer.data(di0 + 39);
    const auto *di0_40 = buffer.data(di0 + 40);
    const auto *di0_41 = buffer.data(di0 + 41);
    const auto *di0_42 = buffer.data(di0 + 42);
    const auto *di0_45 = buffer.data(di0 + 45);
    const auto *di0_47 = buffer.data(di0 + 47);
    const auto *di0_48 = buffer.data(di0 + 48);
    const auto *di0_49 = buffer.data(di0 + 49);
    const auto *di0_51 = buffer.data(di0 + 51);
    const auto *di0_52 = buffer.data(di0 + 52);
    const auto *di0_53 = buffer.data(di0 + 53);
    const auto *di0_55 = buffer.data(di0 + 55);
    const auto *di0_56 = buffer.data(di0 + 56);
    const auto *di0_57 = buffer.data(di0 + 57);
    const auto *di0_58 = buffer.data(di0 + 58);
    const auto *di0_59 = buffer.data(di0 + 59);
    const auto *di0_60 = buffer.data(di0 + 60);
    const auto *di0_61 = buffer.data(di0 + 61);
    const auto *di0_62 = buffer.data(di0 + 62);
    const auto *di0_63 = buffer.data(di0 + 63);
    const auto *di0_64 = buffer.data(di0 + 64);
    const auto *di0_65 = buffer.data(di0 + 65);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_1 = buffer.data(di1 + 1);
    const auto *di1_2 = buffer.data(di1 + 2);
    const auto *di1_3 = buffer.data(di1 + 3);
    const auto *di1_4 = buffer.data(di1 + 4);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_6 = buffer.data(di1 + 6);
    const auto *di1_7 = buffer.data(di1 + 7);
    const auto *di1_8 = buffer.data(di1 + 8);
    const auto *di1_9 = buffer.data(di1 + 9);
    const auto *di1_10 = buffer.data(di1 + 10);
    const auto *di1_11 = buffer.data(di1 + 11);
    const auto *di1_12 = buffer.data(di1 + 12);
    const auto *di1_14 = buffer.data(di1 + 14);
    const auto *di1_15 = buffer.data(di1 + 15);
    const auto *di1_16 = buffer.data(di1 + 16);
    const auto *di1_17 = buffer.data(di1 + 17);
    const auto *di1_18 = buffer.data(di1 + 18);
    const auto *di1_20 = buffer.data(di1 + 20);
    const auto *di1_22 = buffer.data(di1 + 22);
    const auto *di1_23 = buffer.data(di1 + 23);
    const auto *di1_24 = buffer.data(di1 + 24);
    const auto *di1_25 = buffer.data(di1 + 25);
    const auto *di1_26 = buffer.data(di1 + 26);
    const auto *di1_27 = buffer.data(di1 + 27);
    const auto *di1_28 = buffer.data(di1 + 28);
    const auto *di1_29 = buffer.data(di1 + 29);
    const auto *di1_30 = buffer.data(di1 + 30);
    const auto *di1_31 = buffer.data(di1 + 31);
    const auto *di1_32 = buffer.data(di1 + 32);
    const auto *di1_33 = buffer.data(di1 + 33);
    const auto *di1_34 = buffer.data(di1 + 34);
    const auto *di1_35 = buffer.data(di1 + 35);
    const auto *di1_36 = buffer.data(di1 + 36);
    const auto *di1_37 = buffer.data(di1 + 37);
    const auto *di1_38 = buffer.data(di1 + 38);
    const auto *di1_41 = buffer.data(di1 + 41);
    const auto *di1_43 = buffer.data(di1 + 43);
    const auto *di1_44 = buffer.data(di1 + 44);
    const auto *di1_45 = buffer.data(di1 + 45);
    const auto *di1_46 = buffer.data(di1 + 46);
    const auto *di1_47 = buffer.data(di1 + 47);
    const auto *di1_48 = buffer.data(di1 + 48);
    const auto *di1_49 = buffer.data(di1 + 49);
    const auto *di1_50 = buffer.data(di1 + 50);
    const auto *di1_51 = buffer.data(di1 + 51);
    const auto *di1_52 = buffer.data(di1 + 52);
    const auto *di1_53 = buffer.data(di1 + 53);
    const auto *di1_54 = buffer.data(di1 + 54);
    const auto *di1_55 = buffer.data(di1 + 55);
    const auto *di1_56 = buffer.data(di1 + 56);
    const auto *di1_57 = buffer.data(di1 + 57);
    const auto *di1_58 = buffer.data(di1 + 58);
    const auto *di1_59 = buffer.data(di1 + 59);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pk_0, di0_0, di0_1, di1_0, \
                         di1_1, dk_0, dk_1, dk_2, dk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_3 * di0_0[k]
                 - f_4 * di1_0[k]
                 + pb_y[k] * dk_1[k];

        t_2[k] = f_3 * di0_0[k]
                 - f_4 * di1_0[k]
                 + pb_z[k] * dk_2[k];

        t_3[k] = f_5 * di0_1[k]
                 - f_6 * di1_1[k]
                 + pb_y[k] * dk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, di0_2, di0_3, di0_4, di1_2, di1_3, \
                         di1_4, dk_4, dk_5, dk_6, dk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * di0_2[k]
                 - f_6 * di1_2[k]
                 + pb_z[k] * dk_4[k];

        t_5[k] = f_7 * di0_3[k]
                 - f_8 * di1_3[k]
                 + pb_y[k] * dk_5[k];

        t_6[k] = f_3 * di0_4[k]
                 - f_4 * di1_4[k]
                 + pb_y[k] * dk_6[k];

        t_7[k] = f_7 * di0_4[k]
                 - f_8 * di1_4[k]
                 + pb_z[k] * dk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, di0_5, di0_7, di0_8, di1_5, di1_6, \
                         di1_7, dk_8, dk_9, dk_10, dk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * di0_5[k]
                 - f_10 * di1_5[k]
                 + pb_y[k] * dk_8[k];

        t_9[k] = f_5 * di0_7[k]
                 - f_6 * di1_6[k]
                 + pb_y[k] * dk_9[k];

        t_10[k] = f_3 * di0_8[k]
                  - f_4 * di1_7[k]
                  + pb_y[k] * dk_10[k];

        t_11[k] = f_9 * di0_8[k]
                  - f_10 * di1_7[k]
                  + pb_z[k] * dk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, di0_9, di0_11, di0_12, di1_8, di1_9, di1_10, \
                         dk_12, dk_13, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_11 * di0_9[k]
                  - f_12 * di1_8[k]
                  + pb_y[k] * dk_12[k];

        t_13[k] = f_7 * di0_11[k]
                  - f_8 * di1_9[k]
                  + pb_y[k] * dk_13[k];

        t_14[k] = f_5 * di0_12[k]
                  - f_6 * di1_10[k]
                  + pb_y[k] * dk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, di0_13, di0_14, di0_16, di1_11, \
                         di1_12, di1_14, dk_15, dk_16, dk_17, dk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * di0_13[k]
                  - f_4 * di1_11[k]
                  + pb_y[k] * dk_15[k];

        t_16[k] = f_11 * di0_13[k]
                  - f_12 * di1_11[k]
                  + pb_z[k] * dk_16[k];

        t_17[k] = f_1 * di0_14[k]
                  - f_2 * di1_12[k]
                  + pb_y[k] * dk_17[k];

        t_18[k] = f_11 * di0_16[k]
                  - f_12 * di1_14[k]
                  + pb_y[k] * dk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, di0_17, di0_18, di0_19, di1_15, di1_16, \
                         di1_17, dk_19, dk_20, dk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_9 * di0_17[k]
                  - f_10 * di1_15[k]
                  + pb_y[k] * dk_19[k];

        t_20[k] = f_7 * di0_18[k]
                  - f_8 * di1_16[k]
                  + pb_y[k] * dk_20[k];

        t_21[k] = f_5 * di0_19[k]
                  - f_6 * di1_17[k]
                  + pb_y[k] * dk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, pb_y, pb_z, di0_20, di0_22, di1_18, di1_20, \
                         dk_22, dk_23, dk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * di0_20[k]
                  - f_4 * di1_18[k]
                  + pb_y[k] * dk_22[k];

        t_23[k] = f_1 * di0_20[k]
                  - f_2 * di1_18[k]
                  + pb_z[k] * dk_23[k];

        t_24[k] = f_1 * di0_22[k]
                  - f_2 * di1_20[k]
                  + pb_x[k] * dk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, di0_24, di0_25, di0_26, di1_22, di1_23, \
                         di1_24, dk_25, dk_26, dk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * di0_24[k]
                  - f_12 * di1_22[k]
                  + pb_x[k] * dk_25[k];

        t_26[k] = f_11 * di0_25[k]
                  - f_12 * di1_23[k]
                  + pb_x[k] * dk_26[k];

        t_27[k] = f_9 * di0_26[k]
                  - f_10 * di1_24[k]
                  + pb_x[k] * dk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, di0_28, di0_29, di0_31, di1_25, di1_26, \
                         di1_27, dk_28, dk_29, dk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_9 * di0_28[k]
                  - f_10 * di1_25[k]
                  + pb_x[k] * dk_28[k];

        t_29[k] = f_7 * di0_29[k]
                  - f_8 * di1_26[k]
                  + pb_x[k] * dk_29[k];

        t_30[k] = f_7 * di0_31[k]
                  - f_8 * di1_27[k]
                  + pb_x[k] * dk_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, di0_32, di0_33, di0_34, di1_28, di1_29, \
                         di1_30, dk_31, dk_32, dk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * di0_32[k]
                  - f_8 * di1_28[k]
                  + pb_x[k] * dk_31[k];

        t_32[k] = f_5 * di0_33[k]
                  - f_6 * di1_29[k]
                  + pb_x[k] * dk_32[k];

        t_33[k] = f_5 * di0_34[k]
                  - f_6 * di1_30[k]
                  + pb_x[k] * dk_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, di0_35, di0_36, di0_37, di1_31, di1_32, \
                         di1_33, dk_34, dk_35, dk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * di0_35[k]
                  - f_6 * di1_31[k]
                  + pb_x[k] * dk_34[k];

        t_35[k] = f_5 * di0_36[k]
                  - f_6 * di1_32[k]
                  + pb_x[k] * dk_35[k];

        t_36[k] = f_3 * di0_37[k]
                  - f_4 * di1_33[k]
                  + pb_x[k] * dk_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, di0_39, di0_40, di0_41, di1_35, di1_36, \
                         di1_37, dk_37, dk_38, dk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * di0_39[k]
                  - f_4 * di1_35[k]
                  + pb_x[k] * dk_37[k];

        t_38[k] = f_3 * di0_40[k]
                  - f_4 * di1_36[k]
                  + pb_x[k] * dk_38[k];

        t_39[k] = f_3 * di0_41[k]
                  - f_4 * di1_37[k]
                  + pb_x[k] * dk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, pb_y, pb_z, pk_5, di0_37, di0_42, di1_33, \
                         di1_38, dk_40, dk_41, dk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_3 * di0_42[k]
                  - f_4 * di1_38[k]
                  + pb_x[k] * dk_40[k];

        t_41[k] = f_0 * pk_5[k]
                  + f_1 * di0_37[k]
                  - f_2 * di1_33[k]
                  + pb_y[k] * dk_41[k];

        t_42[k] = f_3 * di0_37[k]
                  - f_4 * di1_33[k]
                  + pb_z[k] * dk_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_z, di0_38, di0_39, di0_40, di1_34, di1_35, \
                         di1_36, dk_43, dk_44, dk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * di0_38[k]
                  - f_6 * di1_34[k]
                  + pb_z[k] * dk_43[k];

        t_44[k] = f_7 * di0_39[k]
                  - f_8 * di1_35[k]
                  + pb_z[k] * dk_44[k];

        t_45[k] = f_9 * di0_40[k]
                  - f_10 * di1_36[k]
                  + pb_z[k] * dk_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_z, di0_41, di0_42, di0_45, di1_37, di1_38, \
                         di1_41, dk_46, dk_47, dk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_11 * di0_41[k]
                  - f_12 * di1_37[k]
                  + pb_z[k] * dk_46[k];

        t_47[k] = f_1 * di0_42[k]
                  - f_2 * di1_38[k]
                  + pb_z[k] * dk_47[k];

        t_48[k] = f_1 * di0_45[k]
                  - f_2 * di1_41[k]
                  + pb_x[k] * dk_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, di0_47, di0_48, di0_49, di1_43, di1_44, \
                         di1_45, dk_49, dk_50, dk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * di0_47[k]
                  - f_12 * di1_43[k]
                  + pb_x[k] * dk_49[k];

        t_50[k] = f_11 * di0_48[k]
                  - f_12 * di1_44[k]
                  + pb_x[k] * dk_50[k];

        t_51[k] = f_9 * di0_49[k]
                  - f_10 * di1_45[k]
                  + pb_x[k] * dk_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, di0_51, di0_52, di0_53, di1_46, di1_47, \
                         di1_48, dk_52, dk_53, dk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_9 * di0_51[k]
                  - f_10 * di1_46[k]
                  + pb_x[k] * dk_52[k];

        t_53[k] = f_7 * di0_52[k]
                  - f_8 * di1_47[k]
                  + pb_x[k] * dk_53[k];

        t_54[k] = f_7 * di0_53[k]
                  - f_8 * di1_48[k]
                  + pb_x[k] * dk_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, di0_55, di0_56, di0_57, di1_49, di1_50, \
                         di1_51, dk_55, dk_56, dk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_7 * di0_55[k]
                  - f_8 * di1_49[k]
                  + pb_x[k] * dk_55[k];

        t_56[k] = f_5 * di0_56[k]
                  - f_6 * di1_50[k]
                  + pb_x[k] * dk_56[k];

        t_57[k] = f_5 * di0_57[k]
                  - f_6 * di1_51[k]
                  + pb_x[k] * dk_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, di0_58, di0_59, di0_60, di1_52, di1_53, \
                         di1_54, dk_58, dk_59, dk_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_5 * di0_58[k]
                  - f_6 * di1_52[k]
                  + pb_x[k] * dk_58[k];

        t_59[k] = f_5 * di0_59[k]
                  - f_6 * di1_53[k]
                  + pb_x[k] * dk_59[k];

        t_60[k] = f_3 * di0_60[k]
                  - f_4 * di1_54[k]
                  + pb_x[k] * dk_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, di0_61, di0_62, di0_63, di1_55, di1_56, \
                         di1_57, dk_61, dk_62, dk_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * di0_61[k]
                  - f_4 * di1_55[k]
                  + pb_x[k] * dk_61[k];

        t_62[k] = f_3 * di0_62[k]
                  - f_4 * di1_56[k]
                  + pb_x[k] * dk_62[k];

        t_63[k] = f_3 * di0_63[k]
                  - f_4 * di1_57[k]
                  + pb_x[k] * dk_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, di0_60, di0_61, di0_65, di1_54, di1_55, \
                         di1_59, dk_64, dk_65, dk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * di0_65[k]
                  - f_4 * di1_59[k]
                  + pb_x[k] * dk_64[k];

        t_65[k] = f_1 * di0_60[k]
                  - f_2 * di1_54[k]
                  + pb_y[k] * dk_65[k];

        t_66[k] = f_11 * di0_61[k]
                  - f_12 * di1_55[k]
                  + pb_y[k] * dk_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_y, di0_62, di0_63, di0_64, di1_56, di1_57, \
                         di1_58, dk_67, dk_68, dk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_9 * di0_62[k]
                  - f_10 * di1_56[k]
                  + pb_y[k] * dk_67[k];

        t_68[k] = f_7 * di0_63[k]
                  - f_8 * di1_57[k]
                  + pb_y[k] * dk_68[k];

        t_69[k] = f_5 * di0_64[k]
                  - f_6 * di1_58[k]
                  + pb_y[k] * dk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, pb_y, pb_z, pk_14, di0_65, di1_59, dk_70, \
                         dk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * di0_65[k]
                  - f_4 * di1_59[k]
                  + pb_y[k] * dk_70[k];

        t_71[k] = f_0 * pk_14[k]
                  + f_1 * di0_65[k]
                  - f_2 * di1_59[k]
                  + pb_z[k] * dk_71[k];
    }
}

auto
compute_prim_dl_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_25 = buffer.data(di1 + 25);
    const auto *di1_47 = buffer.data(di1 + 47);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_44 = buffer.data(dk + 44);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_1, di1_0, di1_25, dk_0, \
                         dk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_1[k]
                 - f_2 * di1_25[k]
                 + pb_y[k] * dk_24[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_2, di1_47, dk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_2[k]
                 - f_2 * di1_47[k]
                 + pb_z[k] * dk_44[k];
    }
}

auto
compute_prim_dl_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_25 = buffer.data(di0 + 25);
    const auto *di0_47 = buffer.data(di0 + 47);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_24 = buffer.data(di1 + 24);
    const auto *di1_44 = buffer.data(di1 + 44);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_25, di1_0, di1_24, dk_0, \
                         dk_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_25[k]
                 - f_2 * di1_24[k]
                 + pb_y[k] * dk_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_47, di1_44, dk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_47[k]
                 - f_2 * di1_44[k]
                 + pb_z[k] * dk_2[k];
    }
}

auto
compute_prim_dl_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_1 = buffer.data(di0 + 1);
    const auto *di0_2 = buffer.data(di0 + 2);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_14 = buffer.data(di1 + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_1, di1_0, di1_5, dk_0, \
                         dk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_1[k]
                 - f_2 * di1_5[k]
                 + pb_y[k] * dk_5[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_2, di1_14, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_2[k]
                 - f_2 * di1_14[k]
                 + pb_z[k] * dk_14[k];
    }
}

auto
compute_prim_dl_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_5 = buffer.data(di0 + 5);
    const auto *di0_14 = buffer.data(di0 + 14);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_14 = buffer.data(di1 + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_5, di1_0, di1_5, dk_0, \
                         dk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_5[k]
                 - f_2 * di1_5[k]
                 + pb_y[k] * dk_5[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_14, di1_14, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_14[k]
                 - f_2 * di1_14[k]
                 + pb_z[k] * dk_14[k];
    }
}

auto
compute_prim_dl_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pk, const size_t di0, const size_t di1,
                                     const size_t dk, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_5 = buffer.data(di0 + 5);
    const auto *di0_14 = buffer.data(di0 + 14);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_8 = buffer.data(di1 + 8);
    const auto *di1_20 = buffer.data(di1 + 20);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_14 = buffer.data(dk + 14);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_5, di1_0, di1_8, dk_0, \
                         dk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_5[k]
                 - f_2 * di1_8[k]
                 + pb_y[k] * dk_5[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_14, di1_20, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_14[k]
                 - f_2 * di1_20[k]
                 + pb_z[k] * dk_14[k];
    }
}

auto
compute_prim_dl_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pk, const size_t di0, const size_t di1,
                                      const size_t dk, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_8 = buffer.data(di0 + 8);
    const auto *di0_20 = buffer.data(di0 + 20);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_14 = buffer.data(di1 + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_8, di1_0, di1_5, dk_0, \
                         dk_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_8[k]
                 - f_2 * di1_5[k]
                 + pb_y[k] * dk_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_20, di1_14, dk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_20[k]
                 - f_2 * di1_14[k]
                 + pb_z[k] * dk_2[k];
    }
}

auto
compute_prim_dl_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pk, const size_t di0, const size_t di1,
                                      const size_t dk, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);

    const auto *di0_0 = buffer.data(di0 + 0);
    const auto *di0_5 = buffer.data(di0 + 5);
    const auto *di0_14 = buffer.data(di0 + 14);

    const auto *di1_0 = buffer.data(di1 + 0);
    const auto *di1_5 = buffer.data(di1 + 5);
    const auto *di1_14 = buffer.data(di1 + 14);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pk_0, pk_1, di0_0, di0_5, di1_0, di1_5, dk_0, \
                         dk_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k]
                 + f_1 * di0_0[k]
                 - f_2 * di1_0[k]
                 + pb_x[k] * dk_0[k];

        t_1[k] = f_0 * pk_1[k]
                 + f_1 * di0_5[k]
                 - f_2 * di1_5[k]
                 + pb_y[k] * dk_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pk_2, di0_14, di1_14, dk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pk_2[k]
                 + f_1 * di0_14[k]
                 - f_2 * di1_14[k]
                 + pb_z[k] * dk_2[k];
    }
}

}  // namespace simdt2ceri
