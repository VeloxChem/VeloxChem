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


#include "SimdElectronRepulsionVrrRecLD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ld_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 3.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.5 / p;
    const auto f_14 = 2.0 / alpha;
    const auto f_15 = 2.0 * beta / (alpha * p);
    const auto f_16 = 1.5 / alpha;
    const auto f_17 = 1.5 * beta / (alpha * p);
    const auto f_18 = 2.0 / p;
    const auto f_19 = 1.5 / p;

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

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_1 = buffer.data(id0 + 1);
    const auto *id0_2 = buffer.data(id0 + 2);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_7 = buffer.data(id0 + 7);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_13 = buffer.data(id0 + 13);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_15 = buffer.data(id0 + 15);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_19 = buffer.data(id0 + 19);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_25 = buffer.data(id0 + 25);
    const auto *id0_26 = buffer.data(id0 + 26);
    const auto *id0_27 = buffer.data(id0 + 27);
    const auto *id0_28 = buffer.data(id0 + 28);
    const auto *id0_29 = buffer.data(id0 + 29);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_31 = buffer.data(id0 + 31);
    const auto *id0_32 = buffer.data(id0 + 32);
    const auto *id0_33 = buffer.data(id0 + 33);
    const auto *id0_34 = buffer.data(id0 + 34);
    const auto *id0_35 = buffer.data(id0 + 35);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_1 = buffer.data(id1 + 1);
    const auto *id1_2 = buffer.data(id1 + 2);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_4 = buffer.data(id1 + 4);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_7 = buffer.data(id1 + 7);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_9 = buffer.data(id1 + 9);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_11 = buffer.data(id1 + 11);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_13 = buffer.data(id1 + 13);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_15 = buffer.data(id1 + 15);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_19 = buffer.data(id1 + 19);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_22 = buffer.data(id1 + 22);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_25 = buffer.data(id1 + 25);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_27 = buffer.data(id1 + 27);
    const auto *id1_28 = buffer.data(id1 + 28);
    const auto *id1_29 = buffer.data(id1 + 29);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_31 = buffer.data(id1 + 31);
    const auto *id1_32 = buffer.data(id1 + 32);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_34 = buffer.data(id1 + 34);
    const auto *id1_35 = buffer.data(id1 + 35);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_24 = buffer.data(lp + 24);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_33 = buffer.data(lp + 33);
    const auto *lp_34 = buffer.data(lp + 34);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_36 = buffer.data(lp + 36);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_39 = buffer.data(lp + 39);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_48 = buffer.data(lp + 48);
    const auto *lp_49 = buffer.data(lp + 49);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_51 = buffer.data(lp + 51);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_54 = buffer.data(lp + 54);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_57 = buffer.data(lp + 57);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_60 = buffer.data(lp + 60);
    const auto *lp_61 = buffer.data(lp + 61);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_63 = buffer.data(lp + 63);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_65 = buffer.data(lp + 65);
    const auto *lp_66 = buffer.data(lp + 66);
    const auto *lp_67 = buffer.data(lp + 67);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_69 = buffer.data(lp + 69);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_72 = buffer.data(lp + 72);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_75 = buffer.data(lp + 75);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_78 = buffer.data(lp + 78);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_81 = buffer.data(lp + 81);
    const auto *lp_82 = buffer.data(lp + 82);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_84 = buffer.data(lp + 84);
    const auto *lp_85 = buffer.data(lp + 85);
    const auto *lp_86 = buffer.data(lp + 86);
    const auto *lp_87 = buffer.data(lp + 87);
    const auto *lp_88 = buffer.data(lp + 88);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_90 = buffer.data(lp + 90);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_93 = buffer.data(lp + 93);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_96 = buffer.data(lp + 96);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_99 = buffer.data(lp + 99);
    const auto *lp_100 = buffer.data(lp + 100);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, \
                         lp_0, lp_1, lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = pb_y[k] * lp_0[k];

        t_2[k] = pb_z[k] * lp_0[k];

        t_3[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_4[k] = pb_y[k] * lp_2[k];

        t_5[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, kp_1, kp_3, kd_0, \
                         kd_1, kd_2, lp_3, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * kd_0[k];

        t_7[k] = f_3 * kp_3[k]
                 + pb_x[k] * lp_4[k];

        t_8[k] = pb_z[k] * lp_3[k];

        t_9[k] = f_4 * kp_1[k]
                 + pa_y[k] * kd_1[k];

        t_10[k] = pb_z[k] * lp_4[k];

        t_11[k] = pa_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, kp_2, kp_4, \
                         kd_0, kd_1, kd_2, lp_5, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * kd_0[k];

        t_13[k] = pb_y[k] * lp_5[k];

        t_14[k] = f_3 * kp_4[k]
                  + pb_x[k] * lp_6[k];

        t_15[k] = pa_z[k] * kd_1[k];

        t_16[k] = pb_y[k] * lp_6[k];

        t_17[k] = f_4 * kp_2[k]
                  + pa_z[k] * kd_2[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, id0_0, id1_0, kp_5, kd_3, lp_7, \
                         lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * id0_0[k]
                  - f_6 * id1_0[k]
                  + pa_y[k] * kd_3[k];

        t_19[k] = f_7 * kp_5[k]
                  + pb_x[k] * lp_8[k];

        t_20[k] = pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, id0_4, id1_4, kd_6, kd_11, \
                         ls0_1, ls1_1, lp_8, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * id0_4[k]
                  - f_9 * id1_4[k]
                  + pa_x[k] * kd_11[k];

        t_22[k] = pb_z[k] * lp_8[k];

        t_23[k] = f_1 * ls0_1[k]
                  - f_2 * ls1_1[k]
                  + pb_z[k] * lp_9[k];

        t_24[k] = pa_y[k] * kd_6[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, kp_4, kd_4, kd_5, \
                         kd_7, kd_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * kd_4[k];

        t_26[k] = pa_y[k] * kd_7[k];

        t_27[k] = pa_z[k] * kd_5[k];

        t_28[k] = f_10 * kp_4[k]
                  + pb_y[k] * lp_10[k];

        t_29[k] = pa_y[k] * kd_8[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, id0_0, id1_0, kp_9, kd_6, \
                         ls0_2, ls1_2, lp_11, lp_12, lp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_0[k]
                  - f_6 * id1_0[k]
                  + pa_z[k] * kd_6[k];

        t_31[k] = pb_y[k] * lp_11[k];

        t_32[k] = f_7 * kp_9[k]
                  + pb_x[k] * lp_13[k];

        t_33[k] = f_1 * ls0_2[k]
                  - f_2 * ls1_2[k]
                  + pb_y[k] * lp_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, id0_1, id0_6, id1_1, id1_6, kd_9, \
                         kd_16, lp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * lp_13[k];

        t_35[k] = f_8 * id0_6[k]
                  - f_9 * id1_6[k]
                  + pa_x[k] * kd_16[k];

        t_36[k] = f_11 * id0_1[k]
                  - f_12 * id1_1[k]
                  + pa_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, id0_8, id1_8, kp_10, kd_19, \
                         lp_14, lp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * kp_10[k]
                  + pb_x[k] * lp_15[k];

        t_38[k] = pb_z[k] * lp_14[k];

        t_39[k] = f_14 * id0_8[k]
                  - f_15 * id1_8[k]
                  + pa_x[k] * kd_19[k];

        t_40[k] = pb_z[k] * lp_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, kp_12, kd_9, kd_10, \
                         kd_11, ls0_3, ls1_3, lp_16, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_16[k];

        t_42[k] = pa_z[k] * kd_9[k];

        t_43[k] = pa_z[k] * kd_10[k];

        t_44[k] = f_13 * kp_12[k]
                  + pb_x[k] * lp_17[k];

        t_45[k] = pa_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, kp_6, kp_7, kp_13, \
                         kd_12, kd_13, lp_17, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * kp_7[k]
                  + pb_y[k] * lp_17[k];

        t_47[k] = f_4 * kp_6[k]
                  + pa_z[k] * kd_12[k];

        t_48[k] = pa_y[k] * kd_13[k];

        t_49[k] = f_13 * kp_13[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, kp_8, kp_9, kd_14, kd_15, kd_16, \
                         lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * kd_14[k];

        t_51[k] = f_4 * kp_8[k]
                  + pa_y[k] * kd_15[k];

        t_52[k] = f_10 * kp_9[k]
                  + pb_y[k] * lp_19[k];

        t_53[k] = pa_y[k] * kd_16[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, id0_2, id1_2, kp_16, kd_13, \
                         ls0_4, ls1_4, lp_20, lp_21, lp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * id0_2[k]
                  - f_12 * id1_2[k]
                  + pa_z[k] * kd_13[k];

        t_55[k] = pb_y[k] * lp_20[k];

        t_56[k] = f_13 * kp_16[k]
                  + pb_x[k] * lp_22[k];

        t_57[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_21[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, id0_3, id0_11, id1_3, id1_11, \
                         kd_17, kd_25, lp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * lp_22[k];

        t_59[k] = f_14 * id0_11[k]
                  - f_15 * id1_11[k]
                  + pa_x[k] * kd_25[k];

        t_60[k] = f_16 * id0_3[k]
                  - f_17 * id1_3[k]
                  + pa_y[k] * kd_17[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, id0_13, id1_13, kp_17, \
                         kd_28, lp_23, lp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_18 * kp_17[k]
                  + pb_x[k] * lp_24[k];

        t_62[k] = pb_z[k] * lp_23[k];

        t_63[k] = f_16 * id0_13[k]
                  - f_17 * id1_13[k]
                  + pa_x[k] * kd_28[k];

        t_64[k] = pb_z[k] * lp_24[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_z, kp_19, kd_17, kd_18, \
                         kd_19, ls0_5, ls1_5, lp_25, lp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_25[k];

        t_66[k] = pa_z[k] * kd_17[k];

        t_67[k] = pa_z[k] * kd_18[k];

        t_68[k] = f_18 * kp_19[k]
                  + pb_x[k] * lp_26[k];

        t_69[k] = pa_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_y, id0_5, id1_5, kp_11, kp_12, \
                         kd_20, kd_21, lp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_19 * kp_12[k]
                  + pb_y[k] * lp_26[k];

        t_71[k] = f_4 * kp_11[k]
                  + pa_z[k] * kd_20[k];

        t_72[k] = f_5 * id0_5[k]
                  - f_6 * id1_5[k]
                  + pa_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, pb_y, id0_15, id1_15, kp_14, \
                         kp_20, kp_21, kd_31, lp_27, lp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_18 * kp_20[k]
                  + pb_x[k] * lp_27[k];

        t_74[k] = f_18 * kp_21[k]
                  + pb_x[k] * lp_28[k];

        t_75[k] = f_16 * id0_15[k]
                  - f_17 * id1_15[k]
                  + pa_x[k] * kd_31[k];

        t_76[k] = f_4 * kp_14[k]
                  + pb_y[k] * lp_28[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, id0_16, id1_16, kp_22, \
                         kd_22, kd_23, kd_32, lp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_16 * id0_16[k]
                  - f_17 * id1_16[k]
                  + pa_x[k] * kd_32[k];

        t_78[k] = pa_y[k] * kd_22[k];

        t_79[k] = f_18 * kp_22[k]
                  + pb_x[k] * lp_29[k];

        t_80[k] = pa_y[k] * kd_23[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, id0_5, id1_5, kp_15, kp_16, \
                         kd_22, kd_24, kd_25, lp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * kp_15[k]
                  + pa_y[k] * kd_24[k];

        t_82[k] = f_10 * kp_16[k]
                  + pb_y[k] * lp_30[k];

        t_83[k] = pa_y[k] * kd_25[k];

        t_84[k] = f_16 * id0_5[k]
                  - f_17 * id1_5[k]
                  + pa_z[k] * kd_22[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, kp_25, ls0_6, ls1_6, lp_31, \
                         lp_32, lp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * lp_31[k];

        t_86[k] = f_18 * kp_25[k]
                  + pb_x[k] * lp_33[k];

        t_87[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_32[k];

        t_88[k] = pb_y[k] * lp_33[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pa_y, pb_x, id0_7, id0_19, id1_7, id1_19, \
                         kp_26, kd_26, kd_37, lp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_16 * id0_19[k]
                  - f_17 * id1_19[k]
                  + pa_x[k] * kd_37[k];

        t_90[k] = f_14 * id0_7[k]
                  - f_15 * id1_7[k]
                  + pa_y[k] * kd_26[k];

        t_91[k] = f_19 * kp_26[k]
                  + pb_x[k] * lp_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, id0_20, id1_20, kd_40, ls0_7, \
                         ls1_7, lp_34, lp_35, lp_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * lp_34[k];

        t_93[k] = f_11 * id0_20[k]
                  - f_12 * id1_20[k]
                  + pa_x[k] * kd_40[k];

        t_94[k] = pb_z[k] * lp_35[k];

        t_95[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_36[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, kp_19, kp_28, kd_26, \
                         kd_27, kd_28, lp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * kd_26[k];

        t_97[k] = pa_z[k] * kd_27[k];

        t_98[k] = f_19 * kp_28[k]
                  + pb_x[k] * lp_37[k];

        t_99[k] = pa_z[k] * kd_28[k];

        t_100[k] = f_18 * kp_19[k]
                   + pb_y[k] * lp_37[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pa_z, pb_x, id0_9, id1_9, kp_18, \
                         kp_29, kp_30, kd_29, kd_30, lp_38, lp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_4 * kp_18[k]
                   + pa_z[k] * kd_29[k];

        t_102[k] = f_11 * id0_9[k]
                   - f_12 * id1_9[k]
                   + pa_y[k] * kd_30[k];

        t_103[k] = f_19 * kp_29[k]
                   + pb_x[k] * lp_38[k];

        t_104[k] = f_19 * kp_30[k]
                   + pb_x[k] * lp_39[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_y, id0_21, id0_22, id1_21, id1_22, \
                         kp_21, kd_43, kd_44, lp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_11 * id0_21[k]
                   - f_12 * id1_21[k]
                   + pa_x[k] * kd_43[k];

        t_106[k] = f_19 * kp_21[k]
                   + pb_y[k] * lp_39[k];

        t_107[k] = f_11 * id0_22[k]
                   - f_12 * id1_22[k]
                   + pa_x[k] * kd_44[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_x, id0_10, id1_10, kp_31, kp_32, kd_33, \
                         lp_40, lp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * id0_10[k]
                   - f_6 * id1_10[k]
                   + pa_y[k] * kd_33[k];

        t_109[k] = f_19 * kp_31[k]
                   + pb_x[k] * lp_40[k];

        t_110[k] = f_19 * kp_32[k]
                   + pb_x[k] * lp_41[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, pb_y, id0_23, id0_24, id1_23, \
                         id1_24, kp_23, kd_34, kd_46, kd_47, lp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_11 * id0_23[k]
                   - f_12 * id1_23[k]
                   + pa_x[k] * kd_46[k];

        t_112[k] = f_4 * kp_23[k]
                   + pb_y[k] * lp_41[k];

        t_113[k] = f_11 * id0_24[k]
                   - f_12 * id1_24[k]
                   + pa_x[k] * kd_47[k];

        t_114[k] = pa_y[k] * kd_34[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, kp_24, kp_25, \
                         kp_33, kd_35, kd_36, kd_37, lp_42, lp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_19 * kp_33[k]
                   + pb_x[k] * lp_42[k];

        t_116[k] = pa_y[k] * kd_35[k];

        t_117[k] = f_4 * kp_24[k]
                   + pa_y[k] * kd_36[k];

        t_118[k] = f_10 * kp_25[k]
                   + pb_y[k] * lp_43[k];

        t_119[k] = pa_y[k] * kd_37[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, pb_x, pb_y, id0_10, id1_10, kp_36, \
                         kd_34, ls0_8, ls1_8, lp_44, lp_45, lp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_14 * id0_10[k]
                   - f_15 * id1_10[k]
                   + pa_z[k] * kd_34[k];

        t_121[k] = pb_y[k] * lp_44[k];

        t_122[k] = f_19 * kp_36[k]
                   + pb_x[k] * lp_46[k];

        t_123[k] = f_1 * ls0_8[k]
                   - f_2 * ls1_8[k]
                   + pb_y[k] * lp_45[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pa_y, pb_y, id0_12, id0_25, id1_12, \
                         id1_25, kd_38, kd_52, lp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * lp_46[k];

        t_125[k] = f_11 * id0_25[k]
                   - f_12 * id1_25[k]
                   + pa_x[k] * kd_52[k];

        t_126[k] = f_8 * id0_12[k]
                   - f_9 * id1_12[k]
                   + pa_y[k] * kd_38[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pb_x, pb_z, id0_26, id1_26, kp_37, \
                         kd_55, lp_47, lp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_4 * kp_37[k]
                   + pb_x[k] * lp_48[k];

        t_128[k] = pb_z[k] * lp_47[k];

        t_129[k] = f_5 * id0_26[k]
                   - f_6 * id1_26[k]
                   + pa_x[k] * kd_55[k];

        t_130[k] = pb_z[k] * lp_48[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_z, pb_x, pb_z, kp_38, kd_38, \
                         kd_39, kd_40, ls0_9, ls1_9, lp_49, lp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_1 * ls0_9[k]
                   - f_2 * ls1_9[k]
                   + pb_z[k] * lp_49[k];

        t_132[k] = pa_z[k] * kd_38[k];

        t_133[k] = pa_z[k] * kd_39[k];

        t_134[k] = f_4 * kp_38[k]
                   + pb_x[k] * lp_50[k];

        t_135[k] = pa_z[k] * kd_40[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_y, pa_z, pb_y, id0_14, id1_14, kp_27, kp_28, \
                         kd_41, kd_42, lp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_13 * kp_28[k]
                   + pb_y[k] * lp_50[k];

        t_137[k] = f_4 * kp_27[k]
                   + pa_z[k] * kd_41[k];

        t_138[k] = f_16 * id0_14[k]
                   - f_17 * id1_14[k]
                   + pa_y[k] * kd_42[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, pb_y, id0_28, id1_28, kp_30, \
                         kp_39, kp_40, kd_56, lp_51, lp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_4 * kp_39[k]
                   + pb_x[k] * lp_51[k];

        t_140[k] = f_4 * kp_40[k]
                   + pb_x[k] * lp_52[k];

        t_141[k] = f_5 * id0_28[k]
                   - f_6 * id1_28[k]
                   + pa_x[k] * kd_56[k];

        t_142[k] = f_18 * kp_30[k]
                   + pb_y[k] * lp_52[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pa_y, pb_x, id0_17, id0_29, id1_17, \
                         id1_29, kp_41, kd_45, kd_57, lp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_5 * id0_29[k]
                   - f_6 * id1_29[k]
                   + pa_x[k] * kd_57[k];

        t_144[k] = f_11 * id0_17[k]
                   - f_12 * id1_17[k]
                   + pa_y[k] * kd_45[k];

        t_145[k] = f_4 * kp_41[k]
                   + pb_x[k] * lp_53[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_x, pb_y, id0_30, id0_31, id1_30, \
                         id1_31, kp_32, kp_42, kd_58, kd_59, lp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_4 * kp_42[k]
                   + pb_x[k] * lp_54[k];

        t_147[k] = f_5 * id0_30[k]
                   - f_6 * id1_30[k]
                   + pa_x[k] * kd_58[k];

        t_148[k] = f_19 * kp_32[k]
                   + pb_y[k] * lp_54[k];

        t_149[k] = f_5 * id0_31[k]
                   - f_6 * id1_31[k]
                   + pa_x[k] * kd_59[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_y, pb_x, id0_18, id1_18, kp_43, kp_44, kd_48, \
                         lp_55, lp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * id0_18[k]
                   - f_6 * id1_18[k]
                   + pa_y[k] * kd_48[k];

        t_151[k] = f_4 * kp_43[k]
                   + pb_x[k] * lp_55[k];

        t_152[k] = f_4 * kp_44[k]
                   + pb_x[k] * lp_56[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_x, pa_y, pb_y, id0_32, id0_33, id1_32, \
                         id1_33, kp_34, kd_49, kd_60, kd_61, lp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_5 * id0_32[k]
                   - f_6 * id1_32[k]
                   + pa_x[k] * kd_60[k];

        t_154[k] = f_4 * kp_34[k]
                   + pb_y[k] * lp_56[k];

        t_155[k] = f_5 * id0_33[k]
                   - f_6 * id1_33[k]
                   + pa_x[k] * kd_61[k];

        t_156[k] = pa_y[k] * kd_49[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_y, pb_x, pb_y, kp_35, kp_36, \
                         kp_45, kd_50, kd_51, kd_52, lp_57, lp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_4 * kp_45[k]
                   + pb_x[k] * lp_57[k];

        t_158[k] = pa_y[k] * kd_50[k];

        t_159[k] = f_4 * kp_35[k]
                   + pa_y[k] * kd_51[k];

        t_160[k] = f_10 * kp_36[k]
                   + pb_y[k] * lp_58[k];

        t_161[k] = pa_y[k] * kd_52[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_z, pb_x, pb_y, id0_18, id1_18, kp_46, \
                         kd_49, ls0_10, ls1_10, lp_59, lp_60, lp_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_8 * id0_18[k]
                   - f_9 * id1_18[k]
                   + pa_z[k] * kd_49[k];

        t_163[k] = pb_y[k] * lp_59[k];

        t_164[k] = f_4 * kp_46[k]
                   + pb_x[k] * lp_61[k];

        t_165[k] = f_1 * ls0_10[k]
                   - f_2 * ls1_10[k]
                   + pb_y[k] * lp_60[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_x, pb_x, pb_y, id0_35, id1_35, kp_47, \
                         kp_48, kd_64, kd_65, lp_61, lp_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_y[k] * lp_61[k];

        t_167[k] = f_5 * id0_35[k]
                   - f_6 * id1_35[k]
                   + pa_x[k] * kd_64[k];

        t_168[k] = f_4 * kp_47[k]
                   + pa_x[k] * kd_65[k];

        t_169[k] = f_10 * kp_48[k]
                   + pb_x[k] * lp_63[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, pa_x, pa_z, pb_z, kd_53, \
                         kd_54, kd_66, kd_67, lp_62, lp_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_z[k] * lp_62[k];

        t_171[k] = pa_x[k] * kd_66[k];

        t_172[k] = pb_z[k] * lp_63[k];

        t_173[k] = pa_x[k] * kd_67[k];

        t_174[k] = pa_z[k] * kd_53[k];

        t_175[k] = pa_z[k] * kd_54[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, pa_x, pb_x, kp_50, kp_51, kd_68, \
                         kd_69, kd_70, kd_71, lp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * kp_50[k]
                   + pb_x[k] * lp_64[k];

        t_177[k] = pa_x[k] * kd_68[k];

        t_178[k] = pa_x[k] * kd_69[k];

        t_179[k] = pa_x[k] * kd_70[k];

        t_180[k] = f_4 * kp_51[k]
                   + pa_x[k] * kd_71[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, pa_x, pb_x, kp_52, kp_53, kd_72, \
                         kd_73, kd_74, lp_65, lp_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * kp_52[k]
                   + pb_x[k] * lp_65[k];

        t_182[k] = f_10 * kp_53[k]
                   + pb_x[k] * lp_66[k];

        t_183[k] = pa_x[k] * kd_72[k];

        t_184[k] = pa_x[k] * kd_73[k];

        t_185[k] = pa_x[k] * kd_74[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_x, pb_x, kp_54, kp_55, kp_56, \
                         kd_75, kd_76, kd_77, lp_67, lp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_4 * kp_54[k]
                   + pa_x[k] * kd_75[k];

        t_187[k] = f_10 * kp_55[k]
                   + pb_x[k] * lp_67[k];

        t_188[k] = f_10 * kp_56[k]
                   + pb_x[k] * lp_68[k];

        t_189[k] = pa_x[k] * kd_76[k];

        t_190[k] = pa_x[k] * kd_77[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pa_x, pb_x, kp_57, kp_58, kp_59, \
                         kd_78, kd_79, kd_80, lp_69, lp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = pa_x[k] * kd_78[k];

        t_192[k] = f_4 * kp_57[k]
                   + pa_x[k] * kd_79[k];

        t_193[k] = f_10 * kp_58[k]
                   + pb_x[k] * lp_69[k];

        t_194[k] = f_10 * kp_59[k]
                   + pb_x[k] * lp_70[k];

        t_195[k] = pa_x[k] * kd_80[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_x, pb_x, kp_60, kp_61, kp_62, \
                         kd_81, kd_82, kd_83, lp_71, lp_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * kd_81[k];

        t_197[k] = pa_x[k] * kd_82[k];

        t_198[k] = f_4 * kp_60[k]
                   + pa_x[k] * kd_83[k];

        t_199[k] = f_10 * kp_61[k]
                   + pb_x[k] * lp_71[k];

        t_200[k] = f_10 * kp_62[k]
                   + pb_x[k] * lp_72[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pa_x, pa_y, pb_x, kp_63, \
                         kd_62, kd_63, kd_84, kd_85, kd_86, lp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_x[k] * kd_84[k];

        t_202[k] = pa_x[k] * kd_85[k];

        t_203[k] = pa_x[k] * kd_86[k];

        t_204[k] = pa_y[k] * kd_62[k];

        t_205[k] = f_10 * kp_63[k]
                   + pb_x[k] * lp_73[k];

        t_206[k] = pa_y[k] * kd_63[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_x, pb_y, kp_65, kd_87, kd_88, \
                         kd_89, kd_90, lp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_x[k] * kd_87[k];

        t_208[k] = pa_x[k] * kd_88[k];

        t_209[k] = pa_x[k] * kd_89[k];

        t_210[k] = f_4 * kp_65[k]
                   + pa_x[k] * kd_90[k];

        t_211[k] = pb_y[k] * lp_74[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pb_x, pb_y, kp_67, kd_91, \
                         kd_92, ls0_11, ls1_11, lp_75, lp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_10 * kp_67[k]
                   + pb_x[k] * lp_75[k];

        t_213[k] = pa_x[k] * kd_91[k];

        t_214[k] = pb_y[k] * lp_75[k];

        t_215[k] = pa_x[k] * kd_92[k];

        t_216[k] = f_1 * ls0_11[k]
                   - f_2 * ls1_11[k]
                   + pb_x[k] * lp_76[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, t_222, pa_z, pb_x, pb_y, pb_z, \
                         kp_48, kd_65, ls0_11, ls1_11, lp_77, lp_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_x[k] * lp_77[k];

        t_218[k] = pb_x[k] * lp_78[k];

        t_219[k] = f_0 * kp_48[k]
                   + f_1 * ls0_11[k]
                   - f_2 * ls1_11[k]
                   + pb_y[k] * lp_77[k];

        t_220[k] = pb_z[k] * lp_77[k];

        t_221[k] = f_1 * ls0_11[k]
                   - f_2 * ls1_11[k]
                   + pb_z[k] * lp_78[k];

        t_222[k] = pa_z[k] * kd_65[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_y, kp_49, kp_50, \
                         kd_66, kd_67, lp_79, lp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_x[k] * lp_79[k];

        t_224[k] = pb_x[k] * lp_80[k];

        t_225[k] = pa_z[k] * kd_66[k];

        t_226[k] = f_3 * kp_50[k]
                   + pb_y[k] * lp_80[k];

        t_227[k] = f_4 * kp_49[k]
                   + pa_z[k] * kd_67[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_z, pb_x, id0_26, id1_26, kd_68, \
                         ls0_12, ls1_12, lp_81, lp_82, lp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_1 * ls0_12[k]
                   - f_2 * ls1_12[k]
                   + pb_x[k] * lp_81[k];

        t_229[k] = pb_x[k] * lp_82[k];

        t_230[k] = pb_x[k] * lp_83[k];

        t_231[k] = f_5 * id0_26[k]
                   - f_6 * id1_26[k]
                   + pa_z[k] * kd_68[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_y, pb_x, pb_y, id0_29, id1_29, kp_53, \
                         kd_74, ls0_13, ls1_13, lp_83, lp_84, lp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_7 * kp_53[k]
                   + pb_y[k] * lp_83[k];

        t_233[k] = f_8 * id0_29[k]
                   - f_9 * id1_29[k]
                   + pa_y[k] * kd_74[k];

        t_234[k] = f_1 * ls0_13[k]
                   - f_2 * ls1_13[k]
                   + pb_x[k] * lp_84[k];

        t_235[k] = pb_x[k] * lp_85[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pa_z, pb_x, pb_y, id0_27, id0_31, \
                         id1_27, id1_31, kp_56, kd_72, kd_78, lp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_x[k] * lp_86[k];

        t_237[k] = f_11 * id0_27[k]
                   - f_12 * id1_27[k]
                   + pa_z[k] * kd_72[k];

        t_238[k] = f_13 * kp_56[k]
                   + pb_y[k] * lp_86[k];

        t_239[k] = f_14 * id0_31[k]
                   - f_15 * id1_31[k]
                   + pa_y[k] * kd_78[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pa_z, pb_x, id0_28, id1_28, kd_76, \
                         ls0_14, ls1_14, lp_87, lp_88, lp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_1 * ls0_14[k]
                   - f_2 * ls1_14[k]
                   + pb_x[k] * lp_87[k];

        t_241[k] = pb_x[k] * lp_88[k];

        t_242[k] = pb_x[k] * lp_89[k];

        t_243[k] = f_16 * id0_28[k]
                   - f_17 * id1_28[k]
                   + pa_z[k] * kd_76[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_y, pb_x, pb_y, id0_33, id1_33, kp_59, \
                         kd_82, ls0_15, ls1_15, lp_89, lp_90, lp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_18 * kp_59[k]
                   + pb_y[k] * lp_89[k];

        t_245[k] = f_16 * id0_33[k]
                   - f_17 * id1_33[k]
                   + pa_y[k] * kd_82[k];

        t_246[k] = f_1 * ls0_15[k]
                   - f_2 * ls1_15[k]
                   + pb_x[k] * lp_90[k];

        t_247[k] = pb_x[k] * lp_91[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_y, pa_z, pb_x, pb_y, id0_30, id0_34, \
                         id1_30, id1_34, kp_62, kd_80, kd_86, lp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pb_x[k] * lp_92[k];

        t_249[k] = f_14 * id0_30[k]
                   - f_15 * id1_30[k]
                   + pa_z[k] * kd_80[k];

        t_250[k] = f_19 * kp_62[k]
                   + pb_y[k] * lp_92[k];

        t_251[k] = f_11 * id0_34[k]
                   - f_12 * id1_34[k]
                   + pa_y[k] * kd_86[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_z, pb_x, id0_32, id1_32, kd_84, \
                         ls0_16, ls1_16, lp_93, lp_94, lp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * ls0_16[k]
                   - f_2 * ls1_16[k]
                   + pb_x[k] * lp_93[k];

        t_253[k] = pb_x[k] * lp_94[k];

        t_254[k] = pb_x[k] * lp_95[k];

        t_255[k] = f_8 * id0_32[k]
                   - f_9 * id1_32[k]
                   + pa_z[k] * kd_84[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pa_y, pb_x, pb_y, id0_35, id1_35, \
                         kp_64, kd_89, kd_90, lp_95, lp_96, lp_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_4 * kp_64[k]
                   + pb_y[k] * lp_95[k];

        t_257[k] = f_5 * id0_35[k]
                   - f_6 * id1_35[k]
                   + pa_y[k] * kd_89[k];

        t_258[k] = pa_y[k] * kd_90[k];

        t_259[k] = pb_x[k] * lp_96[k];

        t_260[k] = pb_x[k] * lp_97[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pb_x, pb_y, kp_66, kp_67, kd_91, \
                         kd_92, ls0_17, ls1_17, lp_97, lp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_4 * kp_66[k]
                   + pa_y[k] * kd_91[k];

        t_262[k] = f_10 * kp_67[k]
                   + pb_y[k] * lp_97[k];

        t_263[k] = pa_y[k] * kd_92[k];

        t_264[k] = f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_x[k] * lp_98[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pb_x, pb_y, pb_z, kp_67, ls0_17, \
                         ls1_17, lp_99, lp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = pb_x[k] * lp_99[k];

        t_266[k] = pb_x[k] * lp_100[k];

        t_267[k] = f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_y[k] * lp_99[k];

        t_268[k] = pb_y[k] * lp_100[k];

        t_269[k] = f_0 * kp_67[k]
                   + f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_z[k] * lp_100[k];
    }
}

auto
compute_prim_ld_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 3.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.5 / alpha;
    const auto f_9 = 2.5 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.5 / p;
    const auto f_14 = 2.0 / alpha;
    const auto f_15 = 2.0 * beta / (alpha * p);
    const auto f_16 = 1.5 / alpha;
    const auto f_17 = 1.5 * beta / (alpha * p);
    const auto f_18 = 2.0 / p;
    const auto f_19 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_1 = buffer.data(id0 + 1);
    const auto *id0_2 = buffer.data(id0 + 2);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_7 = buffer.data(id0 + 7);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_13 = buffer.data(id0 + 13);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_15 = buffer.data(id0 + 15);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_19 = buffer.data(id0 + 19);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_25 = buffer.data(id0 + 25);
    const auto *id0_26 = buffer.data(id0 + 26);
    const auto *id0_27 = buffer.data(id0 + 27);
    const auto *id0_28 = buffer.data(id0 + 28);
    const auto *id0_29 = buffer.data(id0 + 29);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_31 = buffer.data(id0 + 31);
    const auto *id0_32 = buffer.data(id0 + 32);
    const auto *id0_33 = buffer.data(id0 + 33);
    const auto *id0_34 = buffer.data(id0 + 34);
    const auto *id0_35 = buffer.data(id0 + 35);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_7 = buffer.data(id1 + 7);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_13 = buffer.data(id1 + 13);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_19 = buffer.data(id1 + 19);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_25 = buffer.data(id1 + 25);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_27 = buffer.data(id1 + 27);
    const auto *id1_29 = buffer.data(id1 + 29);
    const auto *id1_31 = buffer.data(id1 + 31);
    const auto *id1_32 = buffer.data(id1 + 32);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_34 = buffer.data(id1 + 34);
    const auto *id1_35 = buffer.data(id1 + 35);
    const auto *id1_37 = buffer.data(id1 + 37);
    const auto *id1_39 = buffer.data(id1 + 39);
    const auto *id1_41 = buffer.data(id1 + 41);
    const auto *id1_45 = buffer.data(id1 + 45);
    const auto *id1_47 = buffer.data(id1 + 47);
    const auto *id1_49 = buffer.data(id1 + 49);
    const auto *id1_51 = buffer.data(id1 + 51);
    const auto *id1_53 = buffer.data(id1 + 53);
    const auto *id1_55 = buffer.data(id1 + 55);
    const auto *id1_58 = buffer.data(id1 + 58);
    const auto *id1_61 = buffer.data(id1 + 61);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_24 = buffer.data(lp + 24);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_33 = buffer.data(lp + 33);
    const auto *lp_34 = buffer.data(lp + 34);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_36 = buffer.data(lp + 36);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_39 = buffer.data(lp + 39);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_48 = buffer.data(lp + 48);
    const auto *lp_49 = buffer.data(lp + 49);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_51 = buffer.data(lp + 51);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_54 = buffer.data(lp + 54);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_57 = buffer.data(lp + 57);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, kp_0, kd_0, ls0_0, \
                         ls1_0, lp_0, lp_1, lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = pb_z[k] * lp_0[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_3[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];

        t_4[k] = pa_y[k] * kd_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_x, kp_1, kp_3, kp_4, kd_0, \
                         kd_1, kd_2, lp_3, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * kp_3[k]
                 + pb_x[k] * lp_3[k];

        t_6[k] = f_4 * kp_1[k]
                 + pa_y[k] * kd_1[k];

        t_7[k] = pa_y[k] * kd_2[k];

        t_8[k] = pa_z[k] * kd_0[k];

        t_9[k] = f_3 * kp_4[k]
                 + pb_x[k] * lp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, id0_0, id1_0, kp_2, kp_5, \
                         kd_1, kd_2, kd_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * kd_1[k];

        t_11[k] = f_4 * kp_2[k]
                  + pa_z[k] * kd_2[k];

        t_12[k] = f_5 * id0_0[k]
                  - f_6 * id1_0[k]
                  + pa_y[k] * kd_3[k];

        t_13[k] = f_7 * kp_5[k]
                  + pb_x[k] * lp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, id0_4, id1_8, kd_4, kd_8, ls0_1, \
                         ls1_1, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * id0_4[k]
                  - f_9 * id1_8[k]
                  + pa_x[k] * kd_8[k];

        t_15[k] = f_1 * ls0_1[k]
                  - f_2 * ls1_1[k]
                  + pb_z[k] * lp_6[k];

        t_16[k] = pa_z[k] * kd_4[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_x, pb_y, id0_0, id1_0, kp_4, \
                         kp_9, kd_5, kd_6, lp_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * kp_4[k]
                  + pb_y[k] * lp_7[k];

        t_18[k] = pa_y[k] * kd_6[k];

        t_19[k] = f_5 * id0_0[k]
                  - f_6 * id1_0[k]
                  + pa_z[k] * kd_5[k];

        t_20[k] = f_7 * kp_9[k]
                  + pb_x[k] * lp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, id0_1, id0_6, id1_3, id1_12, \
                         kd_7, kd_12, ls0_2, ls1_2, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * ls0_2[k]
                  - f_2 * ls1_2[k]
                  + pb_y[k] * lp_8[k];

        t_22[k] = f_8 * id0_6[k]
                  - f_9 * id1_12[k]
                  + pa_x[k] * kd_12[k];

        t_23[k] = f_11 * id0_1[k]
                  - f_12 * id1_3[k]
                  + pa_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_z, id0_8, id1_14, kp_10, kd_14, \
                         ls0_3, ls1_3, lp_10, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_13 * kp_10[k]
                  + pb_x[k] * lp_10[k];

        t_25[k] = f_14 * id0_8[k]
                  - f_15 * id1_14[k]
                  + pa_x[k] * kd_14[k];

        t_26[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, kp_6, kp_7, kd_7, \
                         kd_8, kd_9, kd_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * kd_7[k];

        t_28[k] = pa_z[k] * kd_8[k];

        t_29[k] = f_4 * kp_7[k]
                  + pb_y[k] * lp_12[k];

        t_30[k] = f_4 * kp_6[k]
                  + pa_z[k] * kd_9[k];

        t_31[k] = pa_y[k] * kd_10[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pa_z, pb_y, id0_2, id1_5, kp_8, kp_9, \
                         kd_10, kd_11, kd_12, lp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * kp_8[k]
                  + pa_y[k] * kd_11[k];

        t_33[k] = f_10 * kp_9[k]
                  + pb_y[k] * lp_13[k];

        t_34[k] = pa_y[k] * kd_12[k];

        t_35[k] = f_11 * id0_2[k]
                  - f_12 * id1_5[k]
                  + pa_z[k] * kd_10[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_x, pb_y, id0_11, id1_19, kp_15, kd_19, \
                         ls0_4, ls1_4, lp_14, lp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * kp_15[k]
                  + pb_x[k] * lp_15[k];

        t_37[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_14[k];

        t_38[k] = f_14 * id0_11[k]
                  - f_15 * id1_19[k]
                  + pa_x[k] * kd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_y, pb_x, id0_3, id0_13, id1_7, id1_21, \
                         kp_16, kd_13, kd_21, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_16 * id0_3[k]
                  - f_17 * id1_7[k]
                  + pa_y[k] * kd_13[k];

        t_40[k] = f_18 * kp_16[k]
                  + pb_x[k] * lp_16[k];

        t_41[k] = f_16 * id0_13[k]
                  - f_17 * id1_21[k]
                  + pa_x[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, kp_12, kd_13, kd_14, ls0_5, \
                         ls1_5, lp_17, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_17[k];

        t_43[k] = pa_z[k] * kd_13[k];

        t_44[k] = pa_z[k] * kd_14[k];

        t_45[k] = f_19 * kp_12[k]
                  + pb_y[k] * lp_18[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pa_y, pa_z, id0_5, id0_15, id1_10, id1_24, \
                         kp_11, kd_15, kd_16, kd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * kp_11[k]
                  + pa_z[k] * kd_15[k];

        t_47[k] = f_5 * id0_5[k]
                  - f_6 * id1_10[k]
                  + pa_y[k] * kd_16[k];

        t_48[k] = f_16 * id0_15[k]
                  - f_17 * id1_24[k]
                  + pa_x[k] * kd_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_y, pb_y, id0_16, id1_25, kp_13, \
                         kp_14, kd_17, kd_18, kd_25, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_4 * kp_13[k]
                  + pb_y[k] * lp_19[k];

        t_50[k] = f_16 * id0_16[k]
                  - f_17 * id1_25[k]
                  + pa_x[k] * kd_25[k];

        t_51[k] = pa_y[k] * kd_17[k];

        t_52[k] = f_4 * kp_14[k]
                  + pa_y[k] * kd_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_x, pb_y, id0_5, id1_10, kp_15, \
                         kp_22, kd_17, kd_19, lp_20, lp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * kp_15[k]
                  + pb_y[k] * lp_20[k];

        t_54[k] = pa_y[k] * kd_19[k];

        t_55[k] = f_16 * id0_5[k]
                  - f_17 * id1_10[k]
                  + pa_z[k] * kd_17[k];

        t_56[k] = f_18 * kp_22[k]
                  + pb_x[k] * lp_22[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pa_y, pb_y, id0_7, id0_19, id1_13, id1_29, \
                         kd_20, kd_29, ls0_6, ls1_6, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_21[k];

        t_58[k] = f_16 * id0_19[k]
                  - f_17 * id1_29[k]
                  + pa_x[k] * kd_29[k];

        t_59[k] = f_14 * id0_7[k]
                  - f_15 * id1_13[k]
                  + pa_y[k] * kd_20[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pb_x, pb_z, id0_20, id1_31, kp_23, kd_31, \
                         ls0_7, ls1_7, lp_23, lp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_19 * kp_23[k]
                  + pb_x[k] * lp_23[k];

        t_61[k] = f_11 * id0_20[k]
                  - f_12 * id1_31[k]
                  + pa_x[k] * kd_31[k];

        t_62[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_24[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_y, kp_17, kp_18, kd_20, kd_21, \
                         kd_22, lp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * kd_20[k];

        t_64[k] = pa_z[k] * kd_21[k];

        t_65[k] = f_18 * kp_18[k]
                  + pb_y[k] * lp_25[k];

        t_66[k] = f_4 * kp_17[k]
                  + pa_z[k] * kd_22[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_x, pa_y, pb_y, id0_9, id0_21, id1_16, id1_32, \
                         kp_19, kd_23, kd_34, lp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * id0_9[k]
                  - f_12 * id1_16[k]
                  + pa_y[k] * kd_23[k];

        t_68[k] = f_11 * id0_21[k]
                  - f_12 * id1_32[k]
                  + pa_x[k] * kd_34[k];

        t_69[k] = f_19 * kp_19[k]
                  + pb_y[k] * lp_26[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_x, pa_y, id0_10, id0_22, id0_23, id1_17, id1_33, \
                         id1_34, kd_26, kd_35, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_11 * id0_22[k]
                  - f_12 * id1_33[k]
                  + pa_x[k] * kd_35[k];

        t_71[k] = f_5 * id0_10[k]
                  - f_6 * id1_17[k]
                  + pa_y[k] * kd_26[k];

        t_72[k] = f_11 * id0_23[k]
                  - f_12 * id1_34[k]
                  + pa_x[k] * kd_37[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pa_y, pb_y, id0_24, id1_35, kp_20, \
                         kp_21, kd_27, kd_28, kd_38, lp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * kp_20[k]
                  + pb_y[k] * lp_27[k];

        t_74[k] = f_11 * id0_24[k]
                  - f_12 * id1_35[k]
                  + pa_x[k] * kd_38[k];

        t_75[k] = pa_y[k] * kd_27[k];

        t_76[k] = f_4 * kp_21[k]
                  + pa_y[k] * kd_28[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pa_z, pb_x, pb_y, id0_10, id1_17, \
                         kp_22, kp_30, kd_27, kd_29, lp_28, lp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_10 * kp_22[k]
                  + pb_y[k] * lp_28[k];

        t_78[k] = pa_y[k] * kd_29[k];

        t_79[k] = f_14 * id0_10[k]
                  - f_15 * id1_17[k]
                  + pa_z[k] * kd_27[k];

        t_80[k] = f_19 * kp_30[k]
                  + pb_x[k] * lp_30[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, pa_y, pb_y, id0_12, id0_25, id1_20, id1_37, \
                         kd_30, kd_42, ls0_8, ls1_8, lp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_29[k];

        t_82[k] = f_11 * id0_25[k]
                  - f_12 * id1_37[k]
                  + pa_x[k] * kd_42[k];

        t_83[k] = f_8 * id0_12[k]
                  - f_9 * id1_20[k]
                  + pa_y[k] * kd_30[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_x, pb_x, pb_z, id0_26, id1_39, kp_31, kd_44, \
                         ls0_9, ls1_9, lp_31, lp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * kp_31[k]
                  + pb_x[k] * lp_31[k];

        t_85[k] = f_5 * id0_26[k]
                  - f_6 * id1_39[k]
                  + pa_x[k] * kd_44[k];

        t_86[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_32[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_z, pb_y, kp_24, kp_25, kd_30, kd_31, \
                         kd_32, lp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * kd_30[k];

        t_88[k] = pa_z[k] * kd_31[k];

        t_89[k] = f_13 * kp_25[k]
                  + pb_y[k] * lp_33[k];

        t_90[k] = f_4 * kp_24[k]
                  + pa_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pa_x, pa_y, pb_y, id0_14, id0_28, id1_23, id1_45, \
                         kp_26, kd_33, kd_45, lp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_16 * id0_14[k]
                  - f_17 * id1_23[k]
                  + pa_y[k] * kd_33[k];

        t_92[k] = f_5 * id0_28[k]
                  - f_6 * id1_45[k]
                  + pa_x[k] * kd_45[k];

        t_93[k] = f_18 * kp_26[k]
                  + pb_y[k] * lp_34[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_x, pa_y, id0_17, id0_29, id0_30, id1_26, id1_47, \
                         id1_49, kd_36, kd_46, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * id0_29[k]
                  - f_6 * id1_47[k]
                  + pa_x[k] * kd_46[k];

        t_95[k] = f_11 * id0_17[k]
                  - f_12 * id1_26[k]
                  + pa_y[k] * kd_36[k];

        t_96[k] = f_5 * id0_30[k]
                  - f_6 * id1_49[k]
                  + pa_x[k] * kd_47[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pa_x, pa_y, pb_y, id0_18, id0_31, id1_27, id1_51, \
                         kp_27, kd_39, kd_48, lp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_19 * kp_27[k]
                  + pb_y[k] * lp_35[k];

        t_98[k] = f_5 * id0_31[k]
                  - f_6 * id1_51[k]
                  + pa_x[k] * kd_48[k];

        t_99[k] = f_5 * id0_18[k]
                  - f_6 * id1_27[k]
                  + pa_y[k] * kd_39[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pa_y, pb_y, id0_32, id0_33, id1_53, \
                         id1_55, kp_28, kd_40, kd_49, kd_50, lp_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * id0_32[k]
                   - f_6 * id1_53[k]
                   + pa_x[k] * kd_49[k];

        t_101[k] = f_4 * kp_28[k]
                   + pb_y[k] * lp_36[k];

        t_102[k] = f_5 * id0_33[k]
                   - f_6 * id1_55[k]
                   + pa_x[k] * kd_50[k];

        t_103[k] = pa_y[k] * kd_40[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pa_z, pb_y, id0_18, id1_27, kp_29, \
                         kp_30, kd_40, kd_41, kd_42, lp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_4 * kp_29[k]
                   + pa_y[k] * kd_41[k];

        t_105[k] = f_10 * kp_30[k]
                   + pb_y[k] * lp_37[k];

        t_106[k] = pa_y[k] * kd_42[k];

        t_107[k] = f_8 * id0_18[k]
                   - f_9 * id1_27[k]
                   + pa_z[k] * kd_40[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_x, pb_x, pb_y, id0_35, id1_61, kp_32, kd_52, \
                         ls0_10, ls1_10, lp_38, lp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_4 * kp_32[k]
                   + pb_x[k] * lp_39[k];

        t_109[k] = f_1 * ls0_10[k]
                   - f_2 * ls1_10[k]
                   + pb_y[k] * lp_38[k];

        t_110[k] = f_5 * id0_35[k]
                   - f_6 * id1_61[k]
                   + pa_x[k] * kd_52[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pa_x, pa_z, pb_x, kp_33, kp_34, \
                         kd_43, kd_53, kd_54, kd_55, lp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_4 * kp_33[k]
                   + pa_x[k] * kd_53[k];

        t_112[k] = f_10 * kp_34[k]
                   + pb_x[k] * lp_40[k];

        t_113[k] = pa_x[k] * kd_54[k];

        t_114[k] = pa_x[k] * kd_55[k];

        t_115[k] = pa_z[k] * kd_43[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, pa_x, kp_37, kd_57, kd_58, \
                         kd_59, kd_60, kd_61, kd_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_x[k] * kd_57[k];

        t_117[k] = pa_x[k] * kd_58[k];

        t_118[k] = f_4 * kp_37[k]
                   + pa_x[k] * kd_59[k];

        t_119[k] = pa_x[k] * kd_60[k];

        t_120[k] = pa_x[k] * kd_61[k];

        t_121[k] = pa_x[k] * kd_62[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, pa_x, kp_39, kp_41, kd_63, \
                         kd_64, kd_65, kd_66, kd_67, kd_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * kp_39[k]
                   + pa_x[k] * kd_63[k];

        t_123[k] = pa_x[k] * kd_64[k];

        t_124[k] = pa_x[k] * kd_65[k];

        t_125[k] = pa_x[k] * kd_66[k];

        t_126[k] = f_4 * kp_41[k]
                   + pa_x[k] * kd_67[k];

        t_127[k] = pa_x[k] * kd_68[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pa_x, kp_43, kd_69, kd_70, \
                         kd_71, kd_72, kd_73, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pa_x[k] * kd_69[k];

        t_129[k] = pa_x[k] * kd_70[k];

        t_130[k] = f_4 * kp_43[k]
                   + pa_x[k] * kd_71[k];

        t_131[k] = pa_x[k] * kd_72[k];

        t_132[k] = pa_x[k] * kd_73[k];

        t_133[k] = pa_x[k] * kd_74[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pa_y, pb_x, kp_46, kp_48, \
                         kd_51, kd_75, kd_76, kd_78, lp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_y[k] * kd_51[k];

        t_135[k] = pa_x[k] * kd_75[k];

        t_136[k] = pa_x[k] * kd_76[k];

        t_137[k] = f_4 * kp_46[k]
                   + pa_x[k] * kd_78[k];

        t_138[k] = f_10 * kp_48[k]
                   + pb_x[k] * lp_41[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_x, pb_y, pb_z, kp_34, \
                         kd_79, kd_80, ls0_11, ls1_11, lp_42, lp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_x[k] * kd_79[k];

        t_140[k] = pa_x[k] * kd_80[k];

        t_141[k] = f_1 * ls0_11[k]
                   - f_2 * ls1_11[k]
                   + pb_x[k] * lp_42[k];

        t_142[k] = f_0 * kp_34[k]
                   + f_1 * ls0_11[k]
                   - f_2 * ls1_11[k]
                   + pb_y[k] * lp_43[k];

        t_143[k] = pb_z[k] * lp_43[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_z, pb_y, pb_z, kp_36, kd_53, kd_54, \
                         ls0_11, ls1_11, lp_44, lp_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_1 * ls0_11[k]
                   - f_2 * ls1_11[k]
                   + pb_z[k] * lp_44[k];

        t_145[k] = pa_z[k] * kd_53[k];

        t_146[k] = pa_z[k] * kd_54[k];

        t_147[k] = f_3 * kp_36[k]
                   + pb_y[k] * lp_45[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_z, pb_x, id0_26, id1_39, kp_35, kd_55, kd_56, \
                         ls0_12, ls1_12, lp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_4 * kp_35[k]
                   + pa_z[k] * kd_55[k];

        t_149[k] = f_1 * ls0_12[k]
                   - f_2 * ls1_12[k]
                   + pb_x[k] * lp_46[k];

        t_150[k] = f_5 * id0_26[k]
                   - f_6 * id1_39[k]
                   + pa_z[k] * kd_56[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_y, pb_x, pb_y, id0_29, id1_47, kp_38, kd_62, \
                         ls0_13, ls1_13, lp_47, lp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_7 * kp_38[k]
                   + pb_y[k] * lp_47[k];

        t_152[k] = f_8 * id0_29[k]
                   - f_9 * id1_47[k]
                   + pa_y[k] * kd_62[k];

        t_153[k] = f_1 * ls0_13[k]
                   - f_2 * ls1_13[k]
                   + pb_x[k] * lp_48[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_y, pa_z, pb_y, id0_27, id0_31, id1_41, \
                         id1_51, kp_40, kd_60, kd_66, lp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_11 * id0_27[k]
                   - f_12 * id1_41[k]
                   + pa_z[k] * kd_60[k];

        t_155[k] = f_13 * kp_40[k]
                   + pb_y[k] * lp_49[k];

        t_156[k] = f_14 * id0_31[k]
                   - f_15 * id1_51[k]
                   + pa_y[k] * kd_66[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pa_z, pb_x, pb_y, id0_28, id1_45, kp_42, kd_64, \
                         ls0_14, ls1_14, lp_50, lp_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_1 * ls0_14[k]
                   - f_2 * ls1_14[k]
                   + pb_x[k] * lp_50[k];

        t_158[k] = f_16 * id0_28[k]
                   - f_17 * id1_45[k]
                   + pa_z[k] * kd_64[k];

        t_159[k] = f_18 * kp_42[k]
                   + pb_y[k] * lp_51[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_y, pa_z, pb_x, id0_30, id0_33, id1_49, \
                         id1_55, kd_68, kd_70, ls0_15, ls1_15, lp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * id0_33[k]
                   - f_17 * id1_55[k]
                   + pa_y[k] * kd_70[k];

        t_161[k] = f_1 * ls0_15[k]
                   - f_2 * ls1_15[k]
                   + pb_x[k] * lp_52[k];

        t_162[k] = f_14 * id0_30[k]
                   - f_15 * id1_49[k]
                   + pa_z[k] * kd_68[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_y, pb_x, pb_y, id0_34, id1_58, kp_44, kd_74, \
                         ls0_16, ls1_16, lp_53, lp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_19 * kp_44[k]
                   + pb_y[k] * lp_53[k];

        t_164[k] = f_11 * id0_34[k]
                   - f_12 * id1_58[k]
                   + pa_y[k] * kd_74[k];

        t_165[k] = f_1 * ls0_16[k]
                   - f_2 * ls1_16[k]
                   + pb_x[k] * lp_54[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_y, pa_z, pb_y, id0_32, id0_35, id1_53, \
                         id1_61, kp_45, kd_72, kd_77, kd_78, lp_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_8 * id0_32[k]
                   - f_9 * id1_53[k]
                   + pa_z[k] * kd_72[k];

        t_167[k] = f_4 * kp_45[k]
                   + pb_y[k] * lp_55[k];

        t_168[k] = f_5 * id0_35[k]
                   - f_6 * id1_61[k]
                   + pa_y[k] * kd_77[k];

        t_169[k] = pa_y[k] * kd_78[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_y, pb_x, pb_y, kp_47, kp_48, kd_79, \
                         kd_80, ls0_17, ls1_17, lp_56, lp_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_4 * kp_47[k]
                   + pa_y[k] * kd_79[k];

        t_171[k] = f_10 * kp_48[k]
                   + pb_y[k] * lp_56[k];

        t_172[k] = pa_y[k] * kd_80[k];

        t_173[k] = f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_x[k] * lp_57[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_y, pb_z, kp_48, ls0_17, ls1_17, lp_58, \
                         lp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_y[k] * lp_58[k];

        t_175[k] = pb_y[k] * lp_59[k];

        t_176[k] = f_0 * kp_48[k]
                   + f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_z[k] * lp_59[k];
    }
}

auto
compute_prim_ld_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / alpha;
    const auto f_6 = 2.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 2.0 / alpha;
    const auto f_10 = 2.0 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_1 = buffer.data(id0 + 1);
    const auto *id0_2 = buffer.data(id0 + 2);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_7 = buffer.data(id0 + 7);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_13 = buffer.data(id0 + 13);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_15 = buffer.data(id0 + 15);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_19 = buffer.data(id0 + 19);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_25 = buffer.data(id0 + 25);
    const auto *id0_26 = buffer.data(id0 + 26);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_1 = buffer.data(id1 + 1);
    const auto *id1_2 = buffer.data(id1 + 2);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_4 = buffer.data(id1 + 4);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_7 = buffer.data(id1 + 7);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_9 = buffer.data(id1 + 9);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_11 = buffer.data(id1 + 11);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_13 = buffer.data(id1 + 13);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_15 = buffer.data(id1 + 15);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_19 = buffer.data(id1 + 19);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_22 = buffer.data(id1 + 22);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_25 = buffer.data(id1 + 25);
    const auto *id1_26 = buffer.data(id1 + 26);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, lp_0, lp_1, \
                         lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, id0_0, id0_4, id1_0, id1_4, kd_3, \
                         kd_6, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_y[k] * kd_3[k];

        t_4[k] = f_5 * id0_4[k]
                 - f_6 * id1_4[k]
                 + pa_x[k] * kd_6[k];

        t_5[k] = f_1 * ls0_1[k]
                 - f_2 * ls1_1[k]
                 + pb_z[k] * lp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, id0_0, id0_6, id1_0, id1_6, kd_4, \
                         kd_10, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_z[k] * kd_4[k];

        t_7[k] = f_1 * ls0_2[k]
                 - f_2 * ls1_2[k]
                 + pb_y[k] * lp_4[k];

        t_8[k] = f_5 * id0_6[k]
                 - f_6 * id1_6[k]
                 + pa_x[k] * kd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, id0_1, id0_8, id1_1, id1_8, kd_5, \
                         kd_12, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * id0_1[k]
                 - f_8 * id1_1[k]
                 + pa_y[k] * kd_5[k];

        t_10[k] = f_9 * id0_8[k]
                  - f_10 * id1_8[k]
                  + pa_x[k] * kd_12[k];

        t_11[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, id0_2, id0_10, id1_2, id1_10, \
                         kd_8, kd_16, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * id0_2[k]
                  - f_8 * id1_2[k]
                  + pa_z[k] * kd_8[k];

        t_13[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_14[k] = f_9 * id0_10[k]
                  - f_10 * id1_10[k]
                  + pa_x[k] * kd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, id0_3, id0_12, id1_3, id1_12, \
                         kd_11, kd_18, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_11 * id0_3[k]
                  - f_12 * id1_3[k]
                  + pa_y[k] * kd_11[k];

        t_16[k] = f_11 * id0_12[k]
                  - f_12 * id1_12[k]
                  + pa_x[k] * kd_18[k];

        t_17[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, id0_5, id0_14, id1_5, id1_14, \
                         kd_14, kd_22, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_11 * id0_5[k]
                  - f_12 * id1_5[k]
                  + pa_z[k] * kd_14[k];

        t_19[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_20[k] = f_11 * id0_14[k]
                  - f_12 * id1_14[k]
                  + pa_x[k] * kd_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, id0_7, id0_15, id1_7, id1_15, \
                         kd_17, kd_24, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * id0_7[k]
                  - f_10 * id1_7[k]
                  + pa_y[k] * kd_17[k];

        t_22[k] = f_7 * id0_15[k]
                  - f_8 * id1_15[k]
                  + pa_x[k] * kd_24[k];

        t_23[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, id0_9, id0_16, id1_9, id1_16, \
                         kd_20, kd_28, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * id0_9[k]
                  - f_10 * id1_9[k]
                  + pa_z[k] * kd_20[k];

        t_25[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_26[k] = f_7 * id0_16[k]
                  - f_8 * id1_16[k]
                  + pa_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, id0_11, id0_17, id1_11, id1_17, \
                         kd_23, kd_29, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * id0_11[k]
                  - f_6 * id1_11[k]
                  + pa_y[k] * kd_23[k];

        t_28[k] = f_3 * id0_17[k]
                  - f_4 * id1_17[k]
                  + pa_x[k] * kd_29[k];

        t_29[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_z, pb_y, id0_13, id0_26, id1_13, id1_26, \
                         kd_26, kd_30, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_13[k]
                  - f_6 * id1_13[k]
                  + pa_z[k] * kd_26[k];

        t_31[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_32[k] = f_3 * id0_26[k]
                  - f_4 * id1_26[k]
                  + pa_x[k] * kd_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, kp_12, ls0_11, ls0_12, \
                         ls1_11, ls1_12, lp_13, lp_14, lp_15, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_34[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_35[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_36[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, id0_17, id0_20, id1_17, id1_20, \
                         kd_34, kd_37, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * id0_17[k]
                  - f_4 * id1_17[k]
                  + pa_z[k] * kd_34[k];

        t_38[k] = f_5 * id0_20[k]
                  - f_6 * id1_20[k]
                  + pa_y[k] * kd_37[k];

        t_39[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, id0_18, id0_22, id1_18, id1_22, \
                         kd_36, kd_40, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * id0_18[k]
                  - f_8 * id1_18[k]
                  + pa_z[k] * kd_36[k];

        t_41[k] = f_9 * id0_22[k]
                  - f_10 * id1_22[k]
                  + pa_y[k] * kd_40[k];

        t_42[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, id0_19, id0_24, id1_19, id1_24, \
                         kd_39, kd_43, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * id0_19[k]
                  - f_12 * id1_19[k]
                  + pa_z[k] * kd_39[k];

        t_44[k] = f_11 * id0_24[k]
                  - f_12 * id1_24[k]
                  + pa_y[k] * kd_43[k];

        t_45[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_x[k] * lp_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, id0_21, id0_25, id1_21, id1_25, \
                         kd_42, kd_46, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * id0_21[k]
                  - f_10 * id1_21[k]
                  + pa_z[k] * kd_42[k];

        t_47[k] = f_7 * id0_25[k]
                  - f_8 * id1_25[k]
                  + pa_y[k] * kd_46[k];

        t_48[k] = f_1 * ls0_16[k]
                  - f_2 * ls1_16[k]
                  + pb_x[k] * lp_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_x, id0_23, id0_26, id1_23, id1_26, \
                         kd_45, kd_47, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * id0_23[k]
                  - f_6 * id1_23[k]
                  + pa_z[k] * kd_45[k];

        t_50[k] = f_3 * id0_26[k]
                  - f_4 * id1_26[k]
                  + pa_y[k] * kd_47[k];

        t_51[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_y[k] * lp_22[k];

        t_53[k] = f_0 * kp_20[k]
                  + f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_z[k] * lp_23[k];
    }
}

auto
compute_prim_ld_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / alpha;
    const auto f_6 = 2.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 2.0 / alpha;
    const auto f_10 = 2.0 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_1 = buffer.data(id0 + 1);
    const auto *id0_2 = buffer.data(id0 + 2);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_7 = buffer.data(id0 + 7);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_13 = buffer.data(id0 + 13);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_15 = buffer.data(id0 + 15);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_19 = buffer.data(id0 + 19);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_25 = buffer.data(id0 + 25);
    const auto *id0_26 = buffer.data(id0 + 26);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_9 = buffer.data(id1 + 9);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_28 = buffer.data(id1 + 28);
    const auto *id1_29 = buffer.data(id1 + 29);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_41 = buffer.data(id1 + 41);
    const auto *id1_43 = buffer.data(id1 + 43);
    const auto *id1_45 = buffer.data(id1 + 45);
    const auto *id1_54 = buffer.data(id1 + 54);
    const auto *id1_56 = buffer.data(id1 + 56);
    const auto *id1_59 = buffer.data(id1 + 59);
    const auto *id1_62 = buffer.data(id1 + 62);
    const auto *id1_63 = buffer.data(id1 + 63);
    const auto *id1_65 = buffer.data(id1 + 65);
    const auto *id1_66 = buffer.data(id1 + 66);
    const auto *id1_68 = buffer.data(id1 + 68);
    const auto *id1_69 = buffer.data(id1 + 69);
    const auto *id1_71 = buffer.data(id1 + 71);
    const auto *id1_74 = buffer.data(id1 + 74);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_86 = buffer.data(kd + 86);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, lp_0, lp_1, \
                         lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, id0_0, id0_4, id1_0, id1_10, kd_3, \
                         kd_10, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_y[k] * kd_3[k];

        t_4[k] = f_5 * id0_4[k]
                 - f_6 * id1_10[k]
                 + pa_x[k] * kd_10[k];

        t_5[k] = f_1 * ls0_1[k]
                 - f_2 * ls1_1[k]
                 + pb_z[k] * lp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, id0_0, id0_6, id1_0, id1_16, kd_6, \
                         kd_16, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_z[k] * kd_6[k];

        t_7[k] = f_1 * ls0_2[k]
                 - f_2 * ls1_2[k]
                 + pb_y[k] * lp_4[k];

        t_8[k] = f_5 * id0_6[k]
                 - f_6 * id1_16[k]
                 + pa_x[k] * kd_16[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, id0_1, id0_8, id1_3, id1_18, kd_9, \
                         kd_18, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * id0_1[k]
                 - f_8 * id1_3[k]
                 + pa_y[k] * kd_9[k];

        t_10[k] = f_9 * id0_8[k]
                  - f_10 * id1_18[k]
                  + pa_x[k] * kd_18[k];

        t_11[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, id0_2, id0_10, id1_6, id1_28, \
                         kd_14, kd_27, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * id0_2[k]
                  - f_8 * id1_6[k]
                  + pa_z[k] * kd_14[k];

        t_13[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_14[k] = f_9 * id0_10[k]
                  - f_10 * id1_28[k]
                  + pa_x[k] * kd_27[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, id0_3, id0_12, id1_9, id1_30, \
                         kd_17, kd_29, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_11 * id0_3[k]
                  - f_12 * id1_9[k]
                  + pa_y[k] * kd_17[k];

        t_16[k] = f_11 * id0_12[k]
                  - f_12 * id1_30[k]
                  + pa_x[k] * kd_29[k];

        t_17[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, id0_5, id0_14, id1_14, id1_43, \
                         kd_25, kd_41, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_11 * id0_5[k]
                  - f_12 * id1_14[k]
                  + pa_z[k] * kd_25[k];

        t_19[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_20[k] = f_11 * id0_14[k]
                  - f_12 * id1_43[k]
                  + pa_x[k] * kd_41[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, id0_7, id0_15, id1_17, id1_45, \
                         kd_28, kd_43, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * id0_7[k]
                  - f_10 * id1_17[k]
                  + pa_y[k] * kd_28[k];

        t_22[k] = f_7 * id0_15[k]
                  - f_8 * id1_45[k]
                  + pa_x[k] * kd_43[k];

        t_23[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, id0_9, id0_16, id1_26, id1_54, \
                         kd_39, kd_58, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * id0_9[k]
                  - f_10 * id1_26[k]
                  + pa_z[k] * kd_39[k];

        t_25[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_26[k] = f_7 * id0_16[k]
                  - f_8 * id1_54[k]
                  + pa_x[k] * kd_58[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, id0_11, id0_17, id1_29, id1_56, \
                         kd_42, kd_60, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * id0_11[k]
                  - f_6 * id1_29[k]
                  + pa_y[k] * kd_42[k];

        t_28[k] = f_3 * id0_17[k]
                  - f_4 * id1_56[k]
                  + pa_x[k] * kd_60[k];

        t_29[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_z, pb_y, id0_13, id0_26, id1_41, id1_74, \
                         kd_56, kd_66, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_13[k]
                  - f_6 * id1_41[k]
                  + pa_z[k] * kd_56[k];

        t_31[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_32[k] = f_3 * id0_26[k]
                  - f_4 * id1_74[k]
                  + pa_x[k] * kd_66[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, kp_12, ls0_11, ls0_12, \
                         ls1_11, ls1_12, lp_13, lp_14, lp_15, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_34[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_35[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_36[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, id0_17, id0_20, id1_56, id1_63, \
                         kd_71, kd_75, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * id0_17[k]
                  - f_4 * id1_56[k]
                  + pa_z[k] * kd_71[k];

        t_38[k] = f_5 * id0_20[k]
                  - f_6 * id1_63[k]
                  + pa_y[k] * kd_75[k];

        t_39[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, id0_18, id0_22, id1_59, id1_66, \
                         kd_74, kd_78, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * id0_18[k]
                  - f_8 * id1_59[k]
                  + pa_z[k] * kd_74[k];

        t_41[k] = f_9 * id0_22[k]
                  - f_10 * id1_66[k]
                  + pa_y[k] * kd_78[k];

        t_42[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, id0_19, id0_24, id1_62, id1_69, \
                         kd_77, kd_81, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * id0_19[k]
                  - f_12 * id1_62[k]
                  + pa_z[k] * kd_77[k];

        t_44[k] = f_11 * id0_24[k]
                  - f_12 * id1_69[k]
                  + pa_y[k] * kd_81[k];

        t_45[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_x[k] * lp_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, id0_21, id0_25, id1_65, id1_71, \
                         kd_80, kd_84, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * id0_21[k]
                  - f_10 * id1_65[k]
                  + pa_z[k] * kd_80[k];

        t_47[k] = f_7 * id0_25[k]
                  - f_8 * id1_71[k]
                  + pa_y[k] * kd_84[k];

        t_48[k] = f_1 * ls0_16[k]
                  - f_2 * ls1_16[k]
                  + pb_x[k] * lp_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_x, id0_23, id0_26, id1_68, id1_74, \
                         kd_83, kd_86, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * id0_23[k]
                  - f_6 * id1_68[k]
                  + pa_z[k] * kd_83[k];

        t_50[k] = f_3 * id0_26[k]
                  - f_4 * id1_74[k]
                  + pa_y[k] * kd_86[k];

        t_51[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_y[k] * lp_22[k];

        t_53[k] = f_0 * kp_20[k]
                  + f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_z[k] * lp_23[k];
    }
}

auto
compute_prim_ld_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.5 / alpha;
    const auto f_7 = 2.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / alpha;
    const auto f_11 = 2.0 * beta / (alpha * p);
    const auto f_12 = 1.5 / alpha;
    const auto f_13 = 1.5 * beta / (alpha * p);

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

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_26 = buffer.data(id0 + 26);
    const auto *id0_28 = buffer.data(id0 + 28);
    const auto *id0_29 = buffer.data(id0 + 29);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_35 = buffer.data(id0 + 35);
    const auto *id0_36 = buffer.data(id0 + 36);
    const auto *id0_37 = buffer.data(id0 + 37);
    const auto *id0_38 = buffer.data(id0 + 38);
    const auto *id0_41 = buffer.data(id0 + 41);
    const auto *id0_43 = buffer.data(id0 + 43);
    const auto *id0_45 = buffer.data(id0 + 45);
    const auto *id0_48 = buffer.data(id0 + 48);
    const auto *id0_49 = buffer.data(id0 + 49);
    const auto *id0_51 = buffer.data(id0 + 51);
    const auto *id0_52 = buffer.data(id0 + 52);
    const auto *id0_54 = buffer.data(id0 + 54);
    const auto *id0_56 = buffer.data(id0 + 56);
    const auto *id0_59 = buffer.data(id0 + 59);
    const auto *id0_62 = buffer.data(id0 + 62);
    const auto *id0_63 = buffer.data(id0 + 63);
    const auto *id0_65 = buffer.data(id0 + 65);
    const auto *id0_66 = buffer.data(id0 + 66);
    const auto *id0_68 = buffer.data(id0 + 68);
    const auto *id0_69 = buffer.data(id0 + 69);
    const auto *id0_71 = buffer.data(id0 + 71);
    const auto *id0_74 = buffer.data(id0 + 74);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_7 = buffer.data(id1 + 7);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_13 = buffer.data(id1 + 13);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_19 = buffer.data(id1 + 19);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_25 = buffer.data(id1 + 25);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_27 = buffer.data(id1 + 27);
    const auto *id1_29 = buffer.data(id1 + 29);
    const auto *id1_31 = buffer.data(id1 + 31);
    const auto *id1_32 = buffer.data(id1 + 32);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_34 = buffer.data(id1 + 34);
    const auto *id1_35 = buffer.data(id1 + 35);
    const auto *id1_37 = buffer.data(id1 + 37);
    const auto *id1_39 = buffer.data(id1 + 39);
    const auto *id1_41 = buffer.data(id1 + 41);
    const auto *id1_44 = buffer.data(id1 + 44);
    const auto *id1_45 = buffer.data(id1 + 45);
    const auto *id1_47 = buffer.data(id1 + 47);
    const auto *id1_48 = buffer.data(id1 + 48);
    const auto *id1_50 = buffer.data(id1 + 50);
    const auto *id1_51 = buffer.data(id1 + 51);
    const auto *id1_53 = buffer.data(id1 + 53);
    const auto *id1_56 = buffer.data(id1 + 56);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, kp_0, kd_0, ls0_0, ls1_0, \
                         lp_0, lp_1, lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];

        t_3[k] = pa_y[k] * kd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, id0_0, id1_0, kp_1, kp_2, \
                         kd_0, kd_1, kd_2, kd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * kp_1[k]
                 + pa_y[k] * kd_1[k];

        t_5[k] = pa_y[k] * kd_2[k];

        t_6[k] = pa_z[k] * kd_0[k];

        t_7[k] = pa_z[k] * kd_1[k];

        t_8[k] = f_3 * kp_2[k]
                 + pa_z[k] * kd_2[k];

        t_9[k] = f_4 * id0_0[k]
                 - f_5 * id1_0[k]
                 + pa_y[k] * kd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, id0_10, id1_8, kd_4, \
                         kd_6, kd_8, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * id0_10[k]
                  - f_7 * id1_8[k]
                  + pa_x[k] * kd_8[k];

        t_11[k] = f_1 * ls0_1[k]
                  - f_2 * ls1_1[k]
                  + pb_z[k] * lp_3[k];

        t_12[k] = pa_z[k] * kd_4[k];

        t_13[k] = pa_y[k] * kd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, id0_0, id0_16, id1_0, id1_12, \
                         kd_5, kd_12, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * id0_0[k]
                  - f_5 * id1_0[k]
                  + pa_z[k] * kd_5[k];

        t_15[k] = f_1 * ls0_2[k]
                  - f_2 * ls1_2[k]
                  + pb_y[k] * lp_4[k];

        t_16[k] = f_6 * id0_16[k]
                  - f_7 * id1_12[k]
                  + pa_x[k] * kd_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, id0_3, id0_18, id1_3, id1_14, \
                         kd_7, kd_14, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * id0_3[k]
                  - f_9 * id1_3[k]
                  + pa_y[k] * kd_7[k];

        t_18[k] = f_10 * id0_18[k]
                  - f_11 * id1_14[k]
                  + pa_x[k] * kd_14[k];

        t_19[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, kp_3, kp_4, kd_7, kd_8, \
                         kd_9, kd_11, kd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * kd_7[k];

        t_21[k] = pa_z[k] * kd_8[k];

        t_22[k] = f_3 * kp_3[k]
                  + pa_z[k] * kd_9[k];

        t_23[k] = f_3 * kp_4[k]
                  + pa_y[k] * kd_11[k];

        t_24[k] = pa_y[k] * kd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, id0_6, id0_28, id1_5, id1_19, \
                         kd_10, kd_19, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * id0_6[k]
                  - f_9 * id1_5[k]
                  + pa_z[k] * kd_10[k];

        t_26[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_27[k] = f_10 * id0_28[k]
                  - f_11 * id1_19[k]
                  + pa_x[k] * kd_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_z, id0_9, id0_30, id1_7, id1_21, \
                         kd_13, kd_21, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_12 * id0_9[k]
                  - f_13 * id1_7[k]
                  + pa_y[k] * kd_13[k];

        t_29[k] = f_12 * id0_30[k]
                  - f_13 * id1_21[k]
                  + pa_x[k] * kd_21[k];

        t_30[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, id0_14, id1_10, kp_5, kd_13, \
                         kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_z[k] * kd_13[k];

        t_32[k] = pa_z[k] * kd_14[k];

        t_33[k] = f_3 * kp_5[k]
                  + pa_z[k] * kd_15[k];

        t_34[k] = f_4 * id0_14[k]
                  - f_5 * id1_10[k]
                  + pa_y[k] * kd_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pa_y, id0_36, id0_37, id1_24, id1_25, \
                         kp_6, kd_18, kd_19, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_12 * id0_36[k]
                  - f_13 * id1_24[k]
                  + pa_x[k] * kd_24[k];

        t_36[k] = f_12 * id0_37[k]
                  - f_13 * id1_25[k]
                  + pa_x[k] * kd_25[k];

        t_37[k] = f_3 * kp_6[k]
                  + pa_y[k] * kd_18[k];

        t_38[k] = pa_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_z, pb_y, id0_14, id0_43, id1_10, id1_29, \
                         kd_17, kd_29, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * id0_14[k]
                  - f_13 * id1_10[k]
                  + pa_z[k] * kd_17[k];

        t_40[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_41[k] = f_12 * id0_43[k]
                  - f_13 * id1_29[k]
                  + pa_x[k] * kd_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pa_y, pb_z, id0_17, id0_45, id1_13, id1_31, \
                         kd_20, kd_31, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_10 * id0_17[k]
                  - f_11 * id1_13[k]
                  + pa_y[k] * kd_20[k];

        t_43[k] = f_8 * id0_45[k]
                  - f_9 * id1_31[k]
                  + pa_x[k] * kd_31[k];

        t_44[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pa_z, id0_23, id1_16, kp_7, kd_20, \
                         kd_21, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * kd_20[k];

        t_46[k] = pa_z[k] * kd_21[k];

        t_47[k] = f_3 * kp_7[k]
                  + pa_z[k] * kd_22[k];

        t_48[k] = f_8 * id0_23[k]
                  - f_9 * id1_16[k]
                  + pa_y[k] * kd_23[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_x, pa_y, id0_26, id0_48, id0_49, id1_17, id1_32, \
                         id1_33, kd_26, kd_34, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * id0_48[k]
                  - f_9 * id1_32[k]
                  + pa_x[k] * kd_34[k];

        t_50[k] = f_8 * id0_49[k]
                  - f_9 * id1_33[k]
                  + pa_x[k] * kd_35[k];

        t_51[k] = f_4 * id0_26[k]
                  - f_5 * id1_17[k]
                  + pa_y[k] * kd_26[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, id0_51, id0_52, id1_34, id1_35, \
                         kp_8, kd_28, kd_29, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * id0_51[k]
                  - f_9 * id1_34[k]
                  + pa_x[k] * kd_37[k];

        t_53[k] = f_8 * id0_52[k]
                  - f_9 * id1_35[k]
                  + pa_x[k] * kd_38[k];

        t_54[k] = f_3 * kp_8[k]
                  + pa_y[k] * kd_28[k];

        t_55[k] = pa_y[k] * kd_29[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_x, pa_z, pb_y, id0_26, id0_54, id1_17, id1_37, \
                         kd_27, kd_42, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_10 * id0_26[k]
                  - f_11 * id1_17[k]
                  + pa_z[k] * kd_27[k];

        t_57[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_58[k] = f_8 * id0_54[k]
                  - f_9 * id1_37[k]
                  + pa_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pa_y, pb_z, id0_29, id0_56, id1_20, id1_39, \
                         kd_30, kd_44, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_6 * id0_29[k]
                  - f_7 * id1_20[k]
                  + pa_y[k] * kd_30[k];

        t_60[k] = f_4 * id0_56[k]
                  - f_5 * id1_39[k]
                  + pa_x[k] * kd_44[k];

        t_61[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, id0_35, id1_23, kp_9, kd_30, \
                         kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * kd_30[k];

        t_63[k] = pa_z[k] * kd_31[k];

        t_64[k] = f_3 * kp_9[k]
                  + pa_z[k] * kd_32[k];

        t_65[k] = f_12 * id0_35[k]
                  - f_13 * id1_23[k]
                  + pa_y[k] * kd_33[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_x, pa_y, id0_38, id0_62, id0_63, id1_26, id1_44, \
                         id1_45, kd_36, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * id0_62[k]
                  - f_5 * id1_44[k]
                  + pa_x[k] * kd_45[k];

        t_67[k] = f_4 * id0_63[k]
                  - f_5 * id1_45[k]
                  + pa_x[k] * kd_46[k];

        t_68[k] = f_8 * id0_38[k]
                  - f_9 * id1_26[k]
                  + pa_y[k] * kd_36[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_x, pa_y, id0_41, id0_65, id0_66, id1_27, id1_47, \
                         id1_48, kd_39, kd_47, kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_4 * id0_65[k]
                  - f_5 * id1_47[k]
                  + pa_x[k] * kd_47[k];

        t_70[k] = f_4 * id0_66[k]
                  - f_5 * id1_48[k]
                  + pa_x[k] * kd_48[k];

        t_71[k] = f_4 * id0_41[k]
                  - f_5 * id1_27[k]
                  + pa_y[k] * kd_39[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_x, pa_y, id0_68, id0_69, id1_50, id1_51, \
                         kp_10, kd_41, kd_42, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_4 * id0_68[k]
                  - f_5 * id1_50[k]
                  + pa_x[k] * kd_49[k];

        t_73[k] = f_4 * id0_69[k]
                  - f_5 * id1_51[k]
                  + pa_x[k] * kd_50[k];

        t_74[k] = f_3 * kp_10[k]
                  + pa_y[k] * kd_41[k];

        t_75[k] = pa_y[k] * kd_42[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_x, pa_z, pb_y, id0_41, id0_74, id1_27, id1_56, \
                         kd_40, kd_52, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_6 * id0_41[k]
                  - f_7 * id1_27[k]
                  + pa_z[k] * kd_40[k];

        t_77[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_78[k] = f_4 * id0_74[k]
                  - f_5 * id1_56[k]
                  + pa_x[k] * kd_52[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_x, pa_z, kp_11, kp_14, kp_15, kd_43, \
                         kd_53, kd_54, kd_58, kd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * kp_11[k]
                  + pa_x[k] * kd_53[k];

        t_80[k] = pa_x[k] * kd_54[k];

        t_81[k] = pa_z[k] * kd_43[k];

        t_82[k] = f_3 * kp_14[k]
                  + pa_x[k] * kd_58[k];

        t_83[k] = f_3 * kp_15[k]
                  + pa_x[k] * kd_61[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, kp_16, kp_17, kp_18, kd_64, kd_67, \
                         kd_72, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_3 * kp_16[k]
                  + pa_x[k] * kd_64[k];

        t_85[k] = f_3 * kp_17[k]
                  + pa_x[k] * kd_67[k];

        t_86[k] = f_3 * kp_18[k]
                  + pa_x[k] * kd_72[k];

        t_87[k] = pa_x[k] * kd_74[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_z, pb_x, pb_y, pb_z, kp_12, kd_53, ls0_11, \
                         ls1_11, lp_13, lp_14, lp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_89[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_90[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_91[k] = pa_z[k] * kd_53[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pb_x, id0_56, id1_39, kp_13, kd_54, \
                         kd_55, kd_56, ls0_12, ls1_12, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * kd_54[k];

        t_93[k] = f_3 * kp_13[k]
                  + pa_z[k] * kd_55[k];

        t_94[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];

        t_95[k] = f_4 * id0_56[k]
                  - f_5 * id1_39[k]
                  + pa_z[k] * kd_56[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_y, pa_z, pb_x, id0_59, id0_63, id1_41, id1_45, \
                         kd_59, kd_60, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_6 * id0_63[k]
                  - f_7 * id1_45[k]
                  + pa_y[k] * kd_60[k];

        t_97[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];

        t_98[k] = f_8 * id0_59[k]
                  - f_9 * id1_41[k]
                  + pa_z[k] * kd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_y, pa_z, pb_x, id0_62, id0_66, id1_44, id1_48, \
                         kd_62, kd_63, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * id0_66[k]
                  - f_11 * id1_48[k]
                  + pa_y[k] * kd_63[k];

        t_100[k] = f_1 * ls0_14[k]
                   - f_2 * ls1_14[k]
                   + pb_x[k] * lp_18[k];

        t_101[k] = f_12 * id0_62[k]
                   - f_13 * id1_44[k]
                   + pa_z[k] * kd_62[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pa_z, pb_x, id0_65, id0_69, id1_47, \
                         id1_51, kd_65, kd_66, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * id0_69[k]
                   - f_13 * id1_51[k]
                   + pa_y[k] * kd_66[k];

        t_103[k] = f_1 * ls0_15[k]
                   - f_2 * ls1_15[k]
                   + pb_x[k] * lp_19[k];

        t_104[k] = f_10 * id0_65[k]
                   - f_11 * id1_47[k]
                   + pa_z[k] * kd_65[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_y, pa_z, pb_x, id0_68, id0_71, id1_50, \
                         id1_53, kd_68, kd_69, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_8 * id0_71[k]
                   - f_9 * id1_53[k]
                   + pa_y[k] * kd_69[k];

        t_106[k] = f_1 * ls0_16[k]
                   - f_2 * ls1_16[k]
                   + pb_x[k] * lp_20[k];

        t_107[k] = f_6 * id0_68[k]
                   - f_7 * id1_50[k]
                   + pa_z[k] * kd_68[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pb_x, id0_74, id1_56, kp_19, kd_71, \
                         kd_73, kd_74, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_4 * id0_74[k]
                   - f_5 * id1_56[k]
                   + pa_y[k] * kd_71[k];

        t_109[k] = f_3 * kp_19[k]
                   + pa_y[k] * kd_73[k];

        t_110[k] = pa_y[k] * kd_74[k];

        t_111[k] = f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_112, t_113, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_y[k] * lp_22[k];

        t_113[k] = f_0 * kp_20[k]
                   + f_1 * ls0_17[k]
                   - f_2 * ls1_17[k]
                   + pb_z[k] * lp_23[k];
    }
}

auto
compute_prim_ld_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / alpha;
    const auto f_6 = 2.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 2.0 / alpha;
    const auto f_10 = 2.0 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_1 = buffer.data(id0 + 1);
    const auto *id0_2 = buffer.data(id0 + 2);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_7 = buffer.data(id0 + 7);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_13 = buffer.data(id0 + 13);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_15 = buffer.data(id0 + 15);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_19 = buffer.data(id0 + 19);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_25 = buffer.data(id0 + 25);
    const auto *id0_26 = buffer.data(id0 + 26);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_4 = buffer.data(id1 + 4);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_11 = buffer.data(id1 + 11);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_22 = buffer.data(id1 + 22);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_28 = buffer.data(id1 + 28);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_31 = buffer.data(id1 + 31);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_34 = buffer.data(id1 + 34);
    const auto *id1_36 = buffer.data(id1 + 36);
    const auto *id1_37 = buffer.data(id1 + 37);
    const auto *id1_38 = buffer.data(id1 + 38);
    const auto *id1_41 = buffer.data(id1 + 41);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, lp_0, lp_1, \
                         lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, id0_0, id0_4, id1_0, id1_6, kd_3, \
                         kd_6, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_y[k] * kd_3[k];

        t_4[k] = f_5 * id0_4[k]
                 - f_6 * id1_6[k]
                 + pa_x[k] * kd_6[k];

        t_5[k] = f_1 * ls0_1[k]
                 - f_2 * ls1_1[k]
                 + pb_z[k] * lp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, id0_0, id0_6, id1_0, id1_10, kd_4, \
                         kd_10, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_z[k] * kd_4[k];

        t_7[k] = f_1 * ls0_2[k]
                 - f_2 * ls1_2[k]
                 + pb_y[k] * lp_4[k];

        t_8[k] = f_5 * id0_6[k]
                 - f_6 * id1_10[k]
                 + pa_x[k] * kd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, id0_1, id0_8, id1_3, id1_12, kd_5, \
                         kd_12, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * id0_1[k]
                 - f_8 * id1_3[k]
                 + pa_y[k] * kd_5[k];

        t_10[k] = f_9 * id0_8[k]
                  - f_10 * id1_12[k]
                  + pa_x[k] * kd_12[k];

        t_11[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, id0_2, id0_10, id1_4, id1_16, \
                         kd_8, kd_16, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * id0_2[k]
                  - f_8 * id1_4[k]
                  + pa_z[k] * kd_8[k];

        t_13[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_14[k] = f_9 * id0_10[k]
                  - f_10 * id1_16[k]
                  + pa_x[k] * kd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, id0_3, id0_12, id1_5, id1_18, \
                         kd_11, kd_18, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_11 * id0_3[k]
                  - f_12 * id1_5[k]
                  + pa_y[k] * kd_11[k];

        t_16[k] = f_11 * id0_12[k]
                  - f_12 * id1_18[k]
                  + pa_x[k] * kd_18[k];

        t_17[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, id0_5, id0_14, id1_8, id1_22, \
                         kd_14, kd_22, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_11 * id0_5[k]
                  - f_12 * id1_8[k]
                  + pa_z[k] * kd_14[k];

        t_19[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_20[k] = f_11 * id0_14[k]
                  - f_12 * id1_22[k]
                  + pa_x[k] * kd_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, id0_7, id0_15, id1_11, id1_23, \
                         kd_17, kd_24, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * id0_7[k]
                  - f_10 * id1_11[k]
                  + pa_y[k] * kd_17[k];

        t_22[k] = f_7 * id0_15[k]
                  - f_8 * id1_23[k]
                  + pa_x[k] * kd_24[k];

        t_23[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, id0_9, id0_16, id1_14, id1_24, \
                         kd_20, kd_28, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * id0_9[k]
                  - f_10 * id1_14[k]
                  + pa_z[k] * kd_20[k];

        t_25[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_26[k] = f_7 * id0_16[k]
                  - f_8 * id1_24[k]
                  + pa_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, id0_11, id0_17, id1_17, id1_26, \
                         kd_23, kd_29, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * id0_11[k]
                  - f_6 * id1_17[k]
                  + pa_y[k] * kd_23[k];

        t_28[k] = f_3 * id0_17[k]
                  - f_4 * id1_26[k]
                  + pa_x[k] * kd_29[k];

        t_29[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_z, pb_y, id0_13, id0_26, id1_20, id1_41, \
                         kd_26, kd_30, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_13[k]
                  - f_6 * id1_20[k]
                  + pa_z[k] * kd_26[k];

        t_31[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_32[k] = f_3 * id0_26[k]
                  - f_4 * id1_41[k]
                  + pa_x[k] * kd_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, kp_12, ls0_11, ls0_12, \
                         ls1_11, ls1_12, lp_13, lp_14, lp_15, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_34[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_35[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_36[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, id0_17, id0_20, id1_26, id1_31, \
                         kd_34, kd_37, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * id0_17[k]
                  - f_4 * id1_26[k]
                  + pa_z[k] * kd_34[k];

        t_38[k] = f_5 * id0_20[k]
                  - f_6 * id1_31[k]
                  + pa_y[k] * kd_37[k];

        t_39[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, id0_18, id0_22, id1_28, id1_34, \
                         kd_36, kd_40, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * id0_18[k]
                  - f_8 * id1_28[k]
                  + pa_z[k] * kd_36[k];

        t_41[k] = f_9 * id0_22[k]
                  - f_10 * id1_34[k]
                  + pa_y[k] * kd_40[k];

        t_42[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, id0_19, id0_24, id1_30, id1_37, \
                         kd_39, kd_43, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * id0_19[k]
                  - f_12 * id1_30[k]
                  + pa_z[k] * kd_39[k];

        t_44[k] = f_11 * id0_24[k]
                  - f_12 * id1_37[k]
                  + pa_y[k] * kd_43[k];

        t_45[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_x[k] * lp_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, id0_21, id0_25, id1_33, id1_38, \
                         kd_42, kd_46, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * id0_21[k]
                  - f_10 * id1_33[k]
                  + pa_z[k] * kd_42[k];

        t_47[k] = f_7 * id0_25[k]
                  - f_8 * id1_38[k]
                  + pa_y[k] * kd_46[k];

        t_48[k] = f_1 * ls0_16[k]
                  - f_2 * ls1_16[k]
                  + pb_x[k] * lp_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_x, id0_23, id0_26, id1_36, id1_41, \
                         kd_45, kd_47, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * id0_23[k]
                  - f_6 * id1_36[k]
                  + pa_z[k] * kd_45[k];

        t_50[k] = f_3 * id0_26[k]
                  - f_4 * id1_41[k]
                  + pa_y[k] * kd_47[k];

        t_51[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_y[k] * lp_22[k];

        t_53[k] = f_0 * kp_20[k]
                  + f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_z[k] * lp_23[k];
    }
}

auto
compute_prim_ld_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / alpha;
    const auto f_6 = 2.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 2.0 / alpha;
    const auto f_10 = 2.0 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_26 = buffer.data(id0 + 26);
    const auto *id0_28 = buffer.data(id0 + 28);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_31 = buffer.data(id0 + 31);
    const auto *id0_33 = buffer.data(id0 + 33);
    const auto *id0_34 = buffer.data(id0 + 34);
    const auto *id0_36 = buffer.data(id0 + 36);
    const auto *id0_37 = buffer.data(id0 + 37);
    const auto *id0_38 = buffer.data(id0 + 38);
    const auto *id0_41 = buffer.data(id0 + 41);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_4 = buffer.data(id1 + 4);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_7 = buffer.data(id1 + 7);
    const auto *id1_9 = buffer.data(id1 + 9);
    const auto *id1_11 = buffer.data(id1 + 11);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_13 = buffer.data(id1 + 13);
    const auto *id1_15 = buffer.data(id1 + 15);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_19 = buffer.data(id1 + 19);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_25 = buffer.data(id1 + 25);
    const auto *id1_27 = buffer.data(id1 + 27);
    const auto *id1_29 = buffer.data(id1 + 29);
    const auto *id1_32 = buffer.data(id1 + 32);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_35 = buffer.data(id1 + 35);
    const auto *id1_36 = buffer.data(id1 + 36);
    const auto *id1_38 = buffer.data(id1 + 38);
    const auto *id1_39 = buffer.data(id1 + 39);
    const auto *id1_41 = buffer.data(id1 + 41);
    const auto *id1_44 = buffer.data(id1 + 44);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_50 = buffer.data(kd + 50);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, lp_0, lp_1, \
                         lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, id0_0, id0_6, id1_0, id1_7, kd_3, \
                         kd_7, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_y[k] * kd_3[k];

        t_4[k] = f_5 * id0_6[k]
                 - f_6 * id1_7[k]
                 + pa_x[k] * kd_7[k];

        t_5[k] = f_1 * ls0_1[k]
                 - f_2 * ls1_1[k]
                 + pb_z[k] * lp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, id0_0, id0_10, id1_0, id1_11, kd_4, \
                         kd_11, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_z[k] * kd_4[k];

        t_7[k] = f_1 * ls0_2[k]
                 - f_2 * ls1_2[k]
                 + pb_y[k] * lp_4[k];

        t_8[k] = f_5 * id0_10[k]
                 - f_6 * id1_11[k]
                 + pa_x[k] * kd_11[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, id0_3, id0_12, id1_3, id1_13, \
                         kd_6, kd_13, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * id0_3[k]
                 - f_8 * id1_3[k]
                 + pa_y[k] * kd_6[k];

        t_10[k] = f_9 * id0_12[k]
                  - f_10 * id1_13[k]
                  + pa_x[k] * kd_13[k];

        t_11[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, id0_4, id0_16, id1_4, id1_17, \
                         kd_9, kd_17, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * id0_4[k]
                  - f_8 * id1_4[k]
                  + pa_z[k] * kd_9[k];

        t_13[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_14[k] = f_9 * id0_16[k]
                  - f_10 * id1_17[k]
                  + pa_x[k] * kd_17[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, id0_5, id0_18, id1_6, id1_19, \
                         kd_12, kd_19, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_11 * id0_5[k]
                  - f_12 * id1_6[k]
                  + pa_y[k] * kd_12[k];

        t_16[k] = f_11 * id0_18[k]
                  - f_12 * id1_19[k]
                  + pa_x[k] * kd_19[k];

        t_17[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, id0_8, id0_22, id1_9, id1_23, \
                         kd_15, kd_23, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_11 * id0_8[k]
                  - f_12 * id1_9[k]
                  + pa_z[k] * kd_15[k];

        t_19[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_20[k] = f_11 * id0_22[k]
                  - f_12 * id1_23[k]
                  + pa_x[k] * kd_23[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, id0_11, id0_23, id1_12, id1_24, \
                         kd_18, kd_25, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * id0_11[k]
                  - f_10 * id1_12[k]
                  + pa_y[k] * kd_18[k];

        t_22[k] = f_7 * id0_23[k]
                  - f_8 * id1_24[k]
                  + pa_x[k] * kd_25[k];

        t_23[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, id0_14, id0_24, id1_15, id1_25, \
                         kd_21, kd_29, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * id0_14[k]
                  - f_10 * id1_15[k]
                  + pa_z[k] * kd_21[k];

        t_25[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_26[k] = f_7 * id0_24[k]
                  - f_8 * id1_25[k]
                  + pa_x[k] * kd_29[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, id0_17, id0_26, id1_18, id1_27, \
                         kd_24, kd_30, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * id0_17[k]
                  - f_6 * id1_18[k]
                  + pa_y[k] * kd_24[k];

        t_28[k] = f_3 * id0_26[k]
                  - f_4 * id1_27[k]
                  + pa_x[k] * kd_30[k];

        t_29[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_z, pb_y, id0_20, id0_41, id1_21, id1_44, \
                         kd_27, kd_31, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_20[k]
                  - f_6 * id1_21[k]
                  + pa_z[k] * kd_27[k];

        t_31[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_32[k] = f_3 * id0_41[k]
                  - f_4 * id1_44[k]
                  + pa_x[k] * kd_31[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, kp_12, ls0_11, ls0_12, \
                         ls1_11, ls1_12, lp_13, lp_14, lp_15, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_34[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_35[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_36[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, id0_26, id0_31, id1_27, id1_33, \
                         kd_35, kd_39, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * id0_26[k]
                  - f_4 * id1_27[k]
                  + pa_z[k] * kd_35[k];

        t_38[k] = f_5 * id0_31[k]
                  - f_6 * id1_33[k]
                  + pa_y[k] * kd_39[k];

        t_39[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, id0_28, id0_34, id1_29, id1_36, \
                         kd_38, kd_42, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * id0_28[k]
                  - f_8 * id1_29[k]
                  + pa_z[k] * kd_38[k];

        t_41[k] = f_9 * id0_34[k]
                  - f_10 * id1_36[k]
                  + pa_y[k] * kd_42[k];

        t_42[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, id0_30, id0_37, id1_32, id1_39, \
                         kd_41, kd_45, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * id0_30[k]
                  - f_12 * id1_32[k]
                  + pa_z[k] * kd_41[k];

        t_44[k] = f_11 * id0_37[k]
                  - f_12 * id1_39[k]
                  + pa_y[k] * kd_45[k];

        t_45[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_x[k] * lp_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, id0_33, id0_38, id1_35, id1_41, \
                         kd_44, kd_48, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * id0_33[k]
                  - f_10 * id1_35[k]
                  + pa_z[k] * kd_44[k];

        t_47[k] = f_7 * id0_38[k]
                  - f_8 * id1_41[k]
                  + pa_y[k] * kd_48[k];

        t_48[k] = f_1 * ls0_16[k]
                  - f_2 * ls1_16[k]
                  + pb_x[k] * lp_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_x, id0_36, id0_41, id1_38, id1_44, \
                         kd_47, kd_50, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * id0_36[k]
                  - f_6 * id1_38[k]
                  + pa_z[k] * kd_47[k];

        t_50[k] = f_3 * id0_41[k]
                  - f_4 * id1_44[k]
                  + pa_y[k] * kd_50[k];

        t_51[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_y[k] * lp_22[k];

        t_53[k] = f_0 * kp_20[k]
                  + f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_z[k] * lp_23[k];
    }
}

auto
compute_prim_ld_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / alpha;
    const auto f_6 = 2.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 2.0 / alpha;
    const auto f_10 = 2.0 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_7 = buffer.data(id0 + 7);
    const auto *id0_9 = buffer.data(id0 + 9);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_13 = buffer.data(id0 + 13);
    const auto *id0_15 = buffer.data(id0 + 15);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_19 = buffer.data(id0 + 19);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_25 = buffer.data(id0 + 25);
    const auto *id0_27 = buffer.data(id0 + 27);
    const auto *id0_29 = buffer.data(id0 + 29);
    const auto *id0_32 = buffer.data(id0 + 32);
    const auto *id0_33 = buffer.data(id0 + 33);
    const auto *id0_35 = buffer.data(id0 + 35);
    const auto *id0_36 = buffer.data(id0 + 36);
    const auto *id0_38 = buffer.data(id0 + 38);
    const auto *id0_39 = buffer.data(id0 + 39);
    const auto *id0_41 = buffer.data(id0 + 41);
    const auto *id0_44 = buffer.data(id0 + 44);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_4 = buffer.data(id1 + 4);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_11 = buffer.data(id1 + 11);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_22 = buffer.data(id1 + 22);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_28 = buffer.data(id1 + 28);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_31 = buffer.data(id1 + 31);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_34 = buffer.data(id1 + 34);
    const auto *id1_36 = buffer.data(id1 + 36);
    const auto *id1_37 = buffer.data(id1 + 37);
    const auto *id1_38 = buffer.data(id1 + 38);
    const auto *id1_41 = buffer.data(id1 + 41);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, lp_0, lp_1, \
                         lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, id0_0, id0_7, id1_0, id1_6, kd_3, \
                         kd_6, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_y[k] * kd_3[k];

        t_4[k] = f_5 * id0_7[k]
                 - f_6 * id1_6[k]
                 + pa_x[k] * kd_6[k];

        t_5[k] = f_1 * ls0_1[k]
                 - f_2 * ls1_1[k]
                 + pb_z[k] * lp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, id0_0, id0_11, id1_0, id1_10, kd_4, \
                         kd_10, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_z[k] * kd_4[k];

        t_7[k] = f_1 * ls0_2[k]
                 - f_2 * ls1_2[k]
                 + pb_y[k] * lp_4[k];

        t_8[k] = f_5 * id0_11[k]
                 - f_6 * id1_10[k]
                 + pa_x[k] * kd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, id0_3, id0_13, id1_3, id1_12, \
                         kd_5, kd_12, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * id0_3[k]
                 - f_8 * id1_3[k]
                 + pa_y[k] * kd_5[k];

        t_10[k] = f_9 * id0_13[k]
                  - f_10 * id1_12[k]
                  + pa_x[k] * kd_12[k];

        t_11[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, id0_4, id0_17, id1_4, id1_16, \
                         kd_8, kd_16, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * id0_4[k]
                  - f_8 * id1_4[k]
                  + pa_z[k] * kd_8[k];

        t_13[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_14[k] = f_9 * id0_17[k]
                  - f_10 * id1_16[k]
                  + pa_x[k] * kd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, id0_6, id0_19, id1_5, id1_18, \
                         kd_11, kd_18, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_11 * id0_6[k]
                  - f_12 * id1_5[k]
                  + pa_y[k] * kd_11[k];

        t_16[k] = f_11 * id0_19[k]
                  - f_12 * id1_18[k]
                  + pa_x[k] * kd_18[k];

        t_17[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, id0_9, id0_23, id1_8, id1_22, \
                         kd_14, kd_22, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_11 * id0_9[k]
                  - f_12 * id1_8[k]
                  + pa_z[k] * kd_14[k];

        t_19[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_20[k] = f_11 * id0_23[k]
                  - f_12 * id1_22[k]
                  + pa_x[k] * kd_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, id0_12, id0_24, id1_11, id1_23, \
                         kd_17, kd_24, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * id0_12[k]
                  - f_10 * id1_11[k]
                  + pa_y[k] * kd_17[k];

        t_22[k] = f_7 * id0_24[k]
                  - f_8 * id1_23[k]
                  + pa_x[k] * kd_24[k];

        t_23[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, id0_15, id0_25, id1_14, id1_24, \
                         kd_20, kd_28, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * id0_15[k]
                  - f_10 * id1_14[k]
                  + pa_z[k] * kd_20[k];

        t_25[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_26[k] = f_7 * id0_25[k]
                  - f_8 * id1_24[k]
                  + pa_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, id0_18, id0_27, id1_17, id1_26, \
                         kd_23, kd_29, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * id0_18[k]
                  - f_6 * id1_17[k]
                  + pa_y[k] * kd_23[k];

        t_28[k] = f_3 * id0_27[k]
                  - f_4 * id1_26[k]
                  + pa_x[k] * kd_29[k];

        t_29[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_z, pb_y, id0_21, id0_44, id1_20, id1_41, \
                         kd_26, kd_30, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_21[k]
                  - f_6 * id1_20[k]
                  + pa_z[k] * kd_26[k];

        t_31[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_32[k] = f_3 * id0_44[k]
                  - f_4 * id1_41[k]
                  + pa_x[k] * kd_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, kp_12, ls0_11, ls0_12, \
                         ls1_11, ls1_12, lp_13, lp_14, lp_15, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_34[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_35[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_36[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, id0_27, id0_33, id1_26, id1_31, \
                         kd_34, kd_37, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * id0_27[k]
                  - f_4 * id1_26[k]
                  + pa_z[k] * kd_34[k];

        t_38[k] = f_5 * id0_33[k]
                  - f_6 * id1_31[k]
                  + pa_y[k] * kd_37[k];

        t_39[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, id0_29, id0_36, id1_28, id1_34, \
                         kd_36, kd_40, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * id0_29[k]
                  - f_8 * id1_28[k]
                  + pa_z[k] * kd_36[k];

        t_41[k] = f_9 * id0_36[k]
                  - f_10 * id1_34[k]
                  + pa_y[k] * kd_40[k];

        t_42[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, id0_32, id0_39, id1_30, id1_37, \
                         kd_39, kd_43, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * id0_32[k]
                  - f_12 * id1_30[k]
                  + pa_z[k] * kd_39[k];

        t_44[k] = f_11 * id0_39[k]
                  - f_12 * id1_37[k]
                  + pa_y[k] * kd_43[k];

        t_45[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_x[k] * lp_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, id0_35, id0_41, id1_33, id1_38, \
                         kd_42, kd_46, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * id0_35[k]
                  - f_10 * id1_33[k]
                  + pa_z[k] * kd_42[k];

        t_47[k] = f_7 * id0_41[k]
                  - f_8 * id1_38[k]
                  + pa_y[k] * kd_46[k];

        t_48[k] = f_1 * ls0_16[k]
                  - f_2 * ls1_16[k]
                  + pb_x[k] * lp_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_x, id0_38, id0_44, id1_36, id1_41, \
                         kd_45, kd_47, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * id0_38[k]
                  - f_6 * id1_36[k]
                  + pa_z[k] * kd_45[k];

        t_50[k] = f_3 * id0_44[k]
                  - f_4 * id1_41[k]
                  + pa_y[k] * kd_47[k];

        t_51[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_y[k] * lp_22[k];

        t_53[k] = f_0 * kp_20[k]
                  + f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_z[k] * lp_23[k];
    }
}

auto
compute_prim_ld_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t id0, const size_t id1,
                                     const size_t kp, const size_t kd, const size_t ls0,
                                     const size_t ls1, const size_t lp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / alpha;
    const auto f_6 = 2.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 2.0 / alpha;
    const auto f_10 = 2.0 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_4 = buffer.data(id0 + 4);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_8 = buffer.data(id0 + 8);
    const auto *id0_10 = buffer.data(id0 + 10);
    const auto *id0_11 = buffer.data(id0 + 11);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_14 = buffer.data(id0 + 14);
    const auto *id0_16 = buffer.data(id0 + 16);
    const auto *id0_17 = buffer.data(id0 + 17);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_20 = buffer.data(id0 + 20);
    const auto *id0_22 = buffer.data(id0 + 22);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_24 = buffer.data(id0 + 24);
    const auto *id0_26 = buffer.data(id0 + 26);
    const auto *id0_28 = buffer.data(id0 + 28);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_31 = buffer.data(id0 + 31);
    const auto *id0_33 = buffer.data(id0 + 33);
    const auto *id0_34 = buffer.data(id0 + 34);
    const auto *id0_36 = buffer.data(id0 + 36);
    const auto *id0_37 = buffer.data(id0 + 37);
    const auto *id0_38 = buffer.data(id0 + 38);
    const auto *id0_41 = buffer.data(id0 + 41);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_4 = buffer.data(id1 + 4);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_8 = buffer.data(id1 + 8);
    const auto *id1_10 = buffer.data(id1 + 10);
    const auto *id1_11 = buffer.data(id1 + 11);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_14 = buffer.data(id1 + 14);
    const auto *id1_16 = buffer.data(id1 + 16);
    const auto *id1_17 = buffer.data(id1 + 17);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_20 = buffer.data(id1 + 20);
    const auto *id1_22 = buffer.data(id1 + 22);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_24 = buffer.data(id1 + 24);
    const auto *id1_26 = buffer.data(id1 + 26);
    const auto *id1_28 = buffer.data(id1 + 28);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_31 = buffer.data(id1 + 31);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_34 = buffer.data(id1 + 34);
    const auto *id1_36 = buffer.data(id1 + 36);
    const auto *id1_37 = buffer.data(id1 + 37);
    const auto *id1_38 = buffer.data(id1 + 38);
    const auto *id1_41 = buffer.data(id1 + 41);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_20 = buffer.data(kp + 20);

    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_1 = buffer.data(ls0 + 1);
    const auto *ls0_2 = buffer.data(ls0 + 2);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_4 = buffer.data(ls0 + 4);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_7 = buffer.data(ls0 + 7);
    const auto *ls0_8 = buffer.data(ls0 + 8);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_11 = buffer.data(ls0 + 11);
    const auto *ls0_12 = buffer.data(ls0 + 12);
    const auto *ls0_13 = buffer.data(ls0 + 13);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_16 = buffer.data(ls0 + 16);
    const auto *ls0_17 = buffer.data(ls0 + 17);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_1 = buffer.data(ls1 + 1);
    const auto *ls1_2 = buffer.data(ls1 + 2);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_4 = buffer.data(ls1 + 4);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_7 = buffer.data(ls1 + 7);
    const auto *ls1_8 = buffer.data(ls1 + 8);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_11 = buffer.data(ls1 + 11);
    const auto *ls1_12 = buffer.data(ls1 + 12);
    const auto *ls1_13 = buffer.data(ls1 + 13);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_16 = buffer.data(ls1 + 16);
    const auto *ls1_17 = buffer.data(ls1 + 17);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, kp_0, ls0_0, ls1_0, lp_0, lp_1, \
                         lp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_0[k]
                 + f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_x[k] * lp_0[k];

        t_1[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_y[k] * lp_1[k];

        t_2[k] = f_1 * ls0_0[k]
                 - f_2 * ls1_0[k]
                 + pb_z[k] * lp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, id0_0, id0_6, id1_0, id1_6, kd_3, \
                         kd_6, ls0_1, ls1_1, lp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_y[k] * kd_3[k];

        t_4[k] = f_5 * id0_6[k]
                 - f_6 * id1_6[k]
                 + pa_x[k] * kd_6[k];

        t_5[k] = f_1 * ls0_1[k]
                 - f_2 * ls1_1[k]
                 + pb_z[k] * lp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, id0_0, id0_10, id1_0, id1_10, kd_4, \
                         kd_10, ls0_2, ls1_2, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pa_z[k] * kd_4[k];

        t_7[k] = f_1 * ls0_2[k]
                 - f_2 * ls1_2[k]
                 + pb_y[k] * lp_4[k];

        t_8[k] = f_5 * id0_10[k]
                 - f_6 * id1_10[k]
                 + pa_x[k] * kd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, id0_3, id0_12, id1_3, id1_12, \
                         kd_5, kd_12, ls0_3, ls1_3, lp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * id0_3[k]
                 - f_8 * id1_3[k]
                 + pa_y[k] * kd_5[k];

        t_10[k] = f_9 * id0_12[k]
                  - f_10 * id1_12[k]
                  + pa_x[k] * kd_12[k];

        t_11[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, id0_4, id0_16, id1_4, id1_16, \
                         kd_8, kd_16, ls0_4, ls1_4, lp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * id0_4[k]
                  - f_8 * id1_4[k]
                  + pa_z[k] * kd_8[k];

        t_13[k] = f_1 * ls0_4[k]
                  - f_2 * ls1_4[k]
                  + pb_y[k] * lp_6[k];

        t_14[k] = f_9 * id0_16[k]
                  - f_10 * id1_16[k]
                  + pa_x[k] * kd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, id0_5, id0_18, id1_5, id1_18, \
                         kd_11, kd_18, ls0_5, ls1_5, lp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_11 * id0_5[k]
                  - f_12 * id1_5[k]
                  + pa_y[k] * kd_11[k];

        t_16[k] = f_11 * id0_18[k]
                  - f_12 * id1_18[k]
                  + pa_x[k] * kd_18[k];

        t_17[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_z[k] * lp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, id0_8, id0_22, id1_8, id1_22, \
                         kd_14, kd_22, ls0_6, ls1_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_11 * id0_8[k]
                  - f_12 * id1_8[k]
                  + pa_z[k] * kd_14[k];

        t_19[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_y[k] * lp_8[k];

        t_20[k] = f_11 * id0_22[k]
                  - f_12 * id1_22[k]
                  + pa_x[k] * kd_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, id0_11, id0_23, id1_11, id1_23, \
                         kd_17, kd_24, ls0_7, ls1_7, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * id0_11[k]
                  - f_10 * id1_11[k]
                  + pa_y[k] * kd_17[k];

        t_22[k] = f_7 * id0_23[k]
                  - f_8 * id1_23[k]
                  + pa_x[k] * kd_24[k];

        t_23[k] = f_1 * ls0_7[k]
                  - f_2 * ls1_7[k]
                  + pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, id0_14, id0_24, id1_14, id1_24, \
                         kd_20, kd_28, ls0_8, ls1_8, lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * id0_14[k]
                  - f_10 * id1_14[k]
                  + pa_z[k] * kd_20[k];

        t_25[k] = f_1 * ls0_8[k]
                  - f_2 * ls1_8[k]
                  + pb_y[k] * lp_10[k];

        t_26[k] = f_7 * id0_24[k]
                  - f_8 * id1_24[k]
                  + pa_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, id0_17, id0_26, id1_17, id1_26, \
                         kd_23, kd_29, ls0_9, ls1_9, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * id0_17[k]
                  - f_6 * id1_17[k]
                  + pa_y[k] * kd_23[k];

        t_28[k] = f_3 * id0_26[k]
                  - f_4 * id1_26[k]
                  + pa_x[k] * kd_29[k];

        t_29[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_z[k] * lp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_z, pb_y, id0_20, id0_41, id1_20, id1_41, \
                         kd_26, kd_30, ls0_10, ls1_10, lp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_20[k]
                  - f_6 * id1_20[k]
                  + pa_z[k] * kd_26[k];

        t_31[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_y[k] * lp_12[k];

        t_32[k] = f_3 * id0_41[k]
                  - f_4 * id1_41[k]
                  + pa_x[k] * kd_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, kp_12, ls0_11, ls0_12, \
                         ls1_11, ls1_12, lp_13, lp_14, lp_15, lp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_x[k] * lp_13[k];

        t_34[k] = f_0 * kp_12[k]
                  + f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_y[k] * lp_14[k];

        t_35[k] = f_1 * ls0_11[k]
                  - f_2 * ls1_11[k]
                  + pb_z[k] * lp_15[k];

        t_36[k] = f_1 * ls0_12[k]
                  - f_2 * ls1_12[k]
                  + pb_x[k] * lp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, id0_26, id0_31, id1_26, id1_31, \
                         kd_34, kd_37, ls0_13, ls1_13, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * id0_26[k]
                  - f_4 * id1_26[k]
                  + pa_z[k] * kd_34[k];

        t_38[k] = f_5 * id0_31[k]
                  - f_6 * id1_31[k]
                  + pa_y[k] * kd_37[k];

        t_39[k] = f_1 * ls0_13[k]
                  - f_2 * ls1_13[k]
                  + pb_x[k] * lp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, id0_28, id0_34, id1_28, id1_34, \
                         kd_36, kd_40, ls0_14, ls1_14, lp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * id0_28[k]
                  - f_8 * id1_28[k]
                  + pa_z[k] * kd_36[k];

        t_41[k] = f_9 * id0_34[k]
                  - f_10 * id1_34[k]
                  + pa_y[k] * kd_40[k];

        t_42[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_x[k] * lp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, id0_30, id0_37, id1_30, id1_37, \
                         kd_39, kd_43, ls0_15, ls1_15, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_11 * id0_30[k]
                  - f_12 * id1_30[k]
                  + pa_z[k] * kd_39[k];

        t_44[k] = f_11 * id0_37[k]
                  - f_12 * id1_37[k]
                  + pa_y[k] * kd_43[k];

        t_45[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_x[k] * lp_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, id0_33, id0_38, id1_33, id1_38, \
                         kd_42, kd_46, ls0_16, ls1_16, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * id0_33[k]
                  - f_10 * id1_33[k]
                  + pa_z[k] * kd_42[k];

        t_47[k] = f_7 * id0_38[k]
                  - f_8 * id1_38[k]
                  + pa_y[k] * kd_46[k];

        t_48[k] = f_1 * ls0_16[k]
                  - f_2 * ls1_16[k]
                  + pb_x[k] * lp_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_x, id0_36, id0_41, id1_36, id1_41, \
                         kd_45, kd_47, ls0_17, ls1_17, lp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_5 * id0_36[k]
                  - f_6 * id1_36[k]
                  + pa_z[k] * kd_45[k];

        t_50[k] = f_3 * id0_41[k]
                  - f_4 * id1_41[k]
                  + pa_y[k] * kd_47[k];

        t_51[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_x[k] * lp_21[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, kp_20, ls0_17, ls1_17, lp_22, \
                         lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_y[k] * lp_22[k];

        t_53[k] = f_0 * kp_20[k]
                  + f_1 * ls0_17[k]
                  - f_2 * ls1_17[k]
                  + pb_z[k] * lp_23[k];
    }
}

}  // namespace simdt2ceri
