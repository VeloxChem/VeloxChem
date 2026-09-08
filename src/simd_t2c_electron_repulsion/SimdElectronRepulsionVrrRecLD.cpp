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
    const auto *id0_6 = buffer.data(id0 + 6);
    const auto *id0_12 = buffer.data(id0 + 12);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_35 = buffer.data(id0 + 35);
    const auto *id0_36 = buffer.data(id0 + 36);
    const auto *id0_39 = buffer.data(id0 + 39);
    const auto *id0_48 = buffer.data(id0 + 48);
    const auto *id0_54 = buffer.data(id0 + 54);
    const auto *id0_59 = buffer.data(id0 + 59);
    const auto *id0_60 = buffer.data(id0 + 60);
    const auto *id0_63 = buffer.data(id0 + 63);
    const auto *id0_72 = buffer.data(id0 + 72);
    const auto *id0_75 = buffer.data(id0 + 75);
    const auto *id0_77 = buffer.data(id0 + 77);
    const auto *id0_78 = buffer.data(id0 + 78);
    const auto *id0_84 = buffer.data(id0 + 84);
    const auto *id0_89 = buffer.data(id0 + 89);
    const auto *id0_93 = buffer.data(id0 + 93);
    const auto *id0_105 = buffer.data(id0 + 105);
    const auto *id0_107 = buffer.data(id0 + 107);
    const auto *id0_111 = buffer.data(id0 + 111);
    const auto *id0_113 = buffer.data(id0 + 113);
    const auto *id0_125 = buffer.data(id0 + 125);
    const auto *id0_129 = buffer.data(id0 + 129);
    const auto *id0_135 = buffer.data(id0 + 135);
    const auto *id0_141 = buffer.data(id0 + 141);
    const auto *id0_143 = buffer.data(id0 + 143);
    const auto *id0_147 = buffer.data(id0 + 147);
    const auto *id0_149 = buffer.data(id0 + 149);
    const auto *id0_153 = buffer.data(id0 + 153);
    const auto *id0_155 = buffer.data(id0 + 155);
    const auto *id0_161 = buffer.data(id0 + 161);
    const auto *id0_167 = buffer.data(id0 + 167);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_6 = buffer.data(id1 + 6);
    const auto *id1_12 = buffer.data(id1 + 12);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_35 = buffer.data(id1 + 35);
    const auto *id1_36 = buffer.data(id1 + 36);
    const auto *id1_39 = buffer.data(id1 + 39);
    const auto *id1_48 = buffer.data(id1 + 48);
    const auto *id1_54 = buffer.data(id1 + 54);
    const auto *id1_59 = buffer.data(id1 + 59);
    const auto *id1_60 = buffer.data(id1 + 60);
    const auto *id1_63 = buffer.data(id1 + 63);
    const auto *id1_72 = buffer.data(id1 + 72);
    const auto *id1_75 = buffer.data(id1 + 75);
    const auto *id1_77 = buffer.data(id1 + 77);
    const auto *id1_78 = buffer.data(id1 + 78);
    const auto *id1_84 = buffer.data(id1 + 84);
    const auto *id1_89 = buffer.data(id1 + 89);
    const auto *id1_93 = buffer.data(id1 + 93);
    const auto *id1_105 = buffer.data(id1 + 105);
    const auto *id1_107 = buffer.data(id1 + 107);
    const auto *id1_111 = buffer.data(id1 + 111);
    const auto *id1_113 = buffer.data(id1 + 113);
    const auto *id1_125 = buffer.data(id1 + 125);
    const auto *id1_129 = buffer.data(id1 + 129);
    const auto *id1_135 = buffer.data(id1 + 135);
    const auto *id1_141 = buffer.data(id1 + 141);
    const auto *id1_143 = buffer.data(id1 + 143);
    const auto *id1_147 = buffer.data(id1 + 147);
    const auto *id1_149 = buffer.data(id1 + 149);
    const auto *id1_153 = buffer.data(id1 + 153);
    const auto *id1_155 = buffer.data(id1 + 155);
    const auto *id1_161 = buffer.data(id1 + 161);
    const auto *id1_167 = buffer.data(id1 + 167);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_215 = buffer.data(kd + 215);

    const auto *ls0_0 = buffer.data(ls0 + 0);
    const auto *ls0_3 = buffer.data(ls0 + 3);
    const auto *ls0_5 = buffer.data(ls0 + 5);
    const auto *ls0_6 = buffer.data(ls0 + 6);
    const auto *ls0_9 = buffer.data(ls0 + 9);
    const auto *ls0_10 = buffer.data(ls0 + 10);
    const auto *ls0_14 = buffer.data(ls0 + 14);
    const auto *ls0_15 = buffer.data(ls0 + 15);
    const auto *ls0_20 = buffer.data(ls0 + 20);
    const auto *ls0_21 = buffer.data(ls0 + 21);
    const auto *ls0_27 = buffer.data(ls0 + 27);
    const auto *ls0_36 = buffer.data(ls0 + 36);
    const auto *ls0_38 = buffer.data(ls0 + 38);
    const auto *ls0_39 = buffer.data(ls0 + 39);
    const auto *ls0_40 = buffer.data(ls0 + 40);
    const auto *ls0_41 = buffer.data(ls0 + 41);
    const auto *ls0_42 = buffer.data(ls0 + 42);
    const auto *ls0_44 = buffer.data(ls0 + 44);

    const auto *ls1_0 = buffer.data(ls1 + 0);
    const auto *ls1_3 = buffer.data(ls1 + 3);
    const auto *ls1_5 = buffer.data(ls1 + 5);
    const auto *ls1_6 = buffer.data(ls1 + 6);
    const auto *ls1_9 = buffer.data(ls1 + 9);
    const auto *ls1_10 = buffer.data(ls1 + 10);
    const auto *ls1_14 = buffer.data(ls1 + 14);
    const auto *ls1_15 = buffer.data(ls1 + 15);
    const auto *ls1_20 = buffer.data(ls1 + 20);
    const auto *ls1_21 = buffer.data(ls1 + 21);
    const auto *ls1_27 = buffer.data(ls1 + 27);
    const auto *ls1_36 = buffer.data(ls1 + 36);
    const auto *ls1_38 = buffer.data(ls1 + 38);
    const auto *ls1_39 = buffer.data(ls1 + 39);
    const auto *ls1_40 = buffer.data(ls1 + 40);
    const auto *ls1_41 = buffer.data(ls1 + 41);
    const auto *ls1_42 = buffer.data(ls1 + 42);
    const auto *ls1_44 = buffer.data(ls1 + 44);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_60 = buffer.data(lp + 60);
    const auto *lp_61 = buffer.data(lp + 61);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_63 = buffer.data(lp + 63);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_65 = buffer.data(lp + 65);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_81 = buffer.data(lp + 81);
    const auto *lp_82 = buffer.data(lp + 82);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_84 = buffer.data(lp + 84);
    const auto *lp_85 = buffer.data(lp + 85);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_100 = buffer.data(lp + 100);
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_105 = buffer.data(lp + 105);
    const auto *lp_107 = buffer.data(lp + 107);
    const auto *lp_108 = buffer.data(lp + 108);
    const auto *lp_109 = buffer.data(lp + 109);
    const auto *lp_110 = buffer.data(lp + 110);
    const auto *lp_112 = buffer.data(lp + 112);
    const auto *lp_113 = buffer.data(lp + 113);
    const auto *lp_114 = buffer.data(lp + 114);
    const auto *lp_115 = buffer.data(lp + 115);
    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_117 = buffer.data(lp + 117);
    const auto *lp_118 = buffer.data(lp + 118);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_120 = buffer.data(lp + 120);
    const auto *lp_121 = buffer.data(lp + 121);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_123 = buffer.data(lp + 123);
    const auto *lp_124 = buffer.data(lp + 124);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_126 = buffer.data(lp + 126);
    const auto *lp_127 = buffer.data(lp + 127);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_131 = buffer.data(lp + 131);
    const auto *lp_132 = buffer.data(lp + 132);
    const auto *lp_133 = buffer.data(lp + 133);
    const auto *lp_134 = buffer.data(lp + 134);

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

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, kp_1, kp_4, kd_0, \
                         kd_3, kd_5, lp_3, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * kd_0[k];

        t_7[k] = f_3 * kp_4[k]
                 + pb_x[k] * lp_4[k];

        t_8[k] = pb_z[k] * lp_3[k];

        t_9[k] = f_4 * kp_1[k]
                 + pa_y[k] * kd_3[k];

        t_10[k] = pb_z[k] * lp_4[k];

        t_11[k] = pa_y[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, kp_2, kp_8, \
                         kd_0, kd_3, kd_5, lp_6, lp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * kd_0[k];

        t_13[k] = pb_y[k] * lp_6[k];

        t_14[k] = f_3 * kp_8[k]
                  + pb_x[k] * lp_8[k];

        t_15[k] = pa_z[k] * kd_3[k];

        t_16[k] = pb_y[k] * lp_8[k];

        t_17[k] = f_4 * kp_2[k]
                  + pa_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, id0_0, id1_0, kp_10, kd_6, lp_9, \
                         lp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * id0_0[k]
                  - f_6 * id1_0[k]
                  + pa_y[k] * kd_6[k];

        t_19[k] = f_7 * kp_10[k]
                  + pb_x[k] * lp_10[k];

        t_20[k] = pb_z[k] * lp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, id0_21, id1_21, kd_12, \
                         kd_21, ls0_3, ls1_3, lp_10, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * id0_21[k]
                  - f_9 * id1_21[k]
                  + pa_x[k] * kd_21[k];

        t_22[k] = pb_z[k] * lp_10[k];

        t_23[k] = f_1 * ls0_3[k]
                  - f_2 * ls1_3[k]
                  + pb_z[k] * lp_11[k];

        t_24[k] = pa_y[k] * kd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, kp_8, kd_7, kd_9, \
                         kd_14, kd_17, lp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * kd_7[k];

        t_26[k] = pa_y[k] * kd_14[k];

        t_27[k] = pa_z[k] * kd_9[k];

        t_28[k] = f_10 * kp_8[k]
                  + pb_y[k] * lp_14[k];

        t_29[k] = pa_y[k] * kd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, id0_0, id1_0, kp_17, kd_12, \
                         ls0_5, ls1_5, lp_15, lp_16, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * id0_0[k]
                  - f_6 * id1_0[k]
                  + pa_z[k] * kd_12[k];

        t_31[k] = pb_y[k] * lp_15[k];

        t_32[k] = f_7 * kp_17[k]
                  + pb_x[k] * lp_17[k];

        t_33[k] = f_1 * ls0_5[k]
                  - f_2 * ls1_5[k]
                  + pb_y[k] * lp_16[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, id0_6, id0_35, id1_6, id1_35, \
                         kd_18, kd_35, lp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * lp_17[k];

        t_35[k] = f_8 * id0_35[k]
                  - f_9 * id1_35[k]
                  + pa_x[k] * kd_35[k];

        t_36[k] = f_11 * id0_6[k]
                  - f_12 * id1_6[k]
                  + pa_y[k] * kd_18[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, id0_39, id1_39, kp_19, \
                         kd_39, lp_18, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * kp_19[k]
                  + pb_x[k] * lp_19[k];

        t_38[k] = pb_z[k] * lp_18[k];

        t_39[k] = f_14 * id0_39[k]
                  - f_15 * id1_39[k]
                  + pa_x[k] * kd_39[k];

        t_40[k] = pb_z[k] * lp_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, kp_23, kd_18, kd_19, \
                         kd_21, ls0_6, ls1_6, lp_20, lp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ls0_6[k]
                  - f_2 * ls1_6[k]
                  + pb_z[k] * lp_20[k];

        t_42[k] = pa_z[k] * kd_18[k];

        t_43[k] = pa_z[k] * kd_19[k];

        t_44[k] = f_13 * kp_23[k]
                  + pb_x[k] * lp_23[k];

        t_45[k] = pa_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, kp_11, kp_14, kp_25, \
                         kd_23, kd_30, lp_23, lp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * kp_14[k]
                  + pb_y[k] * lp_23[k];

        t_47[k] = f_4 * kp_11[k]
                  + pa_z[k] * kd_23[k];

        t_48[k] = pa_y[k] * kd_30[k];

        t_49[k] = f_13 * kp_25[k]
                  + pb_x[k] * lp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, kp_16, kp_17, kd_32, kd_33, \
                         kd_35, lp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * kd_32[k];

        t_51[k] = f_4 * kp_16[k]
                  + pa_y[k] * kd_33[k];

        t_52[k] = f_10 * kp_17[k]
                  + pb_y[k] * lp_26[k];

        t_53[k] = pa_y[k] * kd_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, id0_12, id1_12, kp_29, \
                         kd_30, ls0_9, ls1_9, lp_27, lp_28, lp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * id0_12[k]
                  - f_12 * id1_12[k]
                  + pa_z[k] * kd_30[k];

        t_55[k] = pb_y[k] * lp_27[k];

        t_56[k] = f_13 * kp_29[k]
                  + pb_x[k] * lp_29[k];

        t_57[k] = f_1 * ls0_9[k]
                  - f_2 * ls1_9[k]
                  + pb_y[k] * lp_28[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, id0_18, id0_59, id1_18, id1_59, \
                         kd_36, kd_59, lp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * lp_29[k];

        t_59[k] = f_14 * id0_59[k]
                  - f_15 * id1_59[k]
                  + pa_x[k] * kd_59[k];

        t_60[k] = f_16 * id0_18[k]
                  - f_17 * id1_18[k]
                  + pa_y[k] * kd_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, id0_63, id1_63, kp_31, \
                         kd_63, lp_30, lp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_18 * kp_31[k]
                  + pb_x[k] * lp_31[k];

        t_62[k] = pb_z[k] * lp_30[k];

        t_63[k] = f_16 * id0_63[k]
                  - f_17 * id1_63[k]
                  + pa_x[k] * kd_63[k];

        t_64[k] = pb_z[k] * lp_31[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_z, kp_35, kd_36, kd_37, \
                         kd_39, ls0_10, ls1_10, lp_32, lp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ls0_10[k]
                  - f_2 * ls1_10[k]
                  + pb_z[k] * lp_32[k];

        t_66[k] = pa_z[k] * kd_36[k];

        t_67[k] = pa_z[k] * kd_37[k];

        t_68[k] = f_18 * kp_35[k]
                  + pb_x[k] * lp_35[k];

        t_69[k] = pa_z[k] * kd_39[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_y, id0_30, id1_30, kp_20, kp_23, \
                         kd_41, kd_48, lp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_19 * kp_23[k]
                  + pb_y[k] * lp_35[k];

        t_71[k] = f_4 * kp_20[k]
                  + pa_z[k] * kd_41[k];

        t_72[k] = f_5 * id0_30[k]
                  - f_6 * id1_30[k]
                  + pa_y[k] * kd_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, pb_y, id0_75, id1_75, kp_26, \
                         kp_37, kp_38, kd_75, lp_37, lp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_18 * kp_37[k]
                  + pb_x[k] * lp_37[k];

        t_74[k] = f_18 * kp_38[k]
                  + pb_x[k] * lp_38[k];

        t_75[k] = f_16 * id0_75[k]
                  - f_17 * id1_75[k]
                  + pa_x[k] * kd_75[k];

        t_76[k] = f_4 * kp_26[k]
                  + pb_y[k] * lp_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, id0_77, id1_77, kp_40, \
                         kd_54, kd_56, kd_77, lp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_16 * id0_77[k]
                  - f_17 * id1_77[k]
                  + pa_x[k] * kd_77[k];

        t_78[k] = pa_y[k] * kd_54[k];

        t_79[k] = f_18 * kp_40[k]
                  + pb_x[k] * lp_40[k];

        t_80[k] = pa_y[k] * kd_56[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, id0_30, id1_30, kp_28, \
                         kp_29, kd_54, kd_57, kd_59, lp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * kp_28[k]
                  + pa_y[k] * kd_57[k];

        t_82[k] = f_10 * kp_29[k]
                  + pb_y[k] * lp_41[k];

        t_83[k] = pa_y[k] * kd_59[k];

        t_84[k] = f_16 * id0_30[k]
                  - f_17 * id1_30[k]
                  + pa_z[k] * kd_54[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, kp_44, ls0_14, ls1_14, lp_42, \
                         lp_43, lp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * lp_42[k];

        t_86[k] = f_18 * kp_44[k]
                  + pb_x[k] * lp_44[k];

        t_87[k] = f_1 * ls0_14[k]
                  - f_2 * ls1_14[k]
                  + pb_y[k] * lp_43[k];

        t_88[k] = pb_y[k] * lp_44[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pa_y, pb_x, id0_36, id0_89, id1_36, id1_89, \
                         kp_46, kd_60, kd_89, lp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_16 * id0_89[k]
                  - f_17 * id1_89[k]
                  + pa_x[k] * kd_89[k];

        t_90[k] = f_14 * id0_36[k]
                  - f_15 * id1_36[k]
                  + pa_y[k] * kd_60[k];

        t_91[k] = f_19 * kp_46[k]
                  + pb_x[k] * lp_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, id0_93, id1_93, kd_93, ls0_15, \
                         ls1_15, lp_45, lp_46, lp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * lp_45[k];

        t_93[k] = f_11 * id0_93[k]
                  - f_12 * id1_93[k]
                  + pa_x[k] * kd_93[k];

        t_94[k] = pb_z[k] * lp_46[k];

        t_95[k] = f_1 * ls0_15[k]
                  - f_2 * ls1_15[k]
                  + pb_z[k] * lp_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, kp_35, kp_50, kd_60, \
                         kd_61, kd_63, lp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * kd_60[k];

        t_97[k] = pa_z[k] * kd_61[k];

        t_98[k] = f_19 * kp_50[k]
                  + pb_x[k] * lp_50[k];

        t_99[k] = pa_z[k] * kd_63[k];

        t_100[k] = f_18 * kp_35[k]
                   + pb_y[k] * lp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pa_z, pb_x, id0_48, id1_48, kp_32, \
                         kp_52, kp_53, kd_65, kd_72, lp_52, lp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_4 * kp_32[k]
                   + pa_z[k] * kd_65[k];

        t_102[k] = f_11 * id0_48[k]
                   - f_12 * id1_48[k]
                   + pa_y[k] * kd_72[k];

        t_103[k] = f_19 * kp_52[k]
                   + pb_x[k] * lp_52[k];

        t_104[k] = f_19 * kp_53[k]
                   + pb_x[k] * lp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_y, id0_105, id0_107, id1_105, id1_107, \
                         kp_38, kd_105, kd_107, lp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_11 * id0_105[k]
                   - f_12 * id1_105[k]
                   + pa_x[k] * kd_105[k];

        t_106[k] = f_19 * kp_38[k]
                   + pb_y[k] * lp_53[k];

        t_107[k] = f_11 * id0_107[k]
                   - f_12 * id1_107[k]
                   + pa_x[k] * kd_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_x, id0_54, id1_54, kp_55, kp_56, kd_78, \
                         lp_55, lp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * id0_54[k]
                   - f_6 * id1_54[k]
                   + pa_y[k] * kd_78[k];

        t_109[k] = f_19 * kp_55[k]
                   + pb_x[k] * lp_55[k];

        t_110[k] = f_19 * kp_56[k]
                   + pb_x[k] * lp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, pb_y, id0_111, id0_113, \
                         id1_111, id1_113, kp_41, kd_84, kd_111, kd_113, \
                         lp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_11 * id0_111[k]
                   - f_12 * id1_111[k]
                   + pa_x[k] * kd_111[k];

        t_112[k] = f_4 * kp_41[k]
                   + pb_y[k] * lp_56[k];

        t_113[k] = f_11 * id0_113[k]
                   - f_12 * id1_113[k]
                   + pa_x[k] * kd_113[k];

        t_114[k] = pa_y[k] * kd_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, kp_43, kp_44, \
                         kp_58, kd_86, kd_87, kd_89, lp_58, lp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_19 * kp_58[k]
                   + pb_x[k] * lp_58[k];

        t_116[k] = pa_y[k] * kd_86[k];

        t_117[k] = f_4 * kp_43[k]
                   + pa_y[k] * kd_87[k];

        t_118[k] = f_10 * kp_44[k]
                   + pb_y[k] * lp_59[k];

        t_119[k] = pa_y[k] * kd_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, pb_x, pb_y, id0_54, id1_54, kp_62, \
                         kd_84, ls0_20, ls1_20, lp_60, lp_61, lp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_14 * id0_54[k]
                   - f_15 * id1_54[k]
                   + pa_z[k] * kd_84[k];

        t_121[k] = pb_y[k] * lp_60[k];

        t_122[k] = f_19 * kp_62[k]
                   + pb_x[k] * lp_62[k];

        t_123[k] = f_1 * ls0_20[k]
                   - f_2 * ls1_20[k]
                   + pb_y[k] * lp_61[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pa_y, pb_y, id0_60, id0_125, id1_60, \
                         id1_125, kd_90, kd_125, lp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * lp_62[k];

        t_125[k] = f_11 * id0_125[k]
                   - f_12 * id1_125[k]
                   + pa_x[k] * kd_125[k];

        t_126[k] = f_8 * id0_60[k]
                   - f_9 * id1_60[k]
                   + pa_y[k] * kd_90[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pb_x, pb_z, id0_129, id1_129, \
                         kp_64, kd_129, lp_63, lp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_4 * kp_64[k]
                   + pb_x[k] * lp_64[k];

        t_128[k] = pb_z[k] * lp_63[k];

        t_129[k] = f_5 * id0_129[k]
                   - f_6 * id1_129[k]
                   + pa_x[k] * kd_129[k];

        t_130[k] = pb_z[k] * lp_64[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_z, pb_x, pb_z, kp_68, kd_90, \
                         kd_91, kd_93, ls0_21, ls1_21, lp_65, lp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_1 * ls0_21[k]
                   - f_2 * ls1_21[k]
                   + pb_z[k] * lp_65[k];

        t_132[k] = pa_z[k] * kd_90[k];

        t_133[k] = pa_z[k] * kd_91[k];

        t_134[k] = f_4 * kp_68[k]
                   + pb_x[k] * lp_68[k];

        t_135[k] = pa_z[k] * kd_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_y, pa_z, pb_y, id0_72, id1_72, kp_47, kp_50, \
                         kd_95, kd_102, lp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_13 * kp_50[k]
                   + pb_y[k] * lp_68[k];

        t_137[k] = f_4 * kp_47[k]
                   + pa_z[k] * kd_95[k];

        t_138[k] = f_16 * id0_72[k]
                   - f_17 * id1_72[k]
                   + pa_y[k] * kd_102[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, pb_y, id0_141, id1_141, \
                         kp_53, kp_70, kp_71, kd_141, lp_70, lp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_4 * kp_70[k]
                   + pb_x[k] * lp_70[k];

        t_140[k] = f_4 * kp_71[k]
                   + pb_x[k] * lp_71[k];

        t_141[k] = f_5 * id0_141[k]
                   - f_6 * id1_141[k]
                   + pa_x[k] * kd_141[k];

        t_142[k] = f_18 * kp_53[k]
                   + pb_y[k] * lp_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pa_y, pb_x, id0_78, id0_143, id1_78, \
                         id1_143, kp_73, kd_108, kd_143, lp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_5 * id0_143[k]
                   - f_6 * id1_143[k]
                   + pa_x[k] * kd_143[k];

        t_144[k] = f_11 * id0_78[k]
                   - f_12 * id1_78[k]
                   + pa_y[k] * kd_108[k];

        t_145[k] = f_4 * kp_73[k]
                   + pb_x[k] * lp_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_x, pb_y, id0_147, id0_149, \
                         id1_147, id1_149, kp_56, kp_74, kd_147, kd_149, \
                         lp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_4 * kp_74[k]
                   + pb_x[k] * lp_74[k];

        t_147[k] = f_5 * id0_147[k]
                   - f_6 * id1_147[k]
                   + pa_x[k] * kd_147[k];

        t_148[k] = f_19 * kp_56[k]
                   + pb_y[k] * lp_74[k];

        t_149[k] = f_5 * id0_149[k]
                   - f_6 * id1_149[k]
                   + pa_x[k] * kd_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_y, pb_x, id0_84, id1_84, kp_76, kp_77, \
                         kd_114, lp_76, lp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * id0_84[k]
                   - f_6 * id1_84[k]
                   + pa_y[k] * kd_114[k];

        t_151[k] = f_4 * kp_76[k]
                   + pb_x[k] * lp_76[k];

        t_152[k] = f_4 * kp_77[k]
                   + pb_x[k] * lp_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_x, pa_y, pb_y, id0_153, id0_155, \
                         id1_153, id1_155, kp_59, kd_120, kd_153, kd_155, \
                         lp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_5 * id0_153[k]
                   - f_6 * id1_153[k]
                   + pa_x[k] * kd_153[k];

        t_154[k] = f_4 * kp_59[k]
                   + pb_y[k] * lp_77[k];

        t_155[k] = f_5 * id0_155[k]
                   - f_6 * id1_155[k]
                   + pa_x[k] * kd_155[k];

        t_156[k] = pa_y[k] * kd_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_y, pb_x, pb_y, kp_61, kp_62, \
                         kp_79, kd_122, kd_123, kd_125, lp_79, lp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_4 * kp_79[k]
                   + pb_x[k] * lp_79[k];

        t_158[k] = pa_y[k] * kd_122[k];

        t_159[k] = f_4 * kp_61[k]
                   + pa_y[k] * kd_123[k];

        t_160[k] = f_10 * kp_62[k]
                   + pb_y[k] * lp_80[k];

        t_161[k] = pa_y[k] * kd_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_z, pb_x, pb_y, id0_84, id1_84, kp_83, \
                         kd_120, ls0_27, ls1_27, lp_81, lp_82, lp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_8 * id0_84[k]
                   - f_9 * id1_84[k]
                   + pa_z[k] * kd_120[k];

        t_163[k] = pb_y[k] * lp_81[k];

        t_164[k] = f_4 * kp_83[k]
                   + pb_x[k] * lp_83[k];

        t_165[k] = f_1 * ls0_27[k]
                   - f_2 * ls1_27[k]
                   + pb_y[k] * lp_82[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_x, pb_x, pb_y, id0_167, id1_167, \
                         kp_84, kp_85, kd_167, kd_168, lp_83, lp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_y[k] * lp_83[k];

        t_167[k] = f_5 * id0_167[k]
                   - f_6 * id1_167[k]
                   + pa_x[k] * kd_167[k];

        t_168[k] = f_4 * kp_84[k]
                   + pa_x[k] * kd_168[k];

        t_169[k] = f_10 * kp_85[k]
                   + pb_x[k] * lp_85[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, pa_x, pa_z, pb_z, kd_126, \
                         kd_127, kd_171, kd_173, lp_84, lp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_z[k] * lp_84[k];

        t_171[k] = pa_x[k] * kd_171[k];

        t_172[k] = pb_z[k] * lp_85[k];

        t_173[k] = pa_x[k] * kd_173[k];

        t_174[k] = pa_z[k] * kd_126[k];

        t_175[k] = pa_z[k] * kd_127[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, pa_x, pb_x, kp_89, kp_90, kd_177, \
                         kd_178, kd_179, kd_180, lp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * kp_89[k]
                   + pb_x[k] * lp_89[k];

        t_177[k] = pa_x[k] * kd_177[k];

        t_178[k] = pa_x[k] * kd_178[k];

        t_179[k] = pa_x[k] * kd_179[k];

        t_180[k] = f_4 * kp_90[k]
                   + pa_x[k] * kd_180[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, pa_x, pb_x, kp_91, kp_92, kd_183, \
                         kd_184, kd_185, lp_91, lp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * kp_91[k]
                   + pb_x[k] * lp_91[k];

        t_182[k] = f_10 * kp_92[k]
                   + pb_x[k] * lp_92[k];

        t_183[k] = pa_x[k] * kd_183[k];

        t_184[k] = pa_x[k] * kd_184[k];

        t_185[k] = pa_x[k] * kd_185[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_x, pb_x, kp_93, kp_94, kp_95, \
                         kd_186, kd_189, kd_190, lp_94, lp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_4 * kp_93[k]
                   + pa_x[k] * kd_186[k];

        t_187[k] = f_10 * kp_94[k]
                   + pb_x[k] * lp_94[k];

        t_188[k] = f_10 * kp_95[k]
                   + pb_x[k] * lp_95[k];

        t_189[k] = pa_x[k] * kd_189[k];

        t_190[k] = pa_x[k] * kd_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pa_x, pb_x, kp_96, kp_97, kp_98, \
                         kd_191, kd_192, kd_195, lp_97, lp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = pa_x[k] * kd_191[k];

        t_192[k] = f_4 * kp_96[k]
                   + pa_x[k] * kd_192[k];

        t_193[k] = f_10 * kp_97[k]
                   + pb_x[k] * lp_97[k];

        t_194[k] = f_10 * kp_98[k]
                   + pb_x[k] * lp_98[k];

        t_195[k] = pa_x[k] * kd_195[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_x, pb_x, kp_99, kp_100, kp_101, \
                         kd_196, kd_197, kd_198, lp_100, lp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * kd_196[k];

        t_197[k] = pa_x[k] * kd_197[k];

        t_198[k] = f_4 * kp_99[k]
                   + pa_x[k] * kd_198[k];

        t_199[k] = f_10 * kp_100[k]
                   + pb_x[k] * lp_100[k];

        t_200[k] = f_10 * kp_101[k]
                   + pb_x[k] * lp_101[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pa_x, pa_y, pb_x, kp_103, \
                         kd_162, kd_164, kd_201, kd_202, kd_203, \
                         lp_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_x[k] * kd_201[k];

        t_202[k] = pa_x[k] * kd_202[k];

        t_203[k] = pa_x[k] * kd_203[k];

        t_204[k] = pa_y[k] * kd_162[k];

        t_205[k] = f_10 * kp_103[k]
                   + pb_x[k] * lp_103[k];

        t_206[k] = pa_y[k] * kd_164[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_x, pb_y, kp_105, kd_207, \
                         kd_208, kd_209, kd_210, lp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_x[k] * kd_207[k];

        t_208[k] = pa_x[k] * kd_208[k];

        t_209[k] = pa_x[k] * kd_209[k];

        t_210[k] = f_4 * kp_105[k]
                   + pa_x[k] * kd_210[k];

        t_211[k] = pb_y[k] * lp_105[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pb_x, pb_y, kp_107, kd_213, \
                         kd_215, ls0_36, ls1_36, lp_107, lp_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_10 * kp_107[k]
                   + pb_x[k] * lp_107[k];

        t_213[k] = pa_x[k] * kd_213[k];

        t_214[k] = pb_y[k] * lp_107[k];

        t_215[k] = pa_x[k] * kd_215[k];

        t_216[k] = f_1 * ls0_36[k]
                   - f_2 * ls1_36[k]
                   + pb_x[k] * lp_108[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, t_222, pa_z, pb_x, pb_y, pb_z, \
                         kp_85, kd_168, ls0_36, ls1_36, lp_109, \
                         lp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_x[k] * lp_109[k];

        t_218[k] = pb_x[k] * lp_110[k];

        t_219[k] = f_0 * kp_85[k]
                   + f_1 * ls0_36[k]
                   - f_2 * ls1_36[k]
                   + pb_y[k] * lp_109[k];

        t_220[k] = pb_z[k] * lp_109[k];

        t_221[k] = f_1 * ls0_36[k]
                   - f_2 * ls1_36[k]
                   + pb_z[k] * lp_110[k];

        t_222[k] = pa_z[k] * kd_168[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_y, kp_86, kp_89, \
                         kd_171, kd_173, lp_112, lp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_x[k] * lp_112[k];

        t_224[k] = pb_x[k] * lp_113[k];

        t_225[k] = pa_z[k] * kd_171[k];

        t_226[k] = f_3 * kp_89[k]
                   + pb_y[k] * lp_113[k];

        t_227[k] = f_4 * kp_86[k]
                   + pa_z[k] * kd_173[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_z, pb_x, id0_129, id1_129, kd_177, \
                         ls0_38, ls1_38, lp_114, lp_115, lp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_1 * ls0_38[k]
                   - f_2 * ls1_38[k]
                   + pb_x[k] * lp_114[k];

        t_229[k] = pb_x[k] * lp_115[k];

        t_230[k] = pb_x[k] * lp_116[k];

        t_231[k] = f_5 * id0_129[k]
                   - f_6 * id1_129[k]
                   + pa_z[k] * kd_177[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_y, pb_x, pb_y, id0_143, id1_143, \
                         kp_92, kd_185, ls0_39, ls1_39, lp_116, lp_117, \
                         lp_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_7 * kp_92[k]
                   + pb_y[k] * lp_116[k];

        t_233[k] = f_8 * id0_143[k]
                   - f_9 * id1_143[k]
                   + pa_y[k] * kd_185[k];

        t_234[k] = f_1 * ls0_39[k]
                   - f_2 * ls1_39[k]
                   + pb_x[k] * lp_117[k];

        t_235[k] = pb_x[k] * lp_118[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pa_z, pb_x, pb_y, id0_135, id0_149, \
                         id1_135, id1_149, kp_95, kd_183, kd_191, \
                         lp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_x[k] * lp_119[k];

        t_237[k] = f_11 * id0_135[k]
                   - f_12 * id1_135[k]
                   + pa_z[k] * kd_183[k];

        t_238[k] = f_13 * kp_95[k]
                   + pb_y[k] * lp_119[k];

        t_239[k] = f_14 * id0_149[k]
                   - f_15 * id1_149[k]
                   + pa_y[k] * kd_191[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pa_z, pb_x, id0_141, id1_141, kd_189, \
                         ls0_40, ls1_40, lp_120, lp_121, lp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_1 * ls0_40[k]
                   - f_2 * ls1_40[k]
                   + pb_x[k] * lp_120[k];

        t_241[k] = pb_x[k] * lp_121[k];

        t_242[k] = pb_x[k] * lp_122[k];

        t_243[k] = f_16 * id0_141[k]
                   - f_17 * id1_141[k]
                   + pa_z[k] * kd_189[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_y, pb_x, pb_y, id0_155, id1_155, \
                         kp_98, kd_197, ls0_41, ls1_41, lp_122, lp_123, \
                         lp_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_18 * kp_98[k]
                   + pb_y[k] * lp_122[k];

        t_245[k] = f_16 * id0_155[k]
                   - f_17 * id1_155[k]
                   + pa_y[k] * kd_197[k];

        t_246[k] = f_1 * ls0_41[k]
                   - f_2 * ls1_41[k]
                   + pb_x[k] * lp_123[k];

        t_247[k] = pb_x[k] * lp_124[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_y, pa_z, pb_x, pb_y, id0_147, id0_161, \
                         id1_147, id1_161, kp_101, kd_195, kd_203, \
                         lp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pb_x[k] * lp_125[k];

        t_249[k] = f_14 * id0_147[k]
                   - f_15 * id1_147[k]
                   + pa_z[k] * kd_195[k];

        t_250[k] = f_19 * kp_101[k]
                   + pb_y[k] * lp_125[k];

        t_251[k] = f_11 * id0_161[k]
                   - f_12 * id1_161[k]
                   + pa_y[k] * kd_203[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_z, pb_x, id0_153, id1_153, kd_201, \
                         ls0_42, ls1_42, lp_126, lp_127, lp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * ls0_42[k]
                   - f_2 * ls1_42[k]
                   + pb_x[k] * lp_126[k];

        t_253[k] = pb_x[k] * lp_127[k];

        t_254[k] = pb_x[k] * lp_128[k];

        t_255[k] = f_8 * id0_153[k]
                   - f_9 * id1_153[k]
                   + pa_z[k] * kd_201[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pa_y, pb_x, pb_y, id0_167, \
                         id1_167, kp_104, kd_209, kd_210, lp_128, lp_130, \
                         lp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_4 * kp_104[k]
                   + pb_y[k] * lp_128[k];

        t_257[k] = f_5 * id0_167[k]
                   - f_6 * id1_167[k]
                   + pa_y[k] * kd_209[k];

        t_258[k] = pa_y[k] * kd_210[k];

        t_259[k] = pb_x[k] * lp_130[k];

        t_260[k] = pb_x[k] * lp_131[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pb_x, pb_y, kp_106, kp_107, kd_213, \
                         kd_215, ls0_44, ls1_44, lp_131, lp_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_4 * kp_106[k]
                   + pa_y[k] * kd_213[k];

        t_262[k] = f_10 * kp_107[k]
                   + pb_y[k] * lp_131[k];

        t_263[k] = pa_y[k] * kd_215[k];

        t_264[k] = f_1 * ls0_44[k]
                   - f_2 * ls1_44[k]
                   + pb_x[k] * lp_132[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pb_x, pb_y, pb_z, kp_107, ls0_44, \
                         ls1_44, lp_133, lp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = pb_x[k] * lp_133[k];

        t_266[k] = pb_x[k] * lp_134[k];

        t_267[k] = f_1 * ls0_44[k]
                   - f_2 * ls1_44[k]
                   + pb_y[k] * lp_133[k];

        t_268[k] = pb_y[k] * lp_134[k];

        t_269[k] = f_0 * kp_107[k]
                   + f_1 * ls0_44[k]
                   - f_2 * ls1_44[k]
                   + pb_z[k] * lp_134[k];
    }
}

}  // namespace simdt2ceri
