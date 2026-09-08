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


#include "SimdKineticEnergyVrrRecFI.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_fi_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pi_s, const size_t pi,
                                 const size_t dh, const size_t di, const size_t fg_s,
                                 const size_t fi_s, const size_t fg, const size_t fh,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 5.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 3.0 * alpha / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = beta / p;
    const auto f_11 = 4.0 * alpha / p;
    const auto f_12 = 3.0 / p;

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

    const auto *pi_s_49 = buffer.data(pi_s + 49);
    const auto *pi_s_83 = buffer.data(pi_s + 83);

    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_83 = buffer.data(pi + 83);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_167 = buffer.data(di + 167);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_90 = buffer.data(fg_s + 90);
    const auto *fg_s_91 = buffer.data(fg_s + 91);
    const auto *fg_s_93 = buffer.data(fg_s + 93);
    const auto *fg_s_95 = buffer.data(fg_s + 95);
    const auto *fg_s_96 = buffer.data(fg_s + 96);
    const auto *fg_s_98 = buffer.data(fg_s + 98);
    const auto *fg_s_99 = buffer.data(fg_s + 99);
    const auto *fg_s_100 = buffer.data(fg_s + 100);
    const auto *fg_s_101 = buffer.data(fg_s + 101);
    const auto *fg_s_102 = buffer.data(fg_s + 102);
    const auto *fg_s_103 = buffer.data(fg_s + 103);
    const auto *fg_s_104 = buffer.data(fg_s + 104);
    const auto *fg_s_107 = buffer.data(fg_s + 107);
    const auto *fg_s_110 = buffer.data(fg_s + 110);
    const auto *fg_s_114 = buffer.data(fg_s + 114);
    const auto *fg_s_119 = buffer.data(fg_s + 119);
    const auto *fg_s_135 = buffer.data(fg_s + 135);
    const auto *fg_s_137 = buffer.data(fg_s + 137);
    const auto *fg_s_138 = buffer.data(fg_s + 138);
    const auto *fg_s_140 = buffer.data(fg_s + 140);
    const auto *fg_s_141 = buffer.data(fg_s + 141);
    const auto *fg_s_142 = buffer.data(fg_s + 142);
    const auto *fg_s_144 = buffer.data(fg_s + 144);
    const auto *fg_s_145 = buffer.data(fg_s + 145);
    const auto *fg_s_146 = buffer.data(fg_s + 146);
    const auto *fg_s_147 = buffer.data(fg_s + 147);
    const auto *fg_s_148 = buffer.data(fg_s + 148);
    const auto *fg_s_149 = buffer.data(fg_s + 149);

    const auto *fi_s_0 = buffer.data(fi_s + 0);
    const auto *fi_s_1 = buffer.data(fi_s + 1);
    const auto *fi_s_2 = buffer.data(fi_s + 2);
    const auto *fi_s_3 = buffer.data(fi_s + 3);
    const auto *fi_s_4 = buffer.data(fi_s + 4);
    const auto *fi_s_5 = buffer.data(fi_s + 5);
    const auto *fi_s_6 = buffer.data(fi_s + 6);
    const auto *fi_s_7 = buffer.data(fi_s + 7);
    const auto *fi_s_8 = buffer.data(fi_s + 8);
    const auto *fi_s_9 = buffer.data(fi_s + 9);
    const auto *fi_s_10 = buffer.data(fi_s + 10);
    const auto *fi_s_11 = buffer.data(fi_s + 11);
    const auto *fi_s_12 = buffer.data(fi_s + 12);
    const auto *fi_s_13 = buffer.data(fi_s + 13);
    const auto *fi_s_14 = buffer.data(fi_s + 14);
    const auto *fi_s_15 = buffer.data(fi_s + 15);
    const auto *fi_s_16 = buffer.data(fi_s + 16);
    const auto *fi_s_17 = buffer.data(fi_s + 17);
    const auto *fi_s_18 = buffer.data(fi_s + 18);
    const auto *fi_s_19 = buffer.data(fi_s + 19);
    const auto *fi_s_20 = buffer.data(fi_s + 20);
    const auto *fi_s_21 = buffer.data(fi_s + 21);
    const auto *fi_s_22 = buffer.data(fi_s + 22);
    const auto *fi_s_23 = buffer.data(fi_s + 23);
    const auto *fi_s_24 = buffer.data(fi_s + 24);
    const auto *fi_s_25 = buffer.data(fi_s + 25);
    const auto *fi_s_26 = buffer.data(fi_s + 26);
    const auto *fi_s_27 = buffer.data(fi_s + 27);
    const auto *fi_s_28 = buffer.data(fi_s + 28);
    const auto *fi_s_29 = buffer.data(fi_s + 29);
    const auto *fi_s_30 = buffer.data(fi_s + 30);
    const auto *fi_s_31 = buffer.data(fi_s + 31);
    const auto *fi_s_32 = buffer.data(fi_s + 32);
    const auto *fi_s_33 = buffer.data(fi_s + 33);
    const auto *fi_s_34 = buffer.data(fi_s + 34);
    const auto *fi_s_35 = buffer.data(fi_s + 35);
    const auto *fi_s_36 = buffer.data(fi_s + 36);
    const auto *fi_s_37 = buffer.data(fi_s + 37);
    const auto *fi_s_38 = buffer.data(fi_s + 38);
    const auto *fi_s_39 = buffer.data(fi_s + 39);
    const auto *fi_s_40 = buffer.data(fi_s + 40);
    const auto *fi_s_41 = buffer.data(fi_s + 41);
    const auto *fi_s_42 = buffer.data(fi_s + 42);
    const auto *fi_s_43 = buffer.data(fi_s + 43);
    const auto *fi_s_44 = buffer.data(fi_s + 44);
    const auto *fi_s_45 = buffer.data(fi_s + 45);
    const auto *fi_s_46 = buffer.data(fi_s + 46);
    const auto *fi_s_47 = buffer.data(fi_s + 47);
    const auto *fi_s_48 = buffer.data(fi_s + 48);
    const auto *fi_s_49 = buffer.data(fi_s + 49);
    const auto *fi_s_50 = buffer.data(fi_s + 50);
    const auto *fi_s_51 = buffer.data(fi_s + 51);
    const auto *fi_s_52 = buffer.data(fi_s + 52);
    const auto *fi_s_53 = buffer.data(fi_s + 53);
    const auto *fi_s_54 = buffer.data(fi_s + 54);
    const auto *fi_s_55 = buffer.data(fi_s + 55);
    const auto *fi_s_56 = buffer.data(fi_s + 56);
    const auto *fi_s_57 = buffer.data(fi_s + 57);
    const auto *fi_s_58 = buffer.data(fi_s + 58);
    const auto *fi_s_59 = buffer.data(fi_s + 59);
    const auto *fi_s_60 = buffer.data(fi_s + 60);
    const auto *fi_s_61 = buffer.data(fi_s + 61);
    const auto *fi_s_62 = buffer.data(fi_s + 62);
    const auto *fi_s_63 = buffer.data(fi_s + 63);
    const auto *fi_s_64 = buffer.data(fi_s + 64);
    const auto *fi_s_65 = buffer.data(fi_s + 65);
    const auto *fi_s_66 = buffer.data(fi_s + 66);
    const auto *fi_s_67 = buffer.data(fi_s + 67);
    const auto *fi_s_68 = buffer.data(fi_s + 68);
    const auto *fi_s_69 = buffer.data(fi_s + 69);
    const auto *fi_s_70 = buffer.data(fi_s + 70);
    const auto *fi_s_71 = buffer.data(fi_s + 71);
    const auto *fi_s_72 = buffer.data(fi_s + 72);
    const auto *fi_s_73 = buffer.data(fi_s + 73);
    const auto *fi_s_74 = buffer.data(fi_s + 74);
    const auto *fi_s_75 = buffer.data(fi_s + 75);
    const auto *fi_s_76 = buffer.data(fi_s + 76);
    const auto *fi_s_77 = buffer.data(fi_s + 77);
    const auto *fi_s_78 = buffer.data(fi_s + 78);
    const auto *fi_s_79 = buffer.data(fi_s + 79);
    const auto *fi_s_80 = buffer.data(fi_s + 80);
    const auto *fi_s_81 = buffer.data(fi_s + 81);
    const auto *fi_s_82 = buffer.data(fi_s + 82);
    const auto *fi_s_83 = buffer.data(fi_s + 83);
    const auto *fi_s_84 = buffer.data(fi_s + 84);
    const auto *fi_s_85 = buffer.data(fi_s + 85);
    const auto *fi_s_86 = buffer.data(fi_s + 86);
    const auto *fi_s_87 = buffer.data(fi_s + 87);
    const auto *fi_s_88 = buffer.data(fi_s + 88);
    const auto *fi_s_89 = buffer.data(fi_s + 89);
    const auto *fi_s_90 = buffer.data(fi_s + 90);
    const auto *fi_s_91 = buffer.data(fi_s + 91);
    const auto *fi_s_92 = buffer.data(fi_s + 92);
    const auto *fi_s_93 = buffer.data(fi_s + 93);
    const auto *fi_s_94 = buffer.data(fi_s + 94);
    const auto *fi_s_95 = buffer.data(fi_s + 95);
    const auto *fi_s_96 = buffer.data(fi_s + 96);
    const auto *fi_s_97 = buffer.data(fi_s + 97);
    const auto *fi_s_98 = buffer.data(fi_s + 98);
    const auto *fi_s_99 = buffer.data(fi_s + 99);
    const auto *fi_s_100 = buffer.data(fi_s + 100);
    const auto *fi_s_101 = buffer.data(fi_s + 101);
    const auto *fi_s_102 = buffer.data(fi_s + 102);
    const auto *fi_s_103 = buffer.data(fi_s + 103);
    const auto *fi_s_104 = buffer.data(fi_s + 104);
    const auto *fi_s_105 = buffer.data(fi_s + 105);
    const auto *fi_s_106 = buffer.data(fi_s + 106);
    const auto *fi_s_107 = buffer.data(fi_s + 107);
    const auto *fi_s_108 = buffer.data(fi_s + 108);
    const auto *fi_s_109 = buffer.data(fi_s + 109);
    const auto *fi_s_110 = buffer.data(fi_s + 110);
    const auto *fi_s_111 = buffer.data(fi_s + 111);
    const auto *fi_s_112 = buffer.data(fi_s + 112);
    const auto *fi_s_113 = buffer.data(fi_s + 113);
    const auto *fi_s_114 = buffer.data(fi_s + 114);
    const auto *fi_s_115 = buffer.data(fi_s + 115);
    const auto *fi_s_116 = buffer.data(fi_s + 116);
    const auto *fi_s_117 = buffer.data(fi_s + 117);
    const auto *fi_s_118 = buffer.data(fi_s + 118);
    const auto *fi_s_119 = buffer.data(fi_s + 119);
    const auto *fi_s_120 = buffer.data(fi_s + 120);
    const auto *fi_s_121 = buffer.data(fi_s + 121);
    const auto *fi_s_122 = buffer.data(fi_s + 122);
    const auto *fi_s_123 = buffer.data(fi_s + 123);
    const auto *fi_s_124 = buffer.data(fi_s + 124);
    const auto *fi_s_125 = buffer.data(fi_s + 125);
    const auto *fi_s_126 = buffer.data(fi_s + 126);
    const auto *fi_s_127 = buffer.data(fi_s + 127);
    const auto *fi_s_128 = buffer.data(fi_s + 128);
    const auto *fi_s_129 = buffer.data(fi_s + 129);
    const auto *fi_s_130 = buffer.data(fi_s + 130);
    const auto *fi_s_131 = buffer.data(fi_s + 131);
    const auto *fi_s_132 = buffer.data(fi_s + 132);
    const auto *fi_s_133 = buffer.data(fi_s + 133);
    const auto *fi_s_134 = buffer.data(fi_s + 134);
    const auto *fi_s_135 = buffer.data(fi_s + 135);
    const auto *fi_s_136 = buffer.data(fi_s + 136);
    const auto *fi_s_137 = buffer.data(fi_s + 137);
    const auto *fi_s_138 = buffer.data(fi_s + 138);
    const auto *fi_s_139 = buffer.data(fi_s + 139);
    const auto *fi_s_140 = buffer.data(fi_s + 140);
    const auto *fi_s_141 = buffer.data(fi_s + 141);
    const auto *fi_s_142 = buffer.data(fi_s + 142);
    const auto *fi_s_143 = buffer.data(fi_s + 143);
    const auto *fi_s_144 = buffer.data(fi_s + 144);
    const auto *fi_s_145 = buffer.data(fi_s + 145);
    const auto *fi_s_146 = buffer.data(fi_s + 146);
    const auto *fi_s_147 = buffer.data(fi_s + 147);
    const auto *fi_s_148 = buffer.data(fi_s + 148);
    const auto *fi_s_149 = buffer.data(fi_s + 149);
    const auto *fi_s_150 = buffer.data(fi_s + 150);
    const auto *fi_s_151 = buffer.data(fi_s + 151);
    const auto *fi_s_152 = buffer.data(fi_s + 152);
    const auto *fi_s_153 = buffer.data(fi_s + 153);
    const auto *fi_s_154 = buffer.data(fi_s + 154);
    const auto *fi_s_155 = buffer.data(fi_s + 155);
    const auto *fi_s_156 = buffer.data(fi_s + 156);
    const auto *fi_s_157 = buffer.data(fi_s + 157);
    const auto *fi_s_158 = buffer.data(fi_s + 158);
    const auto *fi_s_159 = buffer.data(fi_s + 159);
    const auto *fi_s_160 = buffer.data(fi_s + 160);
    const auto *fi_s_161 = buffer.data(fi_s + 161);
    const auto *fi_s_162 = buffer.data(fi_s + 162);
    const auto *fi_s_163 = buffer.data(fi_s + 163);
    const auto *fi_s_164 = buffer.data(fi_s + 164);
    const auto *fi_s_165 = buffer.data(fi_s + 165);
    const auto *fi_s_166 = buffer.data(fi_s + 166);
    const auto *fi_s_167 = buffer.data(fi_s + 167);
    const auto *fi_s_168 = buffer.data(fi_s + 168);
    const auto *fi_s_169 = buffer.data(fi_s + 169);
    const auto *fi_s_170 = buffer.data(fi_s + 170);
    const auto *fi_s_171 = buffer.data(fi_s + 171);
    const auto *fi_s_172 = buffer.data(fi_s + 172);
    const auto *fi_s_173 = buffer.data(fi_s + 173);
    const auto *fi_s_174 = buffer.data(fi_s + 174);
    const auto *fi_s_175 = buffer.data(fi_s + 175);
    const auto *fi_s_176 = buffer.data(fi_s + 176);
    const auto *fi_s_177 = buffer.data(fi_s + 177);
    const auto *fi_s_178 = buffer.data(fi_s + 178);
    const auto *fi_s_179 = buffer.data(fi_s + 179);
    const auto *fi_s_180 = buffer.data(fi_s + 180);
    const auto *fi_s_181 = buffer.data(fi_s + 181);
    const auto *fi_s_182 = buffer.data(fi_s + 182);
    const auto *fi_s_183 = buffer.data(fi_s + 183);
    const auto *fi_s_184 = buffer.data(fi_s + 184);
    const auto *fi_s_185 = buffer.data(fi_s + 185);
    const auto *fi_s_186 = buffer.data(fi_s + 186);
    const auto *fi_s_187 = buffer.data(fi_s + 187);
    const auto *fi_s_188 = buffer.data(fi_s + 188);
    const auto *fi_s_189 = buffer.data(fi_s + 189);
    const auto *fi_s_190 = buffer.data(fi_s + 190);
    const auto *fi_s_191 = buffer.data(fi_s + 191);
    const auto *fi_s_192 = buffer.data(fi_s + 192);
    const auto *fi_s_193 = buffer.data(fi_s + 193);
    const auto *fi_s_194 = buffer.data(fi_s + 194);
    const auto *fi_s_195 = buffer.data(fi_s + 195);
    const auto *fi_s_196 = buffer.data(fi_s + 196);
    const auto *fi_s_197 = buffer.data(fi_s + 197);
    const auto *fi_s_198 = buffer.data(fi_s + 198);
    const auto *fi_s_199 = buffer.data(fi_s + 199);
    const auto *fi_s_200 = buffer.data(fi_s + 200);
    const auto *fi_s_201 = buffer.data(fi_s + 201);
    const auto *fi_s_202 = buffer.data(fi_s + 202);
    const auto *fi_s_203 = buffer.data(fi_s + 203);
    const auto *fi_s_204 = buffer.data(fi_s + 204);
    const auto *fi_s_205 = buffer.data(fi_s + 205);
    const auto *fi_s_206 = buffer.data(fi_s + 206);
    const auto *fi_s_207 = buffer.data(fi_s + 207);
    const auto *fi_s_208 = buffer.data(fi_s + 208);
    const auto *fi_s_209 = buffer.data(fi_s + 209);
    const auto *fi_s_210 = buffer.data(fi_s + 210);
    const auto *fi_s_211 = buffer.data(fi_s + 211);
    const auto *fi_s_212 = buffer.data(fi_s + 212);
    const auto *fi_s_213 = buffer.data(fi_s + 213);
    const auto *fi_s_214 = buffer.data(fi_s + 214);
    const auto *fi_s_215 = buffer.data(fi_s + 215);
    const auto *fi_s_216 = buffer.data(fi_s + 216);
    const auto *fi_s_217 = buffer.data(fi_s + 217);
    const auto *fi_s_218 = buffer.data(fi_s + 218);
    const auto *fi_s_219 = buffer.data(fi_s + 219);
    const auto *fi_s_220 = buffer.data(fi_s + 220);
    const auto *fi_s_221 = buffer.data(fi_s + 221);
    const auto *fi_s_222 = buffer.data(fi_s + 222);
    const auto *fi_s_223 = buffer.data(fi_s + 223);
    const auto *fi_s_224 = buffer.data(fi_s + 224);
    const auto *fi_s_225 = buffer.data(fi_s + 225);
    const auto *fi_s_226 = buffer.data(fi_s + 226);
    const auto *fi_s_227 = buffer.data(fi_s + 227);
    const auto *fi_s_228 = buffer.data(fi_s + 228);
    const auto *fi_s_229 = buffer.data(fi_s + 229);
    const auto *fi_s_230 = buffer.data(fi_s + 230);
    const auto *fi_s_231 = buffer.data(fi_s + 231);
    const auto *fi_s_232 = buffer.data(fi_s + 232);
    const auto *fi_s_233 = buffer.data(fi_s + 233);
    const auto *fi_s_234 = buffer.data(fi_s + 234);
    const auto *fi_s_235 = buffer.data(fi_s + 235);
    const auto *fi_s_236 = buffer.data(fi_s + 236);
    const auto *fi_s_237 = buffer.data(fi_s + 237);
    const auto *fi_s_238 = buffer.data(fi_s + 238);
    const auto *fi_s_239 = buffer.data(fi_s + 239);
    const auto *fi_s_240 = buffer.data(fi_s + 240);
    const auto *fi_s_241 = buffer.data(fi_s + 241);
    const auto *fi_s_242 = buffer.data(fi_s + 242);
    const auto *fi_s_243 = buffer.data(fi_s + 243);
    const auto *fi_s_244 = buffer.data(fi_s + 244);
    const auto *fi_s_245 = buffer.data(fi_s + 245);
    const auto *fi_s_246 = buffer.data(fi_s + 246);
    const auto *fi_s_247 = buffer.data(fi_s + 247);
    const auto *fi_s_248 = buffer.data(fi_s + 248);
    const auto *fi_s_249 = buffer.data(fi_s + 249);
    const auto *fi_s_250 = buffer.data(fi_s + 250);
    const auto *fi_s_251 = buffer.data(fi_s + 251);
    const auto *fi_s_252 = buffer.data(fi_s + 252);
    const auto *fi_s_253 = buffer.data(fi_s + 253);
    const auto *fi_s_254 = buffer.data(fi_s + 254);
    const auto *fi_s_255 = buffer.data(fi_s + 255);
    const auto *fi_s_256 = buffer.data(fi_s + 256);
    const auto *fi_s_257 = buffer.data(fi_s + 257);
    const auto *fi_s_258 = buffer.data(fi_s + 258);
    const auto *fi_s_259 = buffer.data(fi_s + 259);
    const auto *fi_s_260 = buffer.data(fi_s + 260);
    const auto *fi_s_261 = buffer.data(fi_s + 261);
    const auto *fi_s_262 = buffer.data(fi_s + 262);
    const auto *fi_s_263 = buffer.data(fi_s + 263);
    const auto *fi_s_264 = buffer.data(fi_s + 264);
    const auto *fi_s_265 = buffer.data(fi_s + 265);
    const auto *fi_s_266 = buffer.data(fi_s + 266);
    const auto *fi_s_267 = buffer.data(fi_s + 267);
    const auto *fi_s_268 = buffer.data(fi_s + 268);
    const auto *fi_s_269 = buffer.data(fi_s + 269);
    const auto *fi_s_270 = buffer.data(fi_s + 270);
    const auto *fi_s_271 = buffer.data(fi_s + 271);
    const auto *fi_s_272 = buffer.data(fi_s + 272);
    const auto *fi_s_273 = buffer.data(fi_s + 273);
    const auto *fi_s_274 = buffer.data(fi_s + 274);
    const auto *fi_s_275 = buffer.data(fi_s + 275);
    const auto *fi_s_276 = buffer.data(fi_s + 276);
    const auto *fi_s_277 = buffer.data(fi_s + 277);
    const auto *fi_s_278 = buffer.data(fi_s + 278);
    const auto *fi_s_279 = buffer.data(fi_s + 279);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dh_0, fg_s_0, fi_s_0, fi_s_1, \
                         fi_s_2, fi_s_3, fg_0, fh_0, fh_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_0[k]
                 - f_1 * fg_s_0[k]
                 + f_2 * fi_s_0[k]
                 + f_3 * fg_0[k]
                 + pb_x[k] * fh_0[k];

        t_1[k] = f_2 * fi_s_1[k]
                 + pb_y[k] * fh_0[k];

        t_2[k] = f_2 * fi_s_2[k]
                 + pb_z[k] * fh_0[k];

        t_3[k] = -f_4 * fg_s_0[k]
                 + f_2 * fi_s_3[k]
                 + f_5 * fg_0[k]
                 + pb_y[k] * fh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, fg_s_0, fg_s_1, fi_s_4, fi_s_5, \
                         fi_s_6, fi_s_7, fg_0, fg_1, fh_2, fh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fi_s_4[k]
                 + pb_y[k] * fh_2[k];

        t_5[k] = -f_4 * fg_s_0[k]
                 + f_2 * fi_s_5[k]
                 + f_5 * fg_0[k]
                 + pb_z[k] * fh_2[k];

        t_6[k] = -f_6 * fg_s_1[k]
                 + f_2 * fi_s_6[k]
                 + f_7 * fg_1[k]
                 + pb_y[k] * fh_3[k];

        t_7[k] = f_2 * fi_s_7[k]
                 + pb_z[k] * fh_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, fg_s_2, fg_s_3, fi_s_8, fi_s_9, \
                         fi_s_10, fi_s_11, fg_2, fg_3, fh_5, fh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * fi_s_8[k]
                 + pb_y[k] * fh_5[k];

        t_9[k] = -f_6 * fg_s_2[k]
                 + f_2 * fi_s_9[k]
                 + f_7 * fg_2[k]
                 + pb_z[k] * fh_5[k];

        t_10[k] = -f_8 * fg_s_3[k]
                  + f_2 * fi_s_10[k]
                  + f_0 * fg_3[k]
                  + pb_y[k] * fh_6[k];

        t_11[k] = f_2 * fi_s_11[k]
                  + pb_z[k] * fh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, fg_s_5, fi_s_12, fi_s_13, fi_s_14, \
                         fg_5, fh_8, fh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_4 * fg_s_5[k]
                  + f_2 * fi_s_12[k]
                  + f_5 * fg_5[k]
                  + pb_y[k] * fh_8[k];

        t_13[k] = f_2 * fi_s_13[k]
                  + pb_y[k] * fh_9[k];

        t_14[k] = -f_8 * fg_s_5[k]
                  + f_2 * fi_s_14[k]
                  + f_0 * fg_5[k]
                  + pb_z[k] * fh_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_z, dh_15, dh_17, fi_s_15, fi_s_16, \
                         fi_s_17, fh_10, fh_15, fh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * dh_15[k]
                  + f_2 * fi_s_15[k]
                  + pb_x[k] * fh_15[k];

        t_16[k] = f_2 * fi_s_16[k]
                  + pb_z[k] * fh_10[k];

        t_17[k] = f_0 * dh_17[k]
                  + f_2 * fi_s_17[k]
                  + pb_x[k] * fh_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_y, dh_18, dh_20, fi_s_18, fi_s_19, \
                         fi_s_20, fh_14, fh_18, fh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * dh_18[k]
                  + f_2 * fi_s_18[k]
                  + pb_x[k] * fh_18[k];

        t_19[k] = f_2 * fi_s_19[k]
                  + pb_y[k] * fh_14[k];

        t_20[k] = f_0 * dh_20[k]
                  + f_2 * fi_s_20[k]
                  + pb_x[k] * fh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, fg_s_10, fg_s_12, fi_s_21, fi_s_22, \
                         fi_s_23, fg_10, fg_12, fh_15, fh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * fg_s_10[k]
                  + f_2 * fi_s_21[k]
                  + f_3 * fg_10[k]
                  + pb_y[k] * fh_15[k];

        t_22[k] = f_2 * fi_s_22[k]
                  + pb_z[k] * fh_15[k];

        t_23[k] = -f_8 * fg_s_12[k]
                  + f_2 * fi_s_23[k]
                  + f_0 * fg_12[k]
                  + pb_y[k] * fh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, fg_s_13, fg_s_14, fi_s_24, fi_s_25, fi_s_26, \
                         fg_13, fg_14, fh_18, fh_19, fh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_6 * fg_s_13[k]
                  + f_2 * fi_s_24[k]
                  + f_7 * fg_13[k]
                  + pb_y[k] * fh_18[k];

        t_25[k] = -f_4 * fg_s_14[k]
                  + f_2 * fi_s_25[k]
                  + f_5 * fg_14[k]
                  + pb_y[k] * fh_19[k];

        t_26[k] = f_2 * fi_s_26[k]
                  + pb_y[k] * fh_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, dh_0, di_0, fg_s_14, fi_s_27, \
                         fi_s_28, fi_s_29, fg_14, fh_20, fh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * fg_s_14[k]
                  + f_2 * fi_s_27[k]
                  + f_3 * fg_14[k]
                  + pb_z[k] * fh_20[k];

        t_28[k] = pa_y[k] * di_0[k]
                  + f_2 * fi_s_28[k];

        t_29[k] = f_5 * dh_0[k]
                  + f_2 * fi_s_29[k]
                  + pb_y[k] * fh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_z, dh_1, di_3, di_5, fi_s_30, \
                         fi_s_31, fi_s_32, fi_s_33, fh_21, fh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * fi_s_30[k]
                  + pb_z[k] * fh_21[k];

        t_31[k] = f_7 * dh_1[k]
                  + pa_y[k] * di_3[k]
                  + f_2 * fi_s_31[k];

        t_32[k] = f_2 * fi_s_32[k]
                  + pb_z[k] * fh_22[k];

        t_33[k] = pa_y[k] * di_5[k]
                  + f_2 * fi_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_y, pb_z, dh_3, dh_5, di_6, fi_s_34, \
                         fi_s_35, fi_s_36, fh_24, fh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dh_3[k]
                  + pa_y[k] * di_6[k]
                  + f_2 * fi_s_34[k];

        t_35[k] = f_2 * fi_s_35[k]
                  + pb_z[k] * fh_24[k];

        t_36[k] = f_5 * dh_5[k]
                  + f_2 * fi_s_36[k]
                  + pb_y[k] * fh_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pb_z, dh_6, dh_8, di_9, di_10, di_12, \
                         fi_s_37, fi_s_38, fi_s_39, fi_s_40, fh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * di_9[k]
                  + f_2 * fi_s_37[k];

        t_38[k] = f_9 * dh_6[k]
                  + pa_y[k] * di_10[k]
                  + f_2 * fi_s_38[k];

        t_39[k] = f_2 * fi_s_39[k]
                  + pb_z[k] * fh_27[k];

        t_40[k] = f_7 * dh_8[k]
                  + pa_y[k] * di_12[k]
                  + f_2 * fi_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pb_x, pb_y, dh_9, dh_36, di_14, fi_s_41, \
                         fi_s_42, fi_s_43, fh_30, fh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * dh_9[k]
                  + f_2 * fi_s_41[k]
                  + pb_y[k] * fh_30[k];

        t_42[k] = pa_y[k] * di_14[k]
                  + f_2 * fi_s_42[k];

        t_43[k] = f_7 * dh_36[k]
                  + f_2 * fi_s_43[k]
                  + pb_x[k] * fh_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_z, dh_38, dh_39, fi_s_44, fi_s_45, \
                         fi_s_46, fh_31, fh_38, fh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * fi_s_44[k]
                  + pb_z[k] * fh_31[k];

        t_45[k] = f_7 * dh_38[k]
                  + f_2 * fi_s_45[k]
                  + pb_x[k] * fh_38[k];

        t_46[k] = f_7 * dh_39[k]
                  + f_2 * fi_s_46[k]
                  + pb_x[k] * fh_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_y, pb_x, pi_s_49, pi_49, dh_40, di_20, \
                         di_49, fi_s_47, fi_s_48, fi_s_49, fh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_7 * dh_40[k]
                  + f_2 * fi_s_47[k]
                  + pb_x[k] * fh_40[k];

        t_48[k] = pa_y[k] * di_20[k]
                  + f_2 * fi_s_48[k];

        t_49[k] = -f_10 * pi_s_49[k]
                  + f_5 * pi_49[k]
                  + pa_x[k] * di_49[k]
                  + f_2 * fi_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_z, fg_s_25, fg_s_26, fi_s_50, fi_s_51, fi_s_52, \
                         fg_25, fg_26, fh_36, fh_37, fh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * fi_s_50[k]
                  + pb_z[k] * fh_36[k];

        t_51[k] = -f_4 * fg_s_25[k]
                  + f_2 * fi_s_51[k]
                  + f_5 * fg_25[k]
                  + pb_z[k] * fh_37[k];

        t_52[k] = -f_6 * fg_s_26[k]
                  + f_2 * fi_s_52[k]
                  + f_7 * fg_26[k]
                  + pb_z[k] * fh_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pb_y, pb_z, dh_20, di_27, fg_s_27, fi_s_53, \
                         fi_s_54, fi_s_55, fg_27, fh_39, fh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_8 * fg_s_27[k]
                  + f_2 * fi_s_53[k]
                  + f_0 * fg_27[k]
                  + pb_z[k] * fh_39[k];

        t_54[k] = f_5 * dh_20[k]
                  + f_2 * fi_s_54[k]
                  + pb_y[k] * fh_41[k];

        t_55[k] = pa_y[k] * di_27[k]
                  + f_2 * fi_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_y, pb_z, dh_0, di_0, di_3, fi_s_56, \
                         fi_s_57, fi_s_58, fi_s_59, fh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_z[k] * di_0[k]
                  + f_2 * fi_s_56[k];

        t_57[k] = f_2 * fi_s_57[k]
                  + pb_y[k] * fh_42[k];

        t_58[k] = f_5 * dh_0[k]
                  + f_2 * fi_s_58[k]
                  + pb_z[k] * fh_42[k];

        t_59[k] = pa_z[k] * di_3[k]
                  + f_2 * fi_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_y, dh_2, dh_3, di_5, di_6, di_7, \
                         fi_s_60, fi_s_61, fi_s_62, fi_s_63, fh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * fi_s_60[k]
                  + pb_y[k] * fh_44[k];

        t_61[k] = f_7 * dh_2[k]
                  + pa_z[k] * di_5[k]
                  + f_2 * fi_s_61[k];

        t_62[k] = pa_z[k] * di_6[k]
                  + f_2 * fi_s_62[k];

        t_63[k] = f_5 * dh_3[k]
                  + pa_z[k] * di_7[k]
                  + f_2 * fi_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_y, dh_5, dh_6, di_9, di_10, di_11, \
                         fi_s_64, fi_s_65, fi_s_66, fi_s_67, fh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * fi_s_64[k]
                  + pb_y[k] * fh_47[k];

        t_65[k] = f_0 * dh_5[k]
                  + pa_z[k] * di_9[k]
                  + f_2 * fi_s_65[k];

        t_66[k] = pa_z[k] * di_10[k]
                  + f_2 * fi_s_66[k];

        t_67[k] = f_5 * dh_6[k]
                  + pa_z[k] * di_11[k]
                  + f_2 * fi_s_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, dh_7, dh_9, di_12, di_14, di_15, \
                         fi_s_68, fi_s_69, fi_s_70, fi_s_71, fh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_7 * dh_7[k]
                  + pa_z[k] * di_12[k]
                  + f_2 * fi_s_68[k];

        t_69[k] = f_2 * fi_s_69[k]
                  + pb_y[k] * fh_51[k];

        t_70[k] = f_9 * dh_9[k]
                  + pa_z[k] * di_14[k]
                  + f_2 * fi_s_70[k];

        t_71[k] = pa_z[k] * di_15[k]
                  + f_2 * fi_s_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, dh_58, dh_59, dh_60, fi_s_72, fi_s_73, \
                         fi_s_74, fh_58, fh_59, fh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * dh_58[k]
                  + f_2 * fi_s_72[k]
                  + pb_x[k] * fh_58[k];

        t_73[k] = f_7 * dh_59[k]
                  + f_2 * fi_s_73[k]
                  + pb_x[k] * fh_59[k];

        t_74[k] = f_7 * dh_60[k]
                  + f_2 * fi_s_74[k]
                  + pb_x[k] * fh_60[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_x, pb_y, dh_62, di_21, fi_s_75, fi_s_76, \
                         fi_s_77, fh_56, fh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * fi_s_75[k]
                  + pb_y[k] * fh_56[k];

        t_76[k] = f_7 * dh_62[k]
                  + f_2 * fi_s_76[k]
                  + pb_x[k] * fh_62[k];

        t_77[k] = pa_z[k] * di_21[k]
                  + f_2 * fi_s_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_y, fg_s_41, fg_s_42, fg_s_43, fi_s_78, fi_s_79, \
                         fi_s_80, fg_41, fg_42, fg_43, fh_58, fh_59, \
                         fh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -f_11 * fg_s_41[k]
                  + f_2 * fi_s_78[k]
                  + f_9 * fg_41[k]
                  + pb_y[k] * fh_58[k];

        t_79[k] = -f_8 * fg_s_42[k]
                  + f_2 * fi_s_79[k]
                  + f_0 * fg_42[k]
                  + pb_y[k] * fh_59[k];

        t_80[k] = -f_6 * fg_s_43[k]
                  + f_2 * fi_s_80[k]
                  + f_7 * fg_43[k]
                  + pb_y[k] * fh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, pb_y, pi_s_83, pi_83, di_83, fg_s_44, \
                         fi_s_81, fi_s_82, fi_s_83, fg_44, fh_61, \
                         fh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_4 * fg_s_44[k]
                  + f_2 * fi_s_81[k]
                  + f_5 * fg_44[k]
                  + pb_y[k] * fh_61[k];

        t_82[k] = f_2 * fi_s_82[k]
                  + pb_y[k] * fh_62[k];

        t_83[k] = -f_10 * pi_s_83[k]
                  + f_5 * pi_83[k]
                  + pa_x[k] * di_83[k]
                  + f_2 * fi_s_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_x, pb_y, pb_z, dh_21, dh_63, di_84, fi_s_84, \
                         fi_s_85, fi_s_86, fh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * dh_63[k]
                  + pa_x[k] * di_84[k]
                  + f_2 * fi_s_84[k];

        t_85[k] = f_7 * dh_21[k]
                  + f_2 * fi_s_85[k]
                  + pb_y[k] * fh_63[k];

        t_86[k] = f_2 * fi_s_86[k]
                  + pb_z[k] * fh_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_x, pb_z, dh_66, dh_68, di_87, di_89, fi_s_87, \
                         fi_s_88, fi_s_89, fh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_9 * dh_66[k]
                  + pa_x[k] * di_87[k]
                  + f_2 * fi_s_87[k];

        t_88[k] = f_2 * fi_s_88[k]
                  + pb_z[k] * fh_64[k];

        t_89[k] = f_9 * dh_68[k]
                  + pa_x[k] * di_89[k]
                  + f_2 * fi_s_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_x, pb_y, pb_z, dh_26, dh_69, di_90, fi_s_90, \
                         fi_s_91, fi_s_92, fh_66, fh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * dh_69[k]
                  + pa_x[k] * di_90[k]
                  + f_2 * fi_s_90[k];

        t_91[k] = f_2 * fi_s_91[k]
                  + pb_z[k] * fh_66[k];

        t_92[k] = f_7 * dh_26[k]
                  + f_2 * fi_s_92[k]
                  + pb_y[k] * fh_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_x, pb_z, dh_72, dh_73, di_93, di_94, fi_s_93, \
                         fi_s_94, fi_s_95, fh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * dh_72[k]
                  + pa_x[k] * di_93[k]
                  + f_2 * fi_s_93[k];

        t_94[k] = f_7 * dh_73[k]
                  + pa_x[k] * di_94[k]
                  + f_2 * fi_s_94[k];

        t_95[k] = f_2 * fi_s_95[k]
                  + pb_z[k] * fh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_x, pb_y, dh_30, dh_75, dh_77, di_96, di_98, \
                         fi_s_96, fi_s_97, fi_s_98, fh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_7 * dh_75[k]
                  + pa_x[k] * di_96[k]
                  + f_2 * fi_s_96[k];

        t_97[k] = f_7 * dh_30[k]
                  + f_2 * fi_s_97[k]
                  + pb_y[k] * fh_72[k];

        t_98[k] = f_7 * dh_77[k]
                  + pa_x[k] * di_98[k]
                  + f_2 * fi_s_98[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_z, dh_78, dh_80, fi_s_99, fi_s_100, \
                         fi_s_101, fh_73, fh_78, fh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_5 * dh_78[k]
                  + f_2 * fi_s_99[k]
                  + pb_x[k] * fh_78[k];

        t_100[k] = f_2 * fi_s_100[k]
                   + pb_z[k] * fh_73[k];

        t_101[k] = f_5 * dh_80[k]
                   + f_2 * fi_s_101[k]
                   + pb_x[k] * fh_80[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_x, dh_81, dh_82, dh_83, fi_s_102, fi_s_103, \
                         fi_s_104, fh_81, fh_82, fh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_5 * dh_81[k]
                   + f_2 * fi_s_102[k]
                   + pb_x[k] * fh_81[k];

        t_103[k] = f_5 * dh_82[k]
                   + f_2 * fi_s_103[k]
                   + pb_x[k] * fh_82[k];

        t_104[k] = f_5 * dh_83[k]
                   + f_2 * fi_s_104[k]
                   + pb_x[k] * fh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_z, di_105, di_107, di_108, \
                         fi_s_105, fi_s_106, fi_s_107, fi_s_108, \
                         fh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_x[k] * di_105[k]
                   + f_2 * fi_s_105[k];

        t_106[k] = f_2 * fi_s_106[k]
                   + pb_z[k] * fh_78[k];

        t_107[k] = pa_x[k] * di_107[k]
                   + f_2 * fi_s_107[k];

        t_108[k] = pa_x[k] * di_108[k]
                   + f_2 * fi_s_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_y, di_56, di_109, di_110, \
                         di_111, fi_s_109, fi_s_110, fi_s_111, \
                         fi_s_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * di_109[k]
                   + f_2 * fi_s_109[k];

        t_110[k] = pa_x[k] * di_110[k]
                   + f_2 * fi_s_110[k];

        t_111[k] = pa_x[k] * di_111[k]
                   + f_2 * fi_s_111[k];

        t_112[k] = pa_y[k] * di_56[k]
                   + f_2 * fi_s_112[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pa_z, pb_y, dh_44, di_29, di_31, \
                         di_58, fi_s_113, fi_s_114, fi_s_115, fi_s_116, \
                         fh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * di_29[k]
                   + f_2 * fi_s_113[k];

        t_114[k] = pa_y[k] * di_58[k]
                   + f_2 * fi_s_114[k];

        t_115[k] = pa_z[k] * di_31[k]
                   + f_2 * fi_s_115[k];

        t_116[k] = f_5 * dh_44[k]
                   + f_2 * fi_s_116[k]
                   + pb_y[k] * fh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_y, pa_z, pb_z, dh_24, di_34, di_61, fi_s_117, \
                         fi_s_118, fi_s_119, fh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * di_61[k]
                   + f_2 * fi_s_117[k];

        t_118[k] = pa_z[k] * di_34[k]
                   + f_2 * fi_s_118[k];

        t_119[k] = f_5 * dh_24[k]
                   + f_2 * fi_s_119[k]
                   + pb_z[k] * fh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_y, pa_z, pb_y, dh_47, di_38, di_65, fi_s_120, \
                         fi_s_121, fi_s_122, fh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_5 * dh_47[k]
                   + f_2 * fi_s_120[k]
                   + pb_y[k] * fh_89[k];

        t_121[k] = pa_y[k] * di_65[k]
                   + f_2 * fi_s_121[k];

        t_122[k] = pa_z[k] * di_38[k]
                   + f_2 * fi_s_122[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pa_x, pb_y, pb_z, dh_27, dh_51, dh_96, di_124, \
                         fi_s_123, fi_s_124, fi_s_125, fh_90, fh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_5 * dh_27[k]
                   + f_2 * fi_s_123[k]
                   + pb_z[k] * fh_90[k];

        t_124[k] = f_7 * dh_96[k]
                   + pa_x[k] * di_124[k]
                   + f_2 * fi_s_124[k];

        t_125[k] = f_5 * dh_51[k]
                   + f_2 * fi_s_125[k]
                   + pb_y[k] * fh_93[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_y, pa_z, pb_x, dh_100, di_43, di_70, \
                         fi_s_126, fi_s_127, fi_s_128, fh_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_y[k] * di_70[k]
                   + f_2 * fi_s_126[k];

        t_127[k] = pa_z[k] * di_43[k]
                   + f_2 * fi_s_127[k];

        t_128[k] = f_5 * dh_100[k]
                   + f_2 * fi_s_128[k]
                   + pb_x[k] * fh_100[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pb_x, dh_101, dh_102, dh_103, fi_s_129, \
                         fi_s_130, fi_s_131, fh_101, fh_102, fh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_5 * dh_101[k]
                   + f_2 * fi_s_129[k]
                   + pb_x[k] * fh_101[k];

        t_130[k] = f_5 * dh_102[k]
                   + f_2 * fi_s_130[k]
                   + pb_x[k] * fh_102[k];

        t_131[k] = f_5 * dh_103[k]
                   + f_2 * fi_s_131[k]
                   + pb_x[k] * fh_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pa_y, di_76, di_133, di_134, \
                         di_135, fi_s_132, fi_s_133, fi_s_134, \
                         fi_s_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * di_76[k]
                   + f_2 * fi_s_132[k];

        t_133[k] = pa_x[k] * di_133[k]
                   + f_2 * fi_s_133[k];

        t_134[k] = pa_x[k] * di_134[k]
                   + f_2 * fi_s_134[k];

        t_135[k] = pa_x[k] * di_135[k]
                   + f_2 * fi_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, di_136, di_137, di_138, di_139, \
                         fi_s_136, fi_s_137, fi_s_138, fi_s_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * di_136[k]
                   + f_2 * fi_s_136[k];

        t_137[k] = pa_x[k] * di_137[k]
                   + f_2 * fi_s_137[k];

        t_138[k] = pa_x[k] * di_138[k]
                   + f_2 * fi_s_138[k];

        t_139[k] = pa_x[k] * di_139[k]
                   + f_2 * fi_s_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_x, pb_y, pb_z, dh_42, dh_105, di_140, \
                         fi_s_140, fi_s_141, fi_s_142, fh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_12 * dh_105[k]
                   + pa_x[k] * di_140[k]
                   + f_2 * fi_s_140[k];

        t_141[k] = f_2 * fi_s_141[k]
                   + pb_y[k] * fh_105[k];

        t_142[k] = f_7 * dh_42[k]
                   + f_2 * fi_s_142[k]
                   + pb_z[k] * fh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pb_y, dh_108, dh_110, di_143, di_145, \
                         fi_s_143, fi_s_144, fi_s_145, fh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_9 * dh_108[k]
                   + pa_x[k] * di_143[k]
                   + f_2 * fi_s_143[k];

        t_144[k] = f_2 * fi_s_144[k]
                   + pb_y[k] * fh_107[k];

        t_145[k] = f_9 * dh_110[k]
                   + pa_x[k] * di_145[k]
                   + f_2 * fi_s_145[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_x, pb_y, dh_111, dh_112, di_146, di_147, \
                         fi_s_146, fi_s_147, fi_s_148, fh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * dh_111[k]
                   + pa_x[k] * di_146[k]
                   + f_2 * fi_s_146[k];

        t_147[k] = f_0 * dh_112[k]
                   + pa_x[k] * di_147[k]
                   + f_2 * fi_s_147[k];

        t_148[k] = f_2 * fi_s_148[k]
                   + pb_y[k] * fh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_x, dh_114, dh_115, dh_116, di_149, di_150, \
                         di_151, fi_s_149, fi_s_150, fi_s_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * dh_114[k]
                   + pa_x[k] * di_149[k]
                   + f_2 * fi_s_149[k];

        t_150[k] = f_7 * dh_115[k]
                   + pa_x[k] * di_150[k]
                   + f_2 * fi_s_150[k];

        t_151[k] = f_7 * dh_116[k]
                   + pa_x[k] * di_151[k]
                   + f_2 * fi_s_151[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pa_x, pb_y, dh_117, dh_119, di_152, di_154, \
                         fi_s_152, fi_s_153, fi_s_154, fh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_7 * dh_117[k]
                   + pa_x[k] * di_152[k]
                   + f_2 * fi_s_152[k];

        t_153[k] = f_2 * fi_s_153[k]
                   + pb_y[k] * fh_114[k];

        t_154[k] = f_7 * dh_119[k]
                   + pa_x[k] * di_154[k]
                   + f_2 * fi_s_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pb_x, dh_120, dh_121, dh_122, fi_s_155, \
                         fi_s_156, fi_s_157, fh_120, fh_121, fh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_5 * dh_120[k]
                   + f_2 * fi_s_155[k]
                   + pb_x[k] * fh_120[k];

        t_156[k] = f_5 * dh_121[k]
                   + f_2 * fi_s_156[k]
                   + pb_x[k] * fh_121[k];

        t_157[k] = f_5 * dh_122[k]
                   + f_2 * fi_s_157[k]
                   + pb_x[k] * fh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pb_x, pb_y, dh_123, dh_125, fi_s_158, fi_s_159, \
                         fi_s_160, fh_119, fh_123, fh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_5 * dh_123[k]
                   + f_2 * fi_s_158[k]
                   + pb_x[k] * fh_123[k];

        t_159[k] = f_2 * fi_s_159[k]
                   + pb_y[k] * fh_119[k];

        t_160[k] = f_5 * dh_125[k]
                   + f_2 * fi_s_160[k]
                   + pb_x[k] * fh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, pa_x, di_161, di_162, di_163, \
                         di_164, di_165, fi_s_161, fi_s_162, fi_s_163, fi_s_164, \
                         fi_s_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pa_x[k] * di_161[k]
                   + f_2 * fi_s_161[k];

        t_162[k] = pa_x[k] * di_162[k]
                   + f_2 * fi_s_162[k];

        t_163[k] = pa_x[k] * di_163[k]
                   + f_2 * fi_s_163[k];

        t_164[k] = pa_x[k] * di_164[k]
                   + f_2 * fi_s_164[k];

        t_165[k] = pa_x[k] * di_165[k]
                   + f_2 * fi_s_165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_x, pb_x, pb_y, di_167, fg_s_90, fi_s_166, \
                         fi_s_167, fi_s_168, fg_90, fh_125, fh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_2 * fi_s_166[k]
                   + pb_y[k] * fh_125[k];

        t_167[k] = pa_x[k] * di_167[k]
                   + f_2 * fi_s_167[k];

        t_168[k] = -f_1 * fg_s_90[k]
                   + f_2 * fi_s_168[k]
                   + f_3 * fg_90[k]
                   + pb_x[k] * fh_126[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pb_x, pb_z, fg_s_91, fg_s_93, fi_s_169, \
                         fi_s_170, fi_s_171, fg_91, fg_93, fh_126, fh_127, \
                         fh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -f_11 * fg_s_91[k]
                   + f_2 * fi_s_169[k]
                   + f_9 * fg_91[k]
                   + pb_x[k] * fh_127[k];

        t_170[k] = f_2 * fi_s_170[k]
                   + pb_z[k] * fh_126[k];

        t_171[k] = -f_8 * fg_s_93[k]
                   + f_2 * fi_s_171[k]
                   + f_0 * fg_93[k]
                   + pb_x[k] * fh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, fg_s_95, fg_s_96, fi_s_172, \
                         fi_s_173, fi_s_174, fg_95, fg_96, fh_127, fh_131, \
                         fh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_2 * fi_s_172[k]
                   + pb_z[k] * fh_127[k];

        t_173[k] = -f_8 * fg_s_95[k]
                   + f_2 * fi_s_173[k]
                   + f_0 * fg_95[k]
                   + pb_x[k] * fh_131[k];

        t_174[k] = -f_6 * fg_s_96[k]
                   + f_2 * fi_s_174[k]
                   + f_7 * fg_96[k]
                   + pb_x[k] * fh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, pb_z, fg_s_98, fg_s_99, fi_s_175, \
                         fi_s_176, fi_s_177, fg_98, fg_99, fh_129, fh_134, \
                         fh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_2 * fi_s_175[k]
                   + pb_z[k] * fh_129[k];

        t_176[k] = -f_6 * fg_s_98[k]
                   + f_2 * fi_s_176[k]
                   + f_7 * fg_98[k]
                   + pb_x[k] * fh_134[k];

        t_177[k] = -f_6 * fg_s_99[k]
                   + f_2 * fi_s_177[k]
                   + f_7 * fg_99[k]
                   + pb_x[k] * fh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pb_x, pb_z, fg_s_100, fg_s_102, fi_s_178, \
                         fi_s_179, fi_s_180, fg_100, fg_102, fh_132, fh_136, \
                         fh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = -f_4 * fg_s_100[k]
                   + f_2 * fi_s_178[k]
                   + f_5 * fg_100[k]
                   + pb_x[k] * fh_136[k];

        t_179[k] = f_2 * fi_s_179[k]
                   + pb_z[k] * fh_132[k];

        t_180[k] = -f_4 * fg_s_102[k]
                   + f_2 * fi_s_180[k]
                   + f_5 * fg_102[k]
                   + pb_x[k] * fh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, pb_x, fg_s_103, fg_s_104, fi_s_181, fi_s_182, \
                         fi_s_183, fg_103, fg_104, fh_139, fh_140, \
                         fh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = -f_4 * fg_s_103[k]
                   + f_2 * fi_s_181[k]
                   + f_5 * fg_103[k]
                   + pb_x[k] * fh_139[k];

        t_182[k] = -f_4 * fg_s_104[k]
                   + f_2 * fi_s_182[k]
                   + f_5 * fg_104[k]
                   + pb_x[k] * fh_140[k];

        t_183[k] = f_2 * fi_s_183[k]
                   + pb_x[k] * fh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pb_x, fi_s_184, fi_s_185, \
                         fi_s_186, fi_s_187, fi_s_188, fh_142, fh_143, fh_144, fh_145, \
                         fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_2 * fi_s_184[k]
                   + pb_x[k] * fh_142[k];

        t_185[k] = f_2 * fi_s_185[k]
                   + pb_x[k] * fh_143[k];

        t_186[k] = f_2 * fi_s_186[k]
                   + pb_x[k] * fh_144[k];

        t_187[k] = f_2 * fi_s_187[k]
                   + pb_x[k] * fh_145[k];

        t_188[k] = f_2 * fi_s_188[k]
                   + pb_x[k] * fh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, dh_78, fg_s_100, fi_s_189, fi_s_190, \
                         fi_s_191, fg_100, fh_141, fh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * dh_78[k]
                   - f_1 * fg_s_100[k]
                   + f_2 * fi_s_189[k]
                   + f_3 * fg_100[k]
                   + pb_y[k] * fh_141[k];

        t_190[k] = f_2 * fi_s_190[k]
                   + pb_z[k] * fh_141[k];

        t_191[k] = -f_4 * fg_s_100[k]
                   + f_2 * fi_s_191[k]
                   + f_5 * fg_100[k]
                   + pb_z[k] * fh_142[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pb_z, dh_83, fg_s_101, fg_s_102, fi_s_192, \
                         fi_s_193, fi_s_194, fg_101, fg_102, fh_143, fh_144, \
                         fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -f_6 * fg_s_101[k]
                   + f_2 * fi_s_192[k]
                   + f_7 * fg_101[k]
                   + pb_z[k] * fh_143[k];

        t_193[k] = -f_8 * fg_s_102[k]
                   + f_2 * fi_s_193[k]
                   + f_0 * fg_102[k]
                   + pb_z[k] * fh_144[k];

        t_194[k] = f_0 * dh_83[k]
                   + f_2 * fi_s_194[k]
                   + pb_y[k] * fh_146[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_z, pb_z, di_84, di_85, fg_s_104, fi_s_195, \
                         fi_s_196, fi_s_197, fg_104, fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -f_1 * fg_s_104[k]
                   + f_2 * fi_s_195[k]
                   + f_3 * fg_104[k]
                   + pb_z[k] * fh_146[k];

        t_196[k] = pa_z[k] * di_84[k]
                   + f_2 * fi_s_196[k];

        t_197[k] = pa_z[k] * di_85[k]
                   + f_2 * fi_s_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_z, pb_x, dh_64, di_87, di_88, fg_s_107, \
                         fi_s_198, fi_s_199, fi_s_200, fg_107, fh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = -f_11 * fg_s_107[k]
                   + f_2 * fi_s_198[k]
                   + f_9 * fg_107[k]
                   + pb_x[k] * fh_149[k];

        t_199[k] = pa_z[k] * di_87[k]
                   + f_2 * fi_s_199[k];

        t_200[k] = f_5 * dh_64[k]
                   + pa_z[k] * di_88[k]
                   + f_2 * fi_s_200[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pb_x, dh_66, di_90, di_91, fg_s_110, \
                         fi_s_201, fi_s_202, fi_s_203, fg_110, fh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -f_8 * fg_s_110[k]
                   + f_2 * fi_s_201[k]
                   + f_0 * fg_110[k]
                   + pb_x[k] * fh_152[k];

        t_202[k] = pa_z[k] * di_90[k]
                   + f_2 * fi_s_202[k];

        t_203[k] = f_5 * dh_66[k]
                   + pa_z[k] * di_91[k]
                   + f_2 * fi_s_203[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pb_x, dh_67, di_92, di_94, fg_s_114, \
                         fi_s_204, fi_s_205, fi_s_206, fg_114, fh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_7 * dh_67[k]
                   + pa_z[k] * di_92[k]
                   + f_2 * fi_s_204[k];

        t_205[k] = -f_6 * fg_s_114[k]
                   + f_2 * fi_s_205[k]
                   + f_7 * fg_114[k]
                   + pb_x[k] * fh_156[k];

        t_206[k] = pa_z[k] * di_94[k]
                   + f_2 * fi_s_206[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, dh_69, dh_70, dh_71, di_95, di_96, di_97, \
                         fi_s_207, fi_s_208, fi_s_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_5 * dh_69[k]
                   + pa_z[k] * di_95[k]
                   + f_2 * fi_s_207[k];

        t_208[k] = f_7 * dh_70[k]
                   + pa_z[k] * di_96[k]
                   + f_2 * fi_s_208[k];

        t_209[k] = f_0 * dh_71[k]
                   + pa_z[k] * di_97[k]
                   + f_2 * fi_s_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pb_x, fg_s_119, fi_s_210, fi_s_211, \
                         fi_s_212, fi_s_213, fg_119, fh_161, fh_162, fh_163, \
                         fh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -f_4 * fg_s_119[k]
                   + f_2 * fi_s_210[k]
                   + f_5 * fg_119[k]
                   + pb_x[k] * fh_161[k];

        t_211[k] = f_2 * fi_s_211[k]
                   + pb_x[k] * fh_162[k];

        t_212[k] = f_2 * fi_s_212[k]
                   + pb_x[k] * fh_163[k];

        t_213[k] = f_2 * fi_s_213[k]
                   + pb_x[k] * fh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pb_x, di_105, fi_s_214, fi_s_215, \
                         fi_s_216, fi_s_217, fh_165, fh_166, fh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_2 * fi_s_214[k]
                   + pb_x[k] * fh_165[k];

        t_215[k] = f_2 * fi_s_215[k]
                   + pb_x[k] * fh_166[k];

        t_216[k] = f_2 * fi_s_216[k]
                   + pb_x[k] * fh_167[k];

        t_217[k] = pa_z[k] * di_105[k]
                   + f_2 * fi_s_217[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pa_z, pb_z, dh_78, dh_79, dh_80, di_107, di_108, \
                         fi_s_218, fi_s_219, fi_s_220, fh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_5 * dh_78[k]
                   + f_2 * fi_s_218[k]
                   + pb_z[k] * fh_162[k];

        t_219[k] = f_7 * dh_79[k]
                   + pa_z[k] * di_107[k]
                   + f_2 * fi_s_219[k];

        t_220[k] = f_0 * dh_80[k]
                   + pa_z[k] * di_108[k]
                   + f_2 * fi_s_220[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pa_y, pa_z, pb_y, pi_s_83, pi_83, dh_81, dh_104, \
                         di_109, di_139, fi_s_221, fi_s_222, fi_s_223, \
                         fh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_9 * dh_81[k]
                   + pa_z[k] * di_109[k]
                   + f_2 * fi_s_221[k];

        t_222[k] = f_7 * dh_104[k]
                   + f_2 * fi_s_222[k]
                   + pb_y[k] * fh_167[k];

        t_223[k] = -f_10 * pi_s_83[k]
                   + f_5 * pi_83[k]
                   + pa_y[k] * di_139[k]
                   + f_2 * fi_s_223[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_y, dh_105, dh_106, di_140, di_141, \
                         di_142, di_143, fi_s_224, fi_s_225, fi_s_226, \
                         fi_s_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * di_140[k]
                   + f_2 * fi_s_224[k];

        t_225[k] = f_5 * dh_105[k]
                   + pa_y[k] * di_141[k]
                   + f_2 * fi_s_225[k];

        t_226[k] = pa_y[k] * di_142[k]
                   + f_2 * fi_s_226[k];

        t_227[k] = f_7 * dh_106[k]
                   + pa_y[k] * di_143[k]
                   + f_2 * fi_s_227[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, dh_107, dh_108, dh_109, di_144, \
                         di_145, di_146, di_147, fi_s_228, fi_s_229, fi_s_230, \
                         fi_s_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * dh_107[k]
                   + pa_y[k] * di_144[k]
                   + f_2 * fi_s_228[k];

        t_229[k] = pa_y[k] * di_145[k]
                   + f_2 * fi_s_229[k];

        t_230[k] = f_0 * dh_108[k]
                   + pa_y[k] * di_146[k]
                   + f_2 * fi_s_230[k];

        t_231[k] = f_7 * dh_109[k]
                   + pa_y[k] * di_147[k]
                   + f_2 * fi_s_231[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_y, dh_110, dh_111, dh_112, di_148, \
                         di_149, di_150, di_151, fi_s_232, fi_s_233, fi_s_234, \
                         fi_s_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_5 * dh_110[k]
                   + pa_y[k] * di_148[k]
                   + f_2 * fi_s_232[k];

        t_233[k] = pa_y[k] * di_149[k]
                   + f_2 * fi_s_233[k];

        t_234[k] = f_9 * dh_111[k]
                   + pa_y[k] * di_150[k]
                   + f_2 * fi_s_234[k];

        t_235[k] = f_0 * dh_112[k]
                   + pa_y[k] * di_151[k]
                   + f_2 * fi_s_235[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pb_x, dh_113, dh_114, di_152, \
                         di_153, di_154, fi_s_236, fi_s_237, fi_s_238, fi_s_239, \
                         fh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_7 * dh_113[k]
                   + pa_y[k] * di_152[k]
                   + f_2 * fi_s_236[k];

        t_237[k] = f_5 * dh_114[k]
                   + pa_y[k] * di_153[k]
                   + f_2 * fi_s_237[k];

        t_238[k] = pa_y[k] * di_154[k]
                   + f_2 * fi_s_238[k];

        t_239[k] = f_2 * fi_s_239[k]
                   + pb_x[k] * fh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pb_x, fi_s_240, fi_s_241, \
                         fi_s_242, fi_s_243, fi_s_244, fh_184, fh_185, fh_186, fh_187, \
                         fh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_2 * fi_s_240[k]
                   + pb_x[k] * fh_184[k];

        t_241[k] = f_2 * fi_s_241[k]
                   + pb_x[k] * fh_185[k];

        t_242[k] = f_2 * fi_s_242[k]
                   + pb_x[k] * fh_186[k];

        t_243[k] = f_2 * fi_s_243[k]
                   + pb_x[k] * fh_187[k];

        t_244[k] = f_2 * fi_s_244[k]
                   + pb_x[k] * fh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pa_y, pb_z, dh_99, dh_120, dh_122, di_161, \
                         di_163, fi_s_245, fi_s_246, fi_s_247, fh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_12 * dh_120[k]
                   + pa_y[k] * di_161[k]
                   + f_2 * fi_s_245[k];

        t_246[k] = f_7 * dh_99[k]
                   + f_2 * fi_s_246[k]
                   + pb_z[k] * fh_183[k];

        t_247[k] = f_9 * dh_122[k]
                   + pa_y[k] * di_163[k]
                   + f_2 * fi_s_247[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pa_y, pb_y, dh_123, dh_124, dh_125, di_164, \
                         di_165, fi_s_248, fi_s_249, fi_s_250, fh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_0 * dh_123[k]
                   + pa_y[k] * di_164[k]
                   + f_2 * fi_s_248[k];

        t_249[k] = f_7 * dh_124[k]
                   + pa_y[k] * di_165[k]
                   + f_2 * fi_s_249[k];

        t_250[k] = f_5 * dh_125[k]
                   + f_2 * fi_s_250[k]
                   + pb_y[k] * fh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, pa_y, pb_x, pb_y, di_167, fg_s_135, fi_s_251, \
                         fi_s_252, fi_s_253, fg_135, fh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * di_167[k]
                   + f_2 * fi_s_251[k];

        t_252[k] = -f_1 * fg_s_135[k]
                   + f_2 * fi_s_252[k]
                   + f_3 * fg_135[k]
                   + pb_x[k] * fh_189[k];

        t_253[k] = f_2 * fi_s_253[k]
                   + pb_y[k] * fh_189[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_x, pb_y, fg_s_137, fg_s_138, fi_s_254, \
                         fi_s_255, fi_s_256, fg_137, fg_138, fh_191, \
                         fh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = -f_11 * fg_s_137[k]
                   + f_2 * fi_s_254[k]
                   + f_9 * fg_137[k]
                   + pb_x[k] * fh_191[k];

        t_255[k] = -f_8 * fg_s_138[k]
                   + f_2 * fi_s_255[k]
                   + f_0 * fg_138[k]
                   + pb_x[k] * fh_192[k];

        t_256[k] = f_2 * fi_s_256[k]
                   + pb_y[k] * fh_191[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_x, fg_s_140, fg_s_141, fg_s_142, fi_s_257, \
                         fi_s_258, fi_s_259, fg_140, fg_141, fg_142, fh_194, fh_195, \
                         fh_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -f_8 * fg_s_140[k]
                   + f_2 * fi_s_257[k]
                   + f_0 * fg_140[k]
                   + pb_x[k] * fh_194[k];

        t_258[k] = -f_6 * fg_s_141[k]
                   + f_2 * fi_s_258[k]
                   + f_7 * fg_141[k]
                   + pb_x[k] * fh_195[k];

        t_259[k] = -f_6 * fg_s_142[k]
                   + f_2 * fi_s_259[k]
                   + f_7 * fg_142[k]
                   + pb_x[k] * fh_196[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_x, pb_y, fg_s_144, fg_s_145, fi_s_260, \
                         fi_s_261, fi_s_262, fg_144, fg_145, fh_194, fh_198, \
                         fh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_2 * fi_s_260[k]
                   + pb_y[k] * fh_194[k];

        t_261[k] = -f_6 * fg_s_144[k]
                   + f_2 * fi_s_261[k]
                   + f_7 * fg_144[k]
                   + pb_x[k] * fh_198[k];

        t_262[k] = -f_4 * fg_s_145[k]
                   + f_2 * fi_s_262[k]
                   + f_5 * fg_145[k]
                   + pb_x[k] * fh_199[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_x, pb_y, fg_s_146, fg_s_147, fi_s_263, \
                         fi_s_264, fi_s_265, fg_146, fg_147, fh_198, fh_200, \
                         fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -f_4 * fg_s_146[k]
                   + f_2 * fi_s_263[k]
                   + f_5 * fg_146[k]
                   + pb_x[k] * fh_200[k];

        t_264[k] = -f_4 * fg_s_147[k]
                   + f_2 * fi_s_264[k]
                   + f_5 * fg_147[k]
                   + pb_x[k] * fh_201[k];

        t_265[k] = f_2 * fi_s_265[k]
                   + pb_y[k] * fh_198[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pb_x, fg_s_149, fi_s_266, fi_s_267, \
                         fi_s_268, fi_s_269, fg_149, fh_203, fh_204, fh_205, \
                         fh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -f_4 * fg_s_149[k]
                   + f_2 * fi_s_266[k]
                   + f_5 * fg_149[k]
                   + pb_x[k] * fh_203[k];

        t_267[k] = f_2 * fi_s_267[k]
                   + pb_x[k] * fh_204[k];

        t_268[k] = f_2 * fi_s_268[k]
                   + pb_x[k] * fh_205[k];

        t_269[k] = f_2 * fi_s_269[k]
                   + pb_x[k] * fh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_x, pb_y, fg_s_145, fi_s_270, fi_s_271, \
                         fi_s_272, fi_s_273, fg_145, fh_204, fh_207, fh_208, \
                         fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_2 * fi_s_270[k]
                   + pb_x[k] * fh_207[k];

        t_271[k] = f_2 * fi_s_271[k]
                   + pb_x[k] * fh_208[k];

        t_272[k] = f_2 * fi_s_272[k]
                   + pb_x[k] * fh_209[k];

        t_273[k] = -f_1 * fg_s_145[k]
                   + f_2 * fi_s_273[k]
                   + f_3 * fg_145[k]
                   + pb_y[k] * fh_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pb_y, fg_s_146, fg_s_147, fg_s_148, fi_s_274, \
                         fi_s_275, fi_s_276, fg_146, fg_147, fg_148, fh_205, fh_206, \
                         fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = -f_11 * fg_s_146[k]
                   + f_2 * fi_s_274[k]
                   + f_9 * fg_146[k]
                   + pb_y[k] * fh_205[k];

        t_275[k] = -f_8 * fg_s_147[k]
                   + f_2 * fi_s_275[k]
                   + f_0 * fg_147[k]
                   + pb_y[k] * fh_206[k];

        t_276[k] = -f_6 * fg_s_148[k]
                   + f_2 * fi_s_276[k]
                   + f_7 * fg_148[k]
                   + pb_y[k] * fh_207[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pb_y, pb_z, dh_125, fg_s_149, fi_s_277, \
                         fi_s_278, fi_s_279, fg_149, fh_208, fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -f_4 * fg_s_149[k]
                   + f_2 * fi_s_277[k]
                   + f_5 * fg_149[k]
                   + pb_y[k] * fh_208[k];

        t_278[k] = f_2 * fi_s_278[k]
                   + pb_y[k] * fh_209[k];

        t_279[k] = f_0 * dh_125[k]
                   - f_1 * fg_s_149[k]
                   + f_2 * fi_s_279[k]
                   + f_3 * fg_149[k]
                   + pb_z[k] * fh_209[k];
    }
}

}  // namespace simdkin
